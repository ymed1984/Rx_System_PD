# OMA-to-BER Dynamic Effects ExecPlan

## Purpose

Plan extensions from the current scalar NRZ/OOK Gaussian OMA-to-BER MVP toward dynamic receiver effects:

- Bandwidth penalty.
- Time-domain waveform simulation.
- ISI estimation.
- S-parameter import and convolution.
- Photodiode saturation and power handling.

The plan preserves the current tested scalar model as the baseline. Dynamic features must be added incrementally and validated against simple limiting cases before being used for design conclusions.

## Current Baseline

The existing MVP computes:

```text
OMA, ER
  -> P0, P1, Pavg
  -> I0, I1
  -> level-dependent noise
  -> optimum threshold
  -> Q, BER
```

The baseline is deterministic, scalar, and Gaussian. It does not model waveform memory, ISI, PD/TIA frequency response, or saturation.

Dynamic extensions should not replace this baseline. They should add optional paths that either:

1. Convert dynamic penalties into an effective scalar penalty for the existing BER flow.
2. Run a separate waveform-based analysis with explicit assumptions.

## Design Principles

- Keep units explicit in names and docstrings.
- Keep core numerical functions deterministic and side-effect free.
- Put plotting in examples or `plotting.py`.
- Do not claim standards compliance.
- Do not implement PAM4 until NRZ/OOK dynamic behavior is tested.
- Add small APIs with clear validation before adding larger workflows.
- Prefer proven numerical primitives from `numpy` and `scipy`.

## Recommended Implementation Order

1. Simple bandwidth penalty model.
2. Deterministic NRZ/OOK waveform generation.
3. Linear time-invariant filtering and ISI metrics.
4. S-parameter import and impulse response conversion.
5. Saturation and responsivity compression.
6. Combined dynamic OMA-to-BER workflow.

This order is intentional. Bandwidth penalty can be validated analytically, waveform generation can be tested without analog modeling, filtering can be tested with simple low-pass responses, and only then should measured S-parameters be imported.

## Milestone D1: Simple Bandwidth Penalty Model

### Goal

Add a first-order way to estimate OMA penalty from limited receiver bandwidth without time-domain simulation.

### Proposed files

```text
src/oma_ber/bandwidth.py
tests/test_bandwidth.py
examples/06_bandwidth_penalty.py
```

### Proposed APIs

```python
def first_order_lowpass_mag(
    frequency_hz: float | np.ndarray,
    bandwidth_3db_hz: float,
) -> float | np.ndarray: ...


def bandwidth_penalty_db(
    signal_frequency_hz: float,
    bandwidth_3db_hz: float,
) -> float: ...


def apply_oma_penalty_db(
    oma_dbm: float,
    penalty_db: float,
) -> float: ...
```

### Model

For a first-order low-pass response:

```text
|H(f)| = 1 / sqrt(1 + (f / f3dB)^2)
penalty_dB = -20 * log10(|H(f_signal)|)
OMA_after_penalty_dBm = OMA_before_penalty_dBm - penalty_dB
```

The signal frequency is user-supplied. Do not silently infer it from baud rate in this first pass.

### Validation

- `bandwidth_3db_hz > 0`.
- `signal_frequency_hz >= 0`.
- `penalty_db >= 0`.
- At `f = 0`, penalty is 0 dB.
- At `f = f3dB`, magnitude is approximately `1 / sqrt(2)` and penalty is approximately 3.01 dB.

### Why first

This gives immediate system-level value and keeps the calculation compatible with the existing scalar BER path.

## Milestone D2: NRZ/OOK Waveform Generation

### Goal

Create deterministic optical or current waveforms from bit patterns for later filtering and ISI analysis.

### Proposed files

```text
src/oma_ber/waveform.py
tests/test_waveform.py
examples/07_nrz_waveform.py
```

### Proposed APIs

```python
def samples_per_symbol(sample_rate_hz: float, symbol_rate_baud: float) -> int: ...


def nrz_bits_to_levels(
    bits: np.ndarray,
    low_level: float,
    high_level: float,
    samples_per_symbol: int,
) -> np.ndarray: ...


def prbs_bits(
    num_bits: int,
    seed: int = 1,
) -> np.ndarray: ...
```

### Validation

- `sample_rate_hz > 0`.
- `symbol_rate_baud > 0`.
- `sample_rate_hz / symbol_rate_baud` must produce at least 2 samples/symbol.
- `bits` must contain only 0 and 1.
- `high_level > low_level`.

### Tests

- Known bit sequence maps to expected repeated levels.
- Output length equals `len(bits) * samples_per_symbol`.
- PRBS is deterministic for fixed seed.

### Notes

This milestone still does not add noise or BER. It only creates deterministic waveforms.

## Milestone D3: Linear Filtering and ISI Metrics

### Goal

Apply an LTI response to the waveform and estimate vertical eye closure or sampled-level degradation.

### Proposed files

```text
src/oma_ber/isi.py
tests/test_isi.py
examples/08_bandlimited_eye_metrics.py
```

### Proposed APIs

```python
def first_order_lowpass_impulse_response(
    sample_rate_hz: float,
    bandwidth_3db_hz: float,
    num_taps: int,
) -> np.ndarray: ...


def apply_lti_filter(
    waveform: np.ndarray,
    impulse_response: np.ndarray,
) -> np.ndarray: ...


def sample_at_symbol_centers(
    waveform: np.ndarray,
    samples_per_symbol: int,
    timing_offset_samples: int | None = None,
) -> np.ndarray: ...


def sampled_eye_levels(
    samples: np.ndarray,
    bits: np.ndarray,
) -> dict[str, float]: ...
```

### Metrics

Start with simple deterministic metrics:

```text
mean_zero_level
mean_one_level
min_one_level
max_zero_level
vertical_eye_opening = min_one_level - max_zero_level
isi_penalty_db = 20 * log10(ideal_delta / filtered_delta)
```

Use these as engineering indicators, not standards metrics.

### Validation

- Impulse response length positive.
- Waveform and impulse response are one-dimensional.
- `samples_per_symbol >= 2`.
- `bits` length matches sampled symbols.

### Tests

- Delta impulse leaves waveform unchanged.
- Low-pass impulse reduces high-frequency transitions.
- Eye opening decreases when bandwidth is reduced.

## Milestone D4: S-Parameter Import and Impulse Response

### Goal

Use measured or simulated frequency response to filter waveforms.

### Proposed files

```text
src/oma_ber/sparameters.py
tests/test_sparameters.py
examples/09_sparameter_filtering.py
```

### Dependency policy

Avoid adding a large dependency initially. Start with simple Touchstone `.s1p` / `.s2p` parsing only if needed, or support CSV frequency response first:

```text
frequency_hz, magnitude_db, phase_deg
```

Add `scikit-rf` only if the project clearly needs robust Touchstone support.

### Proposed APIs

```python
def frequency_response_from_csv(
    path: str | Path,
) -> tuple[np.ndarray, np.ndarray]: ...


def impulse_response_from_frequency_response(
    frequency_hz: np.ndarray,
    response_complex: np.ndarray,
    sample_rate_hz: float,
    num_taps: int,
) -> np.ndarray: ...
```

### Validation

- Frequencies must be positive and strictly increasing.
- Response length must match frequency length.
- Sample rate must exceed twice the maximum frequency used for interpolation.
- Interpolation behavior must be documented.

### Tests

- Flat response produces near-delta impulse.
- First-order synthetic response roughly matches analytic low-pass response.
- Invalid frequency ordering raises `ValueError`.

## Milestone D5: Saturation and Responsivity Compression

### Goal

Model photodiode power handling as a static nonlinear effect before waveform filtering.

### Proposed files

```text
src/oma_ber/saturation.py
tests/test_saturation.py
examples/10_saturation_sweep.py
```

### Proposed APIs

```python
def compressed_responsivity_a_per_w(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
    compression_order: float = 1.0,
) -> float | np.ndarray: ...


def saturated_photocurrent_a(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
    compression_order: float = 1.0,
) -> float | np.ndarray: ...
```

### Example model

A simple bounded model:

```text
I = R * P / (1 + (P / Psat)^order)
```

This is a simplified model and must be labeled as such.

### Validation

- `responsivity_a_per_w > 0`.
- `saturation_power_w > 0`.
- `compression_order > 0`.
- `optical_power_w >= 0`.

### Tests

- At low power, current is approximately `R * P`.
- At high power, incremental responsivity decreases.
- Current remains non-negative.

## Milestone D6: Combined Dynamic OMA-to-BER Workflow

### Goal

Connect dynamic penalties back into BER estimation without pretending it is a full standards eye simulation.

### Proposed files

```text
src/oma_ber/dynamic_receiver.py
tests/test_dynamic_receiver.py
examples/11_dynamic_oma_to_ber.py
```

### Proposed APIs

```python
def calculate_ber_with_bandwidth_penalty(
    oma_dbm: float,
    er_db: float,
    pd: Photodiode,
    rx: Receiver,
    signal_frequency_hz: float,
) -> dict: ...


def calculate_waveform_eye_metrics(
    bits: np.ndarray,
    oma_dbm: float,
    er_db: float,
    pd: Photodiode,
    sample_rate_hz: float,
    symbol_rate_baud: float,
    impulse_response: np.ndarray,
) -> dict: ...
```

### Output

Include explicit fields:

```text
oma_dbm
effective_oma_dbm
bandwidth_penalty_db
vertical_eye_opening_a
isi_penalty_db
ber
ber_from_q
```

### Tests

- Infinite or very high bandwidth gives near-zero penalty.
- Reduced bandwidth worsens effective OMA and BER.
- Deterministic bit patterns produce deterministic metrics.

## Data Model Extensions

Possible future dataclasses:

```python
@dataclass(frozen=True)
class WaveformConfig:
    sample_rate_hz: float
    symbol_rate_baud: float
    samples_per_symbol: int


@dataclass(frozen=True)
class SaturationModel:
    saturation_power_w: float
    compression_order: float = 1.0


@dataclass(frozen=True)
class FrequencyResponse:
    frequency_hz: np.ndarray
    response_complex: np.ndarray
```

Do not add these until the corresponding milestone needs them.

## Validation Checklist

Before finalizing each milestone:

- Units are explicit in names and docstrings.
- No dBm/W or dB/linear mixing.
- Scalar BER baseline tests still pass.
- Dynamic penalty worsens or preserves BER; it should not improve BER unless explicitly modeling equalization.
- Examples run from repository root.
- Plotting and file I/O stay out of core numerical modules.
- Simplified models are labeled as simplified.

## Required Checks

For each implementation milestone, run:

```bash
uv run pytest
uv run ruff check .
```

For example milestones, also run the new example script from the repository root.

## Recommended First Implementation Slice

Start with Milestone D1 only:

```text
src/oma_ber/bandwidth.py
tests/test_bandwidth.py
examples/06_bandwidth_penalty.py
```

This adds the lowest-risk bridge from scalar OMA-to-BER to dynamic effects. It lets the current BER flow answer:

```text
How much does finite receiver bandwidth reduce effective OMA?
```

Do not start with full waveform simulation or S-parameter import. Those are higher-risk and need the simpler bandwidth penalty tests as reference points.
