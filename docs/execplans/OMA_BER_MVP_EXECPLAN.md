# OMA-to-BER / PD Analysis Tool MVP ExecPlan

## Purpose

Build a Python package that calculates BER from optical OMA for a Si photonics photodiode receiver model.

The MVP target is NRZ/OOK with Gaussian level-dependent noise. The tool should be useful for system-level PD design discussions by making the following chain explicit:

```text
OMA, ER
  -> optical levels P0, P1, Pavg
  -> PD current levels I0, I1
  -> shot noise, dark-current shot noise, TIA noise, RIN noise
  -> threshold, Q factor, BER
  -> required OMA for target BER
```

This MVP must be numerically reliable, unit-explicit, tested, and easy to extend toward PD bandwidth, capacitance, saturation, PAM4, link budget, and measurement import.

## Non-goals for MVP

The MVP must not attempt to solve all receiver design problems at once.

Do not implement the following before the NRZ/OOK scalar model and tests are complete:

- PAM4 full BER model.
- Time-domain eye simulation.
- S-parameter convolution.
- Real TIA circuit macromodels.
- Full PD process/device TCAD import.
- Standards-specific compliance such as TDECQ.

Small stubs or TODO comments are acceptable only if they do not confuse users or tests.

## Milestone 0: Inspect repository and set up minimal project

### Tasks

1. Inspect current repository state.
2. If no Python package exists, create a `src/oma_ber/` package layout.
3. Create or update `pyproject.toml` with minimal dependencies.
4. Prefer dependencies:
   - `numpy`
   - `scipy`
   - `matplotlib`
   - `pytest`
   - `ruff`
5. Add `README.md` if missing.
6. Add this ExecPlan under `docs/execplans/OMA_BER_MVP_EXECPLAN.md` if not already present.

### Acceptance criteria

- `python -c "import oma_ber"` works when package is installed or when run in the configured environment.
- Test framework can be invoked.
- No unnecessary dependencies are added.

## Milestone 1: Unit conversion module

### Files

- `src/oma_ber/units.py`
- `tests/test_units.py`

### Required functions

```python
def dbm_to_watt(dbm: float) -> float: ...
def watt_to_dbm(watt: float) -> float: ...
def db_to_linear(db: float) -> float: ...
def linear_to_db(value: float) -> float: ...
def rin_db_per_hz_to_linear(rin_db_per_hz: float) -> float: ...
```

### Validation

- `watt_to_dbm` must reject non-positive W.
- `linear_to_db` must reject non-positive values.

### Tests

- `0 dBm -> 1e-3 W`.
- `-10 dBm -> 1e-4 W`.
- `3 dB` approximately maps to factor `~2`.
- round-trip tests for dBm/W and dB/linear.

## Milestone 2: NRZ/OOK optical levels

### Files

- `src/oma_ber/modulation.py`
- `tests/test_modulation.py`

### Required API

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class NRZLevels:
    p0_w: float
    p1_w: float
    pavg_w: float
    oma_w: float
    er_linear: float
    er_db: float


def nrz_levels_from_oma_er(oma_w: float, er_db: float) -> NRZLevels: ...
```

### Equations

```text
ER_linear = 10 ** (ER_dB / 10)
P0 = OMA / (ER_linear - 1)
P1 = ER_linear * P0
Pavg = (P0 + P1) / 2
```

### Validation

- `oma_w` must be positive.
- `ER_linear` must be larger than 1.

### Tests

- Check that `P1 - P0 == OMA` within tolerance.
- Check that `P1 / P0 == ER_linear` within tolerance.
- Invalid ER raises `ValueError`.

## Milestone 3: Photodiode and receiver models

### Files

- `src/oma_ber/photodiode.py`
- `src/oma_ber/receiver.py`
- `tests/test_noise.py`

### Required classes

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class Photodiode:
    responsivity_a_per_w: float
    dark_current_a: float
    bandwidth_3db_hz: float | None = None
    capacitance_f: float | None = None
    saturation_power_w: float | None = None
    return_loss_db: float | None = None
    bias_v: float | None = None

@dataclass(frozen=True)
class Receiver:
    noise_bandwidth_hz: float
    input_current_noise_density_a_per_sqrt_hz: float
    rin_db_per_hz: float | None = None
```

### Validation

- Responsivity must be positive.
- Dark current must be non-negative.
- Noise bandwidth must be positive.
- Input current noise density must be non-negative.
- Optional positive quantities must be positive if provided.

## Milestone 4: Noise model

### Files

- `src/oma_ber/noise.py`
- `tests/test_noise.py`

### Required functions

```python
def shot_noise_rms_a(current_a: float, dark_current_a: float, bandwidth_hz: float) -> float: ...
def tia_noise_rms_a(input_noise_density_a_per_sqrt_hz: float, bandwidth_hz: float) -> float: ...
def rin_noise_rms_a(responsivity_a_per_w: float, optical_power_w: float, rin_db_per_hz: float | None, bandwidth_hz: float) -> float: ...
def total_noise_rms_a(optical_power_w: float, photocurrent_a: float, pd: Photodiode, rx: Receiver) -> float: ...
```

### Constants

Use the exact elementary charge:

```python
Q_E = 1.602176634e-19
```

### Tests

- Shot noise increases with current.
- Shot noise increases with bandwidth.
- TIA noise is zero if input noise density is zero.
- RIN noise is zero if `rin_db_per_hz is None`.
- Total noise is at least as large as each individual contribution.

## Milestone 5: BER model

### Files

- `src/oma_ber/ber.py`
- `tests/test_ber.py`

### Required functions

```python
def qfunc(x: float | np.ndarray) -> float | np.ndarray: ...
def ber_for_threshold(mu0_a: float, mu1_a: float, sigma0_a: float, sigma1_a: float, threshold_a: float) -> float: ...
def optimum_threshold(mu0_a: float, mu1_a: float, sigma0_a: float, sigma1_a: float) -> float: ...
def q_from_levels(mu0_a: float, mu1_a: float, sigma0_a: float, sigma1_a: float) -> float: ...
def ber_from_q(q_rx: float) -> float: ...
def ber_from_gaussian_levels(mu0_a: float, mu1_a: float, sigma0_a: float, sigma1_a: float, optimize_threshold: bool = True) -> dict: ...
```

### Formula

```text
Q(x) = 0.5 * erfc(x / sqrt(2))
BER = 0.5 * [Q((threshold - mu0) / sigma0) + Q((mu1 - threshold) / sigma1)]
Q_rx = (mu1 - mu0) / (sigma1 + sigma0)
BER_Q = 0.5 * erfc(Q_rx / sqrt(2))
```

### Important sign convention

For transmitted 0, an error occurs when the observed current is above the threshold:

```text
P(error | 0) = Q((threshold - mu0) / sigma0)
```

For transmitted 1, an error occurs when the observed current is below the threshold:

```text
P(error | 1) = Q((mu1 - threshold) / sigma1)
```

### Tests

- `qfunc(0) == 0.5`.
- BER decreases when current separation increases.
- BER increases when noise increases.
- Optimum threshold lies between `mu0` and `mu1`.
- Equal-noise optimum threshold is approximately midpoint.

## Milestone 6: Top-level OMA-to-BER calculation

### Files

- `src/oma_ber/__init__.py`
- `src/oma_ber/sweep.py`
- `tests/test_sweep.py`

### Required top-level API

```python
def calculate_ber_from_oma(
    oma_dbm: float,
    er_db: float,
    pd: Photodiode,
    rx: Receiver,
    optimize_threshold: bool = True,
) -> dict: ...
```

The returned dictionary must include at least:

```text
oma_dbm
oma_w
er_db
p0_w
p1_w
pavg_w
i0_a
i1_a
delta_i_a
sigma0_a
sigma1_a
threshold_a
q_rx
ber
ber_from_q
```

### Required sweep API

```python
def sweep_oma(
    oma_dbm_values: np.ndarray,
    er_db: float,
    pd: Photodiode,
    rx: Receiver,
) -> list[dict]: ...


def required_oma_dbm(
    target_ber: float,
    er_db: float,
    pd: Photodiode,
    rx: Receiver,
    oma_dbm_min: float = -40.0,
    oma_dbm_max: float = 10.0,
) -> float: ...
```

### Tests

- BER decreases as OMA increases for a fixed model.
- Required OMA returns a value whose computed BER is close to target.
- If target BER cannot be reached in the specified search range, raise `ValueError` with a clear message.

## Milestone 7: Link-budget helper

### Files

- `src/oma_ber/link_budget.py`
- `tests/test_link_budget.py`

### Required functions

```python
def pd_input_oma_dbm(
    tx_oma_dbm: float,
    losses_db: list[float] | tuple[float, ...],
    margin_db: float = 0.0,
) -> float: ...
```

### Formula

```text
OMA_PD_dBm = OMA_TX_dBm - sum(losses_db) - margin_db
```

### Why this matters

PD design should be evaluated using OMA at the PD input, not only transmitter-side OMA.

## Milestone 8: Examples

### Files

- `examples/01_nrz_oma_to_ber.py`
- `examples/02_sweep_oma.py`
- `examples/03_sweep_responsivity.py`
- `examples/04_required_oma.py`

### Example default parameters

Use physically plausible defaults, but label them as examples, not specifications:

```python
pd = Photodiode(
    responsivity_a_per_w=0.8,
    dark_current_a=1e-9,
    bandwidth_3db_hz=40e9,
)

rx = Receiver(
    noise_bandwidth_hz=25e9,
    input_current_noise_density_a_per_sqrt_hz=10e-12,
    rin_db_per_hz=-150,
)
```

### Plotting rule

Use `matplotlib` only. Do not require seaborn.

## Milestone 9: README

### Required README sections

1. Purpose.
2. Installation.
3. Minimal example.
4. Physical model.
5. Units and conventions.
6. How to run examples.
7. How to run tests.
8. Limitations.
9. Planned extensions.

### Limitations to state clearly

- MVP assumes Gaussian noise.
- MVP does not include time-domain ISI directly.
- MVP does not implement PAM4 yet.
- Noise bandwidth is supplied by the user; it is not automatically derived from baud rate.
- Bandwidth penalty and saturation are planned extensions.

## Milestone 10: Review and finish

### Required final checks

Run:

```bash
uv run pytest
uv run ruff check .
```

If `uv` is unavailable, run the closest equivalents.

### Final response from Codex must include

- Summary of implemented modules.
- Commands run and pass/fail status.
- Any unresolved limitations.
- Any assumptions made.
- Suggested next milestone.

## Future extension plan

After the MVP is correct and tested, add features in this order:

1. `responsivity` sweep and required OMA vs responsivity.
2. `dark_current` sweep and required OMA vs dark current.
3. TIA input noise sweep.
4. Simple bandwidth penalty model.
5. PD capacitance and TIA input capacitance model.
6. Saturation model: responsivity compression and bandwidth degradation.
7. PAM4 outer OMA and 4-level BER.
8. Measurement import: CSV for IV, CV, S21, responsivity, saturation.
9. Process/temperature/bias corners.
10. Monte Carlo/yield view.
