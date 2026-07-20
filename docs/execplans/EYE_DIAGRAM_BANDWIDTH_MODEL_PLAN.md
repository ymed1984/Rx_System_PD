# Eye Diagram Bandwidth Model Plan

## Purpose

This document originally recorded a future implementation plan for improving
the legacy eye-diagram bandwidth model. Implementation was explicitly requested
later, and the physical direction changed from filtering an arbitrary-amplitude
waveform to separating PD-input optical power, PD current, and TIA voltage.

The legacy plan below is retained as design history. New implementation work
must use `src/oma_ber/time_domain/` and the status update in this document.
Phase 3以降の実行順序と合格条件は、
[`EYE_DIAGRAM_PHASE3_AND_LATER_EXECPLAN.md`](EYE_DIAGRAM_PHASE3_AND_LATER_EXECPLAN.md)
を正本とします。

## Physical-Model Status Update

Phase 1 is now implemented:

- Exact `TimeGrid` with integer samples/UI.
- Standards-polynomial PRBS7/9/15/31.
- Optical power [W] -> PD current [A] -> TIA voltage [V].
- Static PD saturation before electrical bandwidth.
- Causal PD and TIA transfer functions with automatic periodic warm-up.
- Fractional sampling and explicit decision delay.
- Clock-centered current and voltage eyes with polarity handling.
- `examples/14_physical_rx_eye.py`.

Phase 2 is now implemented:

- One-sided PSDs for photocurrent shot, dark-current shot, TIA input-current,
  PD shunt thermal, and RIN noise.
- Noise-source routing through `H_PD*Z_TIA` or `Z_TIA` according to the
  physical generation node.
- Equivalent noise bandwidth from the actual discrete impulse response.
- Pattern-dependent variance and adjacent-sample covariance.
- Gaussian-mixture statistical eye density and BER.
- BER-aware voltage-threshold and sampling-phase optimization.
- `examples/15_statistical_rx_eye.py`.

Phase 3 is now implemented:

- Physical CSV complex-response import for dimensionless PD current response
  and TIA transimpedance [V/A].
- Explicit DC/Nyquist coverage, known-delay removal, causality and FIR-tail
  diagnostics.
- Common measured FIR use for the deterministic signal and noise PSD paths.
- `examples/16_measured_response_statistical_eye.py`.

Phase 4 is now implemented:

- Residual Gaussian RJ and random-phase sinusoidal jitter in seconds and UI.
- Deterministic quadrature, periodic fractional sampling, and fixed bit labels.
- One fixed BER-optimum voltage threshold per nominal sampling phase.
- Jitter-averaged BER bathtub and `examples/17_residual_jitter_bathtub.py`.

Random noise waveform synthesis, CDR transfer, complex optical-field
propagation, BER contours, and measured-response passivity enforcement remain
future phases.

## Current Status

The repository already has the following pieces:

- `src/oma_ber/waveform.py`
  - Deterministic NRZ/OOK bit-to-level waveform generation.
- `src/oma_ber/isi.py`
  - First-order low-pass impulse response.
  - LTI filtering.
  - Symbol-center sampling.
  - Sampled-eye metrics.
- `src/oma_ber/eye.py`
  - Eye-trace slicing helpers.
- `src/oma_ber/plotting.py`
  - Deterministic eye-diagram plotting.
- `src/oma_ber/sparameters.py`
  - Frequency-response CSV import.
  - Frequency-response-to-impulse-response conversion.
- `examples/12_eye_diagram.py`
  - Demonstrates deterministic 2 UI eye visualization with first-order
    bandwidth-limited ISI.

The current implementation is limited to deterministic NRZ/OOK waveforms. It
does not include random noise, jitter, CDR behavior, bathtub curves, BER
contours, TDECQ, or PAM4 eye analysis.

## Why First-Order LPF Is Used Now

The first-order low-pass filter is used because it is the smallest model that
connects directly to a photodiode 3 dB bandwidth:

```text
|H(f)| = 1 / sqrt(1 + (f / f3dB)^2)
```

It needs only a small number of parameters:

```text
bandwidth_3db_hz
sample_rate_hz
```

This makes it suitable for visualizing the basic chain:

```text
limited PD bandwidth
  -> slower transitions
  -> intersymbol interference
  -> smaller sampled vertical eye opening
```

It is also consistent with the existing `bandwidth.py` and `isi.py` helpers.

## Limitations Of The First-Order Model

A first-order low-pass filter is not a full Ge photodiode or receiver model.
Real systems may include:

- Ge carrier transit-time effects.
- RC bandwidth from photodiode capacitance and load impedance.
- TIA input impedance and input capacitance.
- Electrode and interconnect loss.
- Package parasitics.
- Return loss and reflections.
- Peaking or equalization.
- Group-delay ripple.
- Bias-dependent bandwidth.
- Saturation-dependent bandwidth degradation.
- Measured S-parameter behavior.

Therefore, first-order LPF results should be treated as a design-intuition view,
not as calibrated device prediction.

## Future Implementation Direction

### Milestone 1: Keep First-Order LPF As Baseline

Keep the current first-order flow as the default example:

```text
NRZ/OOK bits
  -> optical waveform
  -> first-order low-pass impulse response
  -> filtered waveform
  -> 2 UI eye diagram
  -> sampled-eye metrics
```

This baseline is useful for quick comparison against more realistic models.

### Milestone 2: Add S21-Based Eye Example

Add a new example that uses measured or synthetic S21 data:

```text
frequency_hz,magnitude_db,phase_deg
```

Planned flow:

```text
frequency_response_from_csv()
  -> impulse_response_from_frequency_response()
  -> apply_lti_filter()
  -> plot_eye_diagram()
  -> sampled_eye_levels()
```

Suggested file:

```text
examples/13_sparameter_eye_diagram.py
```

This should reuse existing `sparameters.py`, `isi.py`, and `plotting.py`
functions rather than creating a new filtering path.

### Milestone 3: Compare First-Order And S21 Models

Add an example or helper that overlays or compares:

- First-order LPF eye opening.
- S21-derived eye opening.
- Sampled vertical eye opening.
- ISI penalty indicator.

Suggested comparison outputs:

```text
first_order_vertical_eye_opening
s21_vertical_eye_opening
first_order_isi_penalty_db
s21_isi_penalty_db
```

This should clarify when the simple 3 dB bandwidth model is optimistic or
pessimistic relative to measured response.

### Milestone 4: Optional Higher-Order Analytic Models

If a measured S-parameter file is not available, add optional analytic response
models such as:

- Second-order low-pass.
- Butterworth low-pass.
- Bessel low-pass for smoother group delay.
- Peaking response for TIA/equalized front-end exploration.

These should be explicit options, not silent replacements for the first-order
baseline.

Potential API shape:

```python
def analytic_lowpass_impulse_response(
    sample_rate_hz: float,
    bandwidth_3db_hz: float,
    num_taps: int,
    response_kind: str = "first_order",
) -> np.ndarray:
    ...
```

Do not add this abstraction until there is a clear need for multiple analytic
models.

### Milestone 5: Documentation And Guardrails

Update documentation to make model boundaries explicit:

- First-order LPF is a baseline engineering approximation.
- S21-derived filtering is closer to measured linear behavior, but still
  deterministic and linear.
- Neither path includes random noise, jitter, CDR, nonlinearity during
  filtering, or compliance-specific metrics unless separately implemented.

Documentation targets:

- `README.md`
- `docs/MechanismExplanation.md`
- Example docstrings and printed output.

## Tests To Add When Implemented

For an S21-based eye example:

- CSV response loads with required columns.
- Impulse response has expected length.
- Filtered waveform length is preserved.
- Eye diagram plotting returns a matplotlib `Axes`.
- Sampled-eye metrics are finite for an open eye.
- Invalid CSV column names raise `ValueError`.
- Invalid sample rate relative to max frequency raises `ValueError`.

For higher-order analytic models:

- Invalid `response_kind` raises `ValueError`.
- Invalid physical inputs raise `ValueError`.
- First-order mode remains numerically consistent with existing
  `first_order_lowpass_impulse_response`.

## Historical Recommended Near-Term Choice

The next practical step should be an S21-based eye diagram example, not a large
new model framework. The repository already has the necessary low-level pieces,
so the future implementation can be small and focused:

```text
examples/13_sparameter_eye_diagram.py
```

That recommendation has been superseded. PSD-based statistical receiver noise
and BER-aware phase/threshold analysis are now implemented on top of the
physical-unit Phase 1 path. A next milestone should be selected explicitly
between jitter/CDR-aware sampling and measured S-parameter integration. The
latter must first make response kind, reference impedance, DC extrapolation,
causality, passivity, and group-delay alignment explicit.
