# OMA-to-BER / Photodiode Receiver Analysis Tool

Si photonics receiver design discussion for optical NRZ/OOK links.

This repository provides a Python package for converting optical modulation
amplitude (OMA) and extinction ratio (ER) into photodiode current levels,
level-dependent noise, receiver Q, BER, and required OMA. The main calculation
is a simplified scalar Gaussian-noise model, with separate helpers for bandwidth
penalty, deterministic waveform/ISI checks, S-parameter filtering, link budget,
plotting, simplified photodiode saturation, and a physically separated
time-domain receiver-eye path with PSD-based Gaussian statistical analysis.

The physical equations and model assumptions are documented separately in
[MechanismExplanation.md](docs/MechanismExplanation.md). The numerical details
of the Gaussian BER threshold are covered in
[explanation.md](docs/explanation.md).

## Current Scope

Implemented core flow:

```text
OMA, ER
  -> optical 0/1 powers
  -> photodiode current 0/1 levels
     (linear, or tanh-compressed when saturation_power_w is provided)
  -> shot, TIA, RIN, and optional shunt thermal noise
  -> optimum threshold, Q estimate, BER
  -> required OMA for a target BER
```

Additional analysis helpers are available for:

- External-laser and optical-1-referenced OOK modulator modeling.
- Per-wavelength WDM channels with AWG, coupling, and distance-based fiber losses.
- Named passive-loss level diagrams from modulator output to PD input.
- Integrated TX/link/RX analysis and RX-only photodiode comparisons.
- First-order bandwidth penalty.
- Deterministic NRZ/OOK waveform generation.
- Simple sampled-eye/ISI metrics and deterministic 2 UI eye diagrams.
- Exact-time-grid receiver eyes with separate PD-input power [W], PD current
  [A], and TIA voltage [V].
- Standards-polynomial PRBS7/9/15/31 patterns, causal PD/TIA filters,
  automatic warm-up, fractional sampling, and clock-centered current/voltage
  eyes.
- One-sided PSD propagation for photocurrent shot noise, dark-current shot
  noise, TIA input-current noise, PD shunt thermal noise, and RIN.
- Transfer-response-derived equivalent noise bandwidth, pattern-dependent
  output variance, adjacent-sample covariance, Gaussian-mixture eye density,
  and BER-based threshold/phase optimization.
- Frequency-response CSV import and impulse-response conversion.
- Simplified static photodiode saturation models, including tanh Ge PD saturation.
- Link-budget OMA at the photodiode input.
- Matplotlib plotting helpers for examples.

The tool is intended for system-level engineering discussion. It does not claim
standards compliance such as TDECQ.

## Installation

This repository is configured for Python 3.12 or newer and `uv`.

```bash
uv sync
```

The package uses a `src/oma_ber/` layout. Tests are configured to import from
`src` when run from the repository root.

## Minimal Example

```python
from oma_ber import Photodiode, Receiver, calculate_ber_from_oma

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

result = calculate_ber_from_oma(
    oma_dbm=-10.0,
    er_db=6.0,
    pd=pd,
    rx=rx,
)

print(result["q_rx"], result["ber"])
```

The example parameters are illustrative, not specifications.

## WDM TX / Fiber Link / RX Level-Diagram Example

The WDM flow evaluates every wavelength through the requested reference
planes:

```text
External Laser
  -> laser-to-MOD coupling loss
  -> OOK MOD (optical-1-referenced IL)
  -> AWG split/insertion loss
  -> TX fiber coupling loss
  -> optical-fiber propagation loss
  -> RX fiber coupling loss
  -> optional design margin
  -> PD-input P0/P1
  -> PD/TIA Q and BER
```

Run the multi-wavelength example from the repository root:

```bash
uv sync
uv run python examples/13_tx_link_rx_level_diagram.py
```

Edit the per-wavelength inputs in each `WdmChannel`:

```python
channel = WdmChannel(
    name="lambda_1310",
    wavelength_nm=1310.0,
    laser=ExternalLaser(output_power_dbm=-4.0),
    modulator=OOKModulator(insertion_loss_db=3.0, er_db=6.0),
    laser_to_modulator_coupling_loss_db=1.0,
    awg_split_loss_db=2.0,
    tx_fiber_coupling_loss_db=1.0,
    rx_fiber_coupling_loss_db=1.0,
    fiber_attenuation_db_per_km=None,  # use the common fiber default
)

fiber_link = OpticalFiberLink(
    length_km=10.0,
    attenuation_db_per_km=0.35,
    additional_loss_db=0.2,
)
```

Each `WdmChannel` defines wavelength, External Laser power, MOD IL/ER, AWG
loss, and TX/RX coupling losses. `OpticalFiberLink` defines communication
distance and default attenuation. Set `fiber_attenuation_db_per_km` on a
channel when attenuation must differ by wavelength:

```text
fiber_loss_db
  = length_km * attenuation_db_per_km
    + additional_loss_db
```

The example calls `analyze_wdm_link()` and prints CW power before the MOD and
P0/P1/OMA after the MOD for every channel and reference plane. A typical
1310-nm result includes:

```text
External Laser:      -4.00 dBm
MOD output OMA:      -9.26 dBm
Fiber length:        10.0 km
Fiber loss:           3.70 dB
PD-input OMA:        -17.46 dBm
PD-input Q:            4.455
PD-input BER:          4.187e-06
```

`specified_oma_dbm` is an optional consistency value. Calculated levels come
from laser power, pre-MOD coupling, optical-1 MOD IL, and ER; specified OMA
does not override them.

For RX-only PD evaluation, bypass the TX and link and use fixed PD-input
levels:

```python
from oma_ber import compare_photodiodes_at_optical_levels

comparison = compare_photodiodes_at_optical_levels(
    p0_w=pd_input_p0_w,
    p1_w=pd_input_p1_w,
    photodiodes={"baseline": pd_a, "candidate": pd_b},
    rx=rx,
)
```

This keeps P0/P1 and receiver noise settings identical while changing only
the photodiode specification.

## Physically Separated RX Eye Example

The new time-domain path starts at the PD optical-input reference plane and
keeps every physical domain explicit:

```text
PD-input P0/P1 [W]
  -> rectangular NRZ/OOK optical power [W]
  -> static PD conversion and optional saturation [A]
  -> causal dimensionless PD electrical response [A]
  -> causal TIA transimpedance response [V/A]
  -> clock-phase sweep
  -> decision-centered current and voltage eyes
```

Run the Phase 1 example:

```bash
uv run python examples/14_physical_rx_eye.py
uv run python examples/14_physical_rx_eye.py --save
```

The main inputs are deliberately separate:

```python
from oma_ber import Photodiode
from oma_ber.time_domain import (
    TimeGrid,
    analyze_deterministic_eye,
    first_order_lowpass_transfer,
    prbs_bits,
    simulate_pd_tia_waveform,
)

grid = TimeGrid(symbol_rate_baud=25e9, samples_per_symbol=16)
bits = prbs_bits(order=7, num_bits=127)
pd = Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9)

pd_response = first_order_lowpass_transfer(
    sample_rate_hz=grid.sample_rate_hz,
    bandwidth_3db_hz=14e9,
    response_kind="dimensionless",
)
tia_response = first_order_lowpass_transfer(
    sample_rate_hz=grid.sample_rate_hz,
    bandwidth_3db_hz=18e9,
    dc_gain=1.5e3,
    response_kind="transimpedance_ohm",
)

waveforms = simulate_pd_tia_waveform(
    bits=bits,
    p0_w=pd_input_p0_w,
    p1_w=pd_input_p1_w,
    time_grid=grid,
    pd=pd,
    pd_current_response=pd_response,
    tia_transimpedance_response=tia_response,
)
voltage_eye = analyze_deterministic_eye(
    waveform=waveforms.tia_output_voltage_v,
    bits=bits,
    time_grid=grid,
    amplitude_unit="V",
)
```

`TimeGrid` derives `sample_rate_hz = symbol_rate_baud * samples_per_symbol`, so
the samples/UI ratio is never silently rounded. The filter warm-up is a
periodic continuation of the analysis pattern and is removed before eye
measurement. The reported optimum phase maximizes deterministic vertical
opening; it is not yet a jitter- or CDR-aware optimum.

## Phase 2 Statistical RX Eye and BER

Phase 2 keeps the deterministic signal path above and propagates each
one-sided noise PSD from its physical generation node:

```text
photocurrent shot / Idark shot / PD shunt thermal / RIN
  -> H_PD(f) * Z_TIA(f)

TIA input-referred current noise
  -> Z_TIA(f)

pattern-conditioned mean and variance
  -> Gaussian mixture at every sampling phase
  -> BER-optimum voltage threshold and sampling phase
```

Run the statistical example:

```bash
uv run python examples/15_statistical_rx_eye.py
uv run python examples/15_statistical_rx_eye.py --save
```

The main Phase 2 inputs are supplied separately from the deterministic PD/TIA
responses:

```python
from oma_ber.time_domain import (
    TimeDomainNoiseModel,
    analyze_statistical_tia_eye,
    calculate_tia_output_noise,
)

noise = calculate_tia_output_noise(
    waveforms,
    TimeDomainNoiseModel(
        tia_input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150.0,
        include_shunt_thermal_noise=True,
    ),
)
statistical_eye = analyze_statistical_tia_eye(waveforms, noise)

print(noise.pd_tia_noise_bandwidth_hz)
print(noise.tia_noise_bandwidth_hz)
print(statistical_eye.optimum_result.threshold_v)
print(statistical_eye.optimum_result.sampling_phase_ui)
print(statistical_eye.optimum_result.ber)
```

The PD object already carried by `waveforms` supplies responsivity, dark
current, dark-current Fano factor, optional shunt resistance, temperature, and
optional static saturation. The actual discrete PD/TIA responses determine
the noise bandwidth; the Phase 2 path does not reuse the scalar receiver's
manually supplied `noise_bandwidth_hz`.

The density is an analytic Gaussian mixture over the finite input pattern. It
does not synthesize a random noisy waveform and does not include jitter, CDR,
or optical-field dispersion.

## Repository Layout

```text
src/oma_ber/
  units.py          dBm, dB, W, and RIN conversions
  modulation.py     NRZ/OOK optical level calculation
  photodiode.py     photodiode parameter dataclass
  receiver.py       receiver noise parameter dataclass
  noise.py          RMS noise calculations
  ber.py            Gaussian threshold, Q, and BER calculations
  sweep.py          top-level OMA-to-BER and required-OMA helpers
  bandwidth.py      first-order bandwidth penalty helpers
  waveform.py       deterministic NRZ/OOK waveform helpers
  isi.py            LTI filtering and sampled-eye metrics
  sparameters.py    frequency-response CSV import helpers
  saturation.py     simplified static saturation models
  link_budget.py    named-loss optical level diagrams and OMA helper
  transmitter.py    external laser and optical-1-referenced OOK modulator
  system.py         integrated TX/link/RX analysis
  wdm.py            per-wavelength WDM TX/fiber/RX level diagrams
  pd_analysis.py    receiver-boundary photodiode comparisons
  plotting.py       matplotlib plotting helpers
  time_domain/      physical-unit PD/TIA waveform and statistical-eye path
```

## Units and Naming

Function and field names include units where physical units matter:

- `oma_dbm`: OMA in dBm.
- `oma_w`, `p0_w`, `p1_w`, `pavg_w`: optical powers in watts.
- `er_db`: extinction ratio in dB.
- `er_linear`: extinction ratio as a linear power ratio.
- `responsivity_a_per_w`: photodiode responsivity in A/W.
- `dark_current_a`, `i0_a`, `i1_a`, `delta_i_a`: photocurrent-related values in amperes.
- `i0_total_a`, `i1_total_a`, `threshold_total_a`: DC current coordinates in amperes including dark current.
- `dark_current_fano_factor`: dimensionless dark-current shot-noise power multiplier; 1.0 is ideal Poisson noise.
- `sigma0_a`, `sigma1_a`: RMS current noise in amperes.
- `noise_bandwidth_hz`, `bandwidth_3db_hz`: bandwidths in Hz.
- `input_current_noise_density_a_per_sqrt_hz`: TIA input-referred current noise density in A/sqrt(Hz).
- `rin_db_per_hz`: relative intensity noise in dB/Hz.
- `shunt_resistance_ohm`, `temperature_k`: optional photodiode thermal-noise parameters.
- `wavelength_nm`: WDM channel wavelength in nm.
- `length_km`, `attenuation_db_per_km`: fiber distance and attenuation.
- `awg_split_loss_db`, `tx_fiber_coupling_loss_db`,
  `rx_fiber_coupling_loss_db`: per-channel optical power losses in dB.
- `sample_interval_s`, `unit_interval_s`, `decision_delay_s`: time-domain
  timing coordinates in seconds.
- `raw_signal_current_a`, `pd_output_current_a`: pre/post-bandwidth PD current
  waveforms in amperes.
- `tia_output_voltage_v`: TIA output waveform in volts.
- `tia_transimpedance_response`: an electrical response with gain in ohms.
- `pd_tia_noise_bandwidth_hz`, `tia_noise_bandwidth_hz`: equivalent noise
  bandwidths derived from the discrete responses in Hz.
- `threshold_v`: BER-optimum TIA-output decision threshold in volts.
- `sampling_phase_ui`: decision phase within one unit interval.

## Examples

Run examples from the repository root:

```bash
uv run python examples/01_nrz_oma_to_ber.py
uv run python examples/02_sweep_oma.py
uv run python examples/03_sweep_responsivity.py
uv run python examples/04_required_oma.py
uv run python examples/05_plot_oma_sweep.py
uv run python examples/06_bandwidth_penalty.py
uv run python examples/07_nrz_waveform.py
uv run python examples/08_bandlimited_eye_metrics.py
uv run python examples/09_sparameter_filtering.py
uv run python examples/10_saturation_sweep.py
uv run python examples/11_ge_saturation_ber.py
uv run python examples/12_eye_diagram.py
uv run python examples/13_tx_link_rx_level_diagram.py
uv run python examples/14_physical_rx_eye.py
uv run python examples/15_statistical_rx_eye.py
```

Most examples print numerical results. Plotting examples may create figures or
image files.

Examples add `src/` to `sys.path`, so the documented commands work directly
from a fresh repository checkout after `uv sync`. Package users can instead
import the same public APIs from an installed environment.

## Checks

```bash
uv run pytest
uv run ruff check .
```

## Limitations

- The main BER calculation assumes scalar Gaussian noise.
- The `time_domain` Phase 2 path combines deterministic ISI with analytic
  Gaussian noise moments and mixture BER. It does not generate random noise
  samples or include jitter, CDR, bathtub curves, BER contours, or
  standards-specific masks.
- The older `waveform.py` / `isi.py` eye helpers accept arbitrary amplitude
  units and remain simplified compatibility aids. New receiver-eye work should
  use `oma_ber.time_domain` so W, A, and V are not mixed.
- PAM4 is not implemented.
- Noise bandwidth is supplied by the user in the scalar BER path. The Phase 2
  time-domain path instead derives equivalent noise bandwidth from the actual
  discrete PD/TIA responses; neither path infers it from baud rate alone.
- Photodiode capacitance, return loss, bias, and 3 dB bandwidth are stored as
  parameters but are not automatically folded into the scalar BER result.
- If `saturation_power_w` is provided, the scalar BER calculation uses a
  simplified static tanh saturation model for photodiode current. This should be
  treated as an engineering approximation.
- WDM channels are evaluated independently. AWG crosstalk, inter-channel noise,
  fiber dispersion/nonlinearity, and aggregate WDM power effects are not yet
  modeled.
