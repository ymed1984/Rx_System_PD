# OMA-to-BER / Photodiode Receiver Analysis Tool

Si photonics receiver design discussion for optical NRZ/OOK links.

This repository provides a Python package for converting optical modulation
amplitude (OMA) and extinction ratio (ER) into photodiode current levels,
level-dependent noise, receiver Q, BER, and required OMA. The main calculation
is a simplified scalar Gaussian-noise model, with separate helpers for bandwidth
penalty, deterministic waveform/ISI checks, S-parameter filtering, link budget,
plotting, and simplified photodiode saturation.

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
- Time-domain waveform and ISI helpers are deterministic analysis aids, not a
  full statistical waveform simulator.
- PAM4 is not implemented.
- Noise bandwidth is supplied by the user; it is not automatically derived from
  baud rate.
- Photodiode capacitance, return loss, bias, and 3 dB bandwidth are stored as
  parameters but are not automatically folded into the scalar BER result.
- If `saturation_power_w` is provided, the scalar BER calculation uses a
  simplified static tanh saturation model for photodiode current. This should be
  treated as an engineering approximation.
- WDM channels are evaluated independently. AWG crosstalk, inter-channel noise,
  fiber dispersion/nonlinearity, and aggregate WDM power effects are not yet
  modeled.
