# OMA-to-BER / Photodiode Receiver Analysis Tool

Si photonics receiver design discussion for optical NRZ/OOK links.

This repository provides a Python package for converting optical modulation
amplitude (OMA) and extinction ratio (ER) into photodiode current levels,
level-dependent noise, receiver Q, BER, and required OMA. The main calculation
is a simplified scalar Gaussian-noise model, with separate helpers for bandwidth
penalty, deterministic waveform/ISI checks, S-parameter filtering, link budget,
plotting, and simplified photodiode saturation.

The physical equations and model assumptions are documented separately in
[MechanismExplanation.md](MechanismExplanation.md).

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
  link_budget.py    OMA link-budget helper
  plotting.py       matplotlib plotting helpers
```

## Units and Naming

Function and field names include units where physical units matter:

- `oma_dbm`: OMA in dBm.
- `oma_w`, `p0_w`, `p1_w`, `pavg_w`: optical powers in watts.
- `er_db`: extinction ratio in dB.
- `er_linear`: extinction ratio as a linear power ratio.
- `responsivity_a_per_w`: photodiode responsivity in A/W.
- `dark_current_a`, `i0_a`, `i1_a`, `delta_i_a`: currents in amperes.
- `sigma0_a`, `sigma1_a`: RMS current noise in amperes.
- `noise_bandwidth_hz`, `bandwidth_3db_hz`: bandwidths in Hz.
- `input_current_noise_density_a_per_sqrt_hz`: TIA input-referred current noise density in A/sqrt(Hz).
- `rin_db_per_hz`: relative intensity noise in dB/Hz.
- `shunt_resistance_ohm`, `temperature_k`: optional photodiode thermal-noise parameters.

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
```

Most examples print numerical results. Plotting examples may create figures or
image files.

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
