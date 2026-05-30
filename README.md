# OMA-to-BER Tool

## Purpose

This package implements a Python-based optical modulation amplitude (OMA) to bit error ratio (BER) calculator for Si photonics photodiode receiver discussions.

The current MVP is a simplified scalar NRZ/OOK Gaussian-noise model. It converts optical OMA and extinction ratio (ER) into optical 0/1 powers, photodiode currents, level-dependent noise, receiver Q estimate, BER, and required OMA for a target BER.

## Installation

This repository is configured for Python 3.12 or newer and `uv`.

```bash
uv sync
```

The package uses a `src/oma_ber/` layout. Tests are configured to import from `src` when run from the repository root.

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

## Physical Model

OMA is the optical modulation amplitude in watts or dBm. For NRZ/OOK it is the difference between the optical 1 level and optical 0 level:

```text
OMA = P1 - P0
```

ER is the extinction ratio:

```text
ER = P1 / P0
```

Given OMA and ER, the model calculates:

```text
P0 = OMA / (ER - 1)
P1 = ER * OMA / (ER - 1)
Pavg = (P0 + P1) / 2
```

Photodiode responsivity converts optical power to current:

```text
I0 = Rpd * P0
I1 = Rpd * P1
Delta_I = Rpd * OMA
```

The receiver noise model is level-dependent:

```text
shot noise = sqrt(2 * q * (I + I_dark) * Bn)
PD thermal noise = sqrt(4 * k_B * T * Bn / Rsh)
TIA noise = input_current_noise_density * sqrt(Bn)
RIN noise = Rpd * P * sqrt(RIN_linear * Bn)
total noise = sqrt(shot^2 + thermal^2 + TIA^2 + RIN^2)
```

The BER model assumes two Gaussian current distributions and an electrical decision threshold. For threshold `gamma`:

```text
BER = 0.5 * [Q((gamma - mu0) / sigma0) + Q((mu1 - gamma) / sigma1)]
Q(x) = 0.5 * erfc(x / sqrt(2))
```

The convenience receiver Q estimate is:

```text
Q_rx = (mu1 - mu0) / (sigma1 + sigma0)
BER_Q = 0.5 * erfc(Q_rx / sqrt(2))
```

## Units and Conventions

Function and field names include units where physical units matter:

- `oma_dbm`: OMA in dBm.
- `oma_w`, `p0_w`, `p1_w`, `pavg_w`: optical powers in watts.
- `er_db`: extinction ratio in dB.
- `er_linear`: extinction ratio as a linear power ratio.
- `responsivity_a_per_w`: photodiode responsivity in A/W.
- `dark_current_a`, `i0_a`, `i1_a`, `delta_i_a`, `sigma0_a`, `sigma1_a`: currents or RMS current noise in amperes.
- `shunt_resistance_ohm`: photodiode shunt resistance in ohms for thermal noise, when provided.
- `temperature_k`: photodiode temperature in K for thermal noise.
- `noise_bandwidth_hz`, `bandwidth_3db_hz`: bandwidths in Hz.
- `input_current_noise_density_a_per_sqrt_hz`: TIA input-referred current noise density in A/sqrt(Hz).
- `rin_db_per_hz`: relative intensity noise in dB/Hz.

Conversions use:

```text
W = 1e-3 * 10 ** (dBm / 10)
linear = 10 ** (dB / 10)
RIN_linear = 10 ** (RIN_dB_per_Hz / 10)
```

## How to Run Examples

Run examples from the repository root:

```bash
uv run python examples/01_nrz_oma_to_ber.py
uv run python examples/02_sweep_oma.py
uv run python examples/03_sweep_responsivity.py
uv run python examples/04_required_oma.py
```

The examples print numerical results and do not write files.

## How to Run Tests

```bash
uv run pytest
uv run ruff check .
```

## Limitations

- The MVP assumes scalar Gaussian noise.
- The MVP does not include time-domain intersymbol interference directly.
- The MVP does not implement PAM4.
- Noise bandwidth is supplied by the user; it is not automatically derived from baud rate.
- Photodiode 3 dB bandwidth, capacitance, saturation power, return loss, and bias are stored as parameters but are not yet used in the BER calculation.
- Bandwidth penalty and saturation behavior are planned extensions.
- This tool does not claim standards compliance such as TDECQ.

## Planned Extensions

Planned next steps after the MVP:

1. Responsivity sweep and required OMA versus responsivity.
2. Dark-current sweep and required OMA versus dark current.
3. TIA input noise sweep.
4. Simple bandwidth penalty model.
5. Photodiode capacitance and TIA input capacitance model.
6. Saturation model with responsivity compression and bandwidth degradation.
7. PAM4 outer OMA and 4-level BER.
8. Measurement import for IV, CV, S21, responsivity, and saturation CSV data.
9. Process, temperature, and bias corners.
10. Monte Carlo and yield views.
