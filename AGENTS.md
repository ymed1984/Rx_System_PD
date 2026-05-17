# AGENTS.md

## Project purpose

This repository implements a Python-based OMA-to-BER and photodiode system-level analysis tool for Si photonics receiver design.

The first deliverable is an NRZ/OOK Gaussian-noise MVP that converts optical OMA and ER into PD current levels, level-dependent noise, Q factor, BER, and required OMA for a target BER.

The long-term direction is to extend the tool toward PD design discussion: responsivity, dark current, capacitance, bandwidth penalty, saturation/power handling, RIN, TIA noise, link budget, PAM4 levels, process corners, and measurement import.

## Working rules

- Prefer small, high-confidence changes.
- Do not change public APIs without updating tests and examples.
- Do not add large dependencies unless they are clearly justified.
- Keep physical units explicit in names and docstrings, for example `oma_dbm`, `oma_w`, `bandwidth_hz`, `current_a`, `capacitance_f`.
- Never silently mix dBm, dB, W, A, Hz, GHz, or dB/Hz.
- Validate physically invalid inputs and raise clear `ValueError`s.
- Keep numerical functions deterministic and side-effect free where possible.
- Put plotting in examples or `plotting.py`; core calculation modules should not show figures or write files unless explicitly requested.
- Use Python type hints for public functions.
- Use dataclasses for compact physical parameter containers unless a stronger validation layer is explicitly required.

## Preferred Python stack

- Python 3.12 or newer.
- Package layout: `src/oma_ber/`.
- Use `numpy` and `scipy` for numerical calculation.
- Use `matplotlib` for examples.
- Use `pytest` for tests.
- Use `ruff` for linting/formatting.
- Use `mypy` only if the project is already configured for it; otherwise do not block MVP delivery on mypy.

## Required checks before final response

After modifying code, run the relevant checks that are available in the repository. Prefer:

```bash
uv run pytest
uv run ruff check .
```

If `uv` is not configured yet, use the equivalent available commands, for example:

```bash
python -m pytest
python -m ruff check .
```

If a check cannot be run because a dependency or environment is missing, report that explicitly and explain the blocker.

## Numerical correctness expectations

Implement and test the following baseline equations.

### Unit conversion

```text
W = 1e-3 * 10 ** (dBm / 10)
linear = 10 ** (dB / 10)
RIN_linear = 10 ** (RIN_dB_per_Hz / 10)
```

### NRZ/OOK optical levels

For OMA = P1 - P0 and ER = P1 / P0:

```text
P0 = OMA / (ER - 1)
P1 = ER * OMA / (ER - 1)
Pavg = (P0 + P1) / 2
```

ER must be larger than 1 in linear scale.

### PD current levels

```text
I0 = Rpd * P0
I1 = Rpd * P1
Delta_I = I1 - I0 = Rpd * OMA
```

### Noise model

Level-dependent shot noise:

```text
sigma_shot_i = sqrt(2 * q * (I_i + I_dark) * Bn)
```

TIA input-referred current noise:

```text
sigma_tia = i_n * sqrt(Bn)
```

RIN noise:

```text
sigma_rin_i = Rpd * P_i * sqrt(RIN_linear * Bn)
```

Total RMS noise:

```text
sigma_i = sqrt(sigma_shot_i**2 + sigma_tia**2 + sigma_rin_i**2)
```

### BER model

For threshold `gamma`:

```text
BER = 0.5 * [ Q((gamma - mu0) / sigma0) + Q((mu1 - gamma) / sigma1) ]
Q(x) = 0.5 * erfc(x / sqrt(2))
```

Also provide a convenience Q estimate:

```text
Q_rx = (mu1 - mu0) / (sigma1 + sigma0)
BER_Q = 0.5 * erfc(Q_rx / sqrt(2))
```

## Repository structure target

Use this structure unless the existing repository already has a better compatible structure:

```text
oma-ber-tool/
  pyproject.toml
  README.md
  AGENTS.md
  docs/
    execplans/
      OMA_BER_MVP_EXECPLAN.md
  src/
    oma_ber/
      __init__.py
      units.py
      modulation.py
      photodiode.py
      receiver.py
      noise.py
      ber.py
      sweep.py
      plotting.py
      link_budget.py
  examples/
    01_nrz_oma_to_ber.py
    02_sweep_oma.py
    03_sweep_responsivity.py
    04_required_oma.py
  tests/
    test_units.py
    test_modulation.py
    test_noise.py
    test_ber.py
    test_sweep.py
```

## Documentation expectations

- `README.md` must explain the physical meaning of OMA, ER, responsivity, noise bandwidth, TIA noise, RIN, Q, and BER.
- Public functions must have docstrings with units.
- Examples must be executable with a fresh environment.
- Any simplified physical model must be labeled as simplified and not overclaimed.

## Review expectations

Before finalizing a task, review the diff for:

- Unit mistakes.
- dB/dBm/linear conversion mistakes.
- BER expression sign mistakes.
- Threshold direction mistakes.
- Monotonicity regressions, especially BER vs OMA.
- Missing validation for invalid physical inputs.
- Hidden plotting or file I/O in core numerical functions.
- Examples that cannot run from the repository root.

## What not to do

- Do not implement PAM4 before the NRZ/OOK MVP is tested.
- Do not implement time-domain waveform simulation before the scalar Gaussian model is complete.
- Do not hide assumptions inside magic constants.
- Do not make claims of standards compliance unless explicitly implemented and tested.
- Do not make git commits unless the user explicitly asks for commits.
