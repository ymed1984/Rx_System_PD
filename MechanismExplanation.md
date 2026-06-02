# Mechanism Explanation

This document explains the physical model used by the OMA-to-BER tool. The
implementation is intentionally simple and deterministic: it is a scalar
NRZ/OOK Gaussian-noise model for receiver-level discussion, not a full
time-domain compliance model.

## 1. Unit Conversions

Optical powers may be specified in dBm, while ratios such as ER and RIN use dB
forms. The code keeps these conversions explicit:

```text
W = 1e-3 * 10 ** (dBm / 10)
dBm = 10 * log10(W / 1e-3)
linear = 10 ** (dB / 10)
dB = 10 * log10(linear)
RIN_linear_per_Hz = 10 ** (RIN_dB_per_Hz / 10)
```

Invalid logarithmic conversions are rejected when the linear value must be
positive, for example `watt_to_dbm(0)` and `linear_to_db(0)`.

## 2. NRZ/OOK Optical Levels

For NRZ/OOK, OMA is the difference between the optical 1 level and the optical
0 level:

```text
OMA = P1 - P0
```

The extinction ratio is:

```text
ER = P1 / P0
```

Given OMA and ER in linear scale:

```text
P0 = OMA / (ER - 1)
P1 = ER * OMA / (ER - 1)
Pavg = (P0 + P1) / 2
```

`ER` must be larger than 1. `OMA` must be positive. These checks prevent
unphysical cases such as a zero optical modulation amplitude or a 1 level that
is not larger than the 0 level.

## 3. Photodiode Current Levels

The scalar model converts optical power to photocurrent with responsivity:

```text
I0 = Rpd * P0
I1 = Rpd * P1
Delta_I = I1 - I0 = Rpd * OMA
```

Where:

- `Rpd` is photodiode responsivity in A/W.
- `P0` and `P1` are optical powers in W.
- `I0`, `I1`, and `Delta_I` are currents in A.

The photodiode parameter container also stores optional quantities such as
`bandwidth_3db_hz`, `capacitance_f`, `saturation_power_w`, `return_loss_db`,
`bias_v`, `shunt_resistance_ohm`, and `temperature_k`. In the main scalar BER
calculation, responsivity and dark current directly affect the result.

If `saturation_power_w` is not provided, the current conversion remains linear.
If `saturation_power_w` is provided, the main BER calculation uses the
simplified static tanh Ge PD saturation model:

```text
I = Rpd * Psat * tanh(P / Psat)
```

This means `I0`, `I1`, and `Delta_I` become compressed values. The result
dictionary also reports the corresponding linear currents and a current OMA
compression indicator:

```text
compression_dB = 20 * log10(Delta_I_linear / Delta_I_compressed)
```

Because both optical levels approach the same saturation current at very high
input power, increasing OMA does not necessarily keep improving the electrical
eye when this static model is enabled. In strong saturation, `Delta_I` can
shrink and BER can degrade.

Optional shunt resistance and temperature affect the result only when thermal
noise is enabled by providing `shunt_resistance_ohm`.

## 4. Noise Model

The receiver noise is level dependent because shot noise and RIN noise depend
on optical/current level. For each level `i`, the total RMS current noise is:

```text
sigma_i = sqrt(
    sigma_shot_i**2
  + sigma_tia**2
  + sigma_rin_i**2
  + sigma_thermal**2
)
```

All terms are RMS current noise in amperes.

### Shot Noise

Shot noise is calculated from signal photocurrent plus dark current:

```text
sigma_shot_i = sqrt(2 * q * (I_i + I_dark) * Bn)
```

Where:

- `q = 1.602176634e-19 C`.
- `I_i` is the level current in A.
- `I_dark` is dark current in A.
- `Bn` is receiver noise bandwidth in Hz.

Dark current contributes to shot noise even when the optical level is low.

### TIA Input-Referred Current Noise

TIA input-referred noise is modeled as a white current noise density integrated
over the supplied noise bandwidth:

```text
sigma_tia = i_n * sqrt(Bn)
```

Where `i_n` is in A/sqrt(Hz). The model assumes the user has already chosen an
appropriate equivalent noise bandwidth.

### RIN Noise

Relative intensity noise is applied as optical-power-proportional current noise:

```text
sigma_rin_i = Rpd * P_i * sqrt(RIN_linear_per_Hz * Bn)
```

If `rin_db_per_hz` is `None`, the RIN term is zero.

When tanh saturation is enabled, the shot-noise term uses the compressed
photocurrent. The current RIN term remains based on the small-signal
responsivity expression above. A differential saturated RIN model using
`dI/dP` would be a more detailed extension.

### Photodiode Shunt Thermal Noise

If `shunt_resistance_ohm` is provided, the model includes the simplified
Johnson-Nyquist current noise of the photodiode shunt resistance:

```text
sigma_thermal = sqrt(4 * k_B * T * Bn / Rsh)
```

Where:

- `k_B = 1.380649e-23 J/K`.
- `T` is temperature in K.
- `Rsh` is shunt resistance in ohms.

If `shunt_resistance_ohm` is not provided, this term is zero.

## 5. Gaussian BER Model

The model treats the sampled receiver current as two Gaussian distributions:

```text
level 0: N(mu0, sigma0)
level 1: N(mu1, sigma1)
```

For NRZ/OOK:

```text
mu0 = I0
mu1 = I1
```

Given an electrical decision threshold `gamma`, the two conditional error
probabilities are:

```text
P(error | 0) = Q((gamma - mu0) / sigma0)
P(error | 1) = Q((mu1 - gamma) / sigma1)
```

With equal 0/1 priors:

```text
BER = 0.5 * [
    Q((gamma - mu0) / sigma0)
  + Q((mu1 - gamma) / sigma1)
]
```

The Gaussian Q function is:

```text
Q(x) = 0.5 * erfc(x / sqrt(2))
```

By default, the code numerically finds the BER-minimizing threshold between
`mu0` and `mu1`. If threshold optimization is disabled, the midpoint threshold
is used:

```text
gamma = (mu0 + mu1) / 2
```

## 6. Convenience Q Estimate

The code also reports the common receiver Q estimate:

```text
Q_rx = (mu1 - mu0) / (sigma1 + sigma0)
BER_Q = 0.5 * erfc(Q_rx / sqrt(2))
```

This `BER_Q` is a convenience approximation. The reported `ber` from
`calculate_ber_from_oma` is calculated from the Gaussian threshold expression
above.

## 7. Required OMA Search

`required_oma_dbm` searches for the OMA that reaches a target BER:

```text
BER(OMA_required) = target_BER
```

The implementation evaluates BER at the user-supplied lower and upper OMA
bounds, then uses a scalar root finder when the target is bracketed. If the
target BER is already met at the lower bound, the lower bound is returned. If
the target cannot be reached at the upper bound, a `ValueError` is raised.

The search assumes BER generally improves as OMA increases under the scalar
model. Tests cover this monotonic behavior for representative conditions.

## 8. Bandwidth Penalty Helper

The simple bandwidth helper uses the magnitude of a first-order low-pass:

```text
|H(f)| = 1 / sqrt(1 + (f / f3dB)**2)
```

The corresponding OMA penalty is:

```text
penalty_dB = -20 * log10(|H(f)|)
effective_OMA_dBm = OMA_dBm - penalty_dB
```

This is a scalar engineering helper. It does not replace waveform-level ISI or
standards-specific transmitter/receiver penalty metrics.

## 9. Waveform, ISI, and Eye-Diagram Helpers

The waveform helpers generate deterministic NRZ/OOK sample sequences from
0/1 bits. The ISI helpers can convolve a waveform with a supplied impulse
response and sample one point per symbol.

ISI means inter-symbol interference. In this receiver context, it means part of
one symbol response remains into later symbol periods, so the sampled value for
the current bit depends on the previous bit pattern.

Ideally, every symbol-center sample would map cleanly to one of two levels:

```text
0 bit -> P0 or I0
1 bit -> P1 or I1
```

With limited photodiode or receiver bandwidth, transitions are slower. After a
long run of zeros, a following one may not have reached the full high level at
the sampling instant:

```text
...0001
```

After a long run of ones, a following zero may not have settled to the full low
level:

```text
...1110
```

Therefore, samples for the same logical bit can differ depending on the
preceding symbols:

```text
...0001 sample for 1 != ...1111 sample for 1
...1110 sample for 0 != ...0000 sample for 0
```

This pattern dependence closes the eye even without random noise.

The eye-diagram helpers slice the deterministic waveform into overlapping
2 UI traces and plot those traces with matplotlib. This visualizes vertical eye
closure caused by bandwidth-limited ISI, but it does not add random noise,
jitter, CDR behavior, or BER contours.

The sampled-eye metrics include:

- Mean zero level.
- Mean one level.
- Minimum sampled one level.
- Maximum sampled zero level.
- Vertical eye opening.
- A simplified ISI penalty indicator.

The vertical eye opening is calculated from the worst sampled levels:

```text
vertical_eye_opening = min_one_level - max_zero_level
```

If:

```text
min_one_level > max_zero_level
```

the sampled eye remains open. If:

```text
min_one_level <= max_zero_level
```

the eye is closed in this deterministic sampled-eye sense. In that case, a
single threshold cannot separate all sampled zeros and ones even before adding
random noise.

The ISI penalty indicator is:

```text
isi_penalty_db = 20 * log10(mean_level_separation / vertical_eye_opening)
```

It is an engineering indicator of how much the usable worst-case eye opening
has shrunk relative to the average zero/one level separation. For example, if:

```text
mean_level_separation = 1.0
vertical_eye_opening = 0.8
```

then:

```text
isi_penalty_db = 20 * log10(1.0 / 0.8) ~= 1.94 dB
```

If the vertical eye opening is zero or negative, the indicator is infinite.
This is an engineering diagnostic, not a compliance metric.

ISI affects OMA-to-BER interpretation because optical OMA is the transmitted
level difference:

```text
OMA = P1 - P0
```

but the receiver decision uses sampled levels after the bandwidth-limited
response. With strong ISI, the usable sampled opening can be much smaller than
the nominal OMA-derived level separation:

```text
min_one_level - max_zero_level < P1 - P0
```

In current units, this reduces the effective decision margin, lowers the
effective Q, and can worsen BER. The deterministic ISI helpers do not directly
replace the scalar Gaussian BER calculation; they show how bandwidth-limited
memory can reduce the eye opening that a later noise/BER model would use.

## 10. S-Parameter / Frequency Response Helper

`frequency_response_from_csv` reads a CSV file with:

```text
frequency_hz,magnitude_db,phase_deg
```

The magnitude is interpreted as voltage/current amplitude gain in dB, not
optical power dB. The helper converts magnitude and phase into a complex linear
response:

```text
H(f) = 10 ** (magnitude_db / 20) * exp(j * phase_rad)
```

`impulse_response_from_frequency_response` interpolates this response onto a
real-FFT frequency grid and converts it to a real impulse response. This is a
deterministic filtering helper for exploration and examples.

## 11. Saturation Helpers and Ge PD Integration

The saturation module contains simplified static models that convert optical
power into compressed photocurrent or effective responsivity. These models are
not device-physics simulations.

The main scalar OMA-to-BER flow uses the tanh model when
`Photodiode.saturation_power_w` is provided:

```text
I = R * Psat * tanh(P / Psat)
R_eff = I / P
```

Soft-clipping model:

```text
I = R * P / (1 + (P / Psat)**order)**(1 / order)
```

Exponential model:

```text
I = R * Psat * (1 - exp(-P / Psat))
```

Rational model:

```text
R_eff = R / (1 + (P / Psat)**order)
I = P * R_eff
```

These helpers are useful for comparing qualitative compression behavior, but
model parameters should not be interpreted as measured device limits unless
they have been calibrated against data.

## 12. Link Budget Helper

The link-budget helper calculates OMA at the photodiode input:

```text
PD_input_OMA_dBm = TX_OMA_dBm - sum(losses_dB) - margin_dB
```

Losses and margin must be non-negative dB quantities.

## 13. Main Assumptions and Boundaries

- The main BER model is scalar and memoryless.
- Noise distributions are Gaussian.
- 0 and 1 symbols are assumed equally likely in the BER expression.
- Noise bandwidth is a user input.
- No PAM4 BER model is implemented.
- No standards-specific compliance calculation is implemented.
- Waveform, S-parameter, bandwidth, and saturation helpers are separate
  engineering aids unless explicitly combined by user code.
