"""Linear filtering and deterministic ISI metric helpers."""

import numpy as np


def first_order_lowpass_impulse_response(
    sample_rate_hz: float,
    bandwidth_3db_hz: float,
    num_taps: int,
) -> np.ndarray:
    """Create a causal first-order low-pass impulse response.

    sample_rate_hz and bandwidth_3db_hz are in Hz. num_taps is the number of
    discrete-time impulse-response samples. The returned impulse response is
    normalized to unity DC gain, so a constant waveform keeps the same level.
    """
    if sample_rate_hz <= 0:
        msg = "sample_rate_hz must be positive."
        raise ValueError(msg)
    if bandwidth_3db_hz <= 0:
        msg = "bandwidth_3db_hz must be positive."
        raise ValueError(msg)
    if num_taps <= 0:
        msg = "num_taps must be positive."
        raise ValueError(msg)

    time_s = np.arange(num_taps, dtype=float) / sample_rate_hz
    impulse = np.exp(-2 * np.pi * bandwidth_3db_hz * time_s)
    return impulse / np.sum(impulse)


def apply_lti_filter(waveform: np.ndarray, impulse_response: np.ndarray) -> np.ndarray:
    """Apply a causal LTI impulse response to a one-dimensional waveform.

    waveform and impulse_response use arbitrary but consistent amplitude units,
    such as optical power in W or current in A. The returned waveform has the
    same number of samples as waveform.
    """
    waveform_values = np.asarray(waveform, dtype=float)
    impulse_values = np.asarray(impulse_response, dtype=float)

    if waveform_values.ndim != 1:
        msg = "waveform must be a one-dimensional array."
        raise ValueError(msg)
    if impulse_values.ndim != 1:
        msg = "impulse_response must be a one-dimensional array."
        raise ValueError(msg)
    if waveform_values.size == 0:
        msg = "waveform must contain at least one sample."
        raise ValueError(msg)
    if impulse_values.size == 0:
        msg = "impulse_response must contain at least one sample."
        raise ValueError(msg)

    return np.convolve(waveform_values, impulse_values, mode="full")[: waveform_values.size]


def sample_at_symbol_centers(
    waveform: np.ndarray,
    samples_per_symbol: int,
    timing_offset_samples: int | None = None,
) -> np.ndarray:
    """Sample a waveform once per symbol at symbol centers.

    samples_per_symbol is in samples/symbol. timing_offset_samples is the
    sample index within each symbol; if omitted, samples_per_symbol // 2 is
    used.
    """
    waveform_values = np.asarray(waveform, dtype=float)
    if waveform_values.ndim != 1:
        msg = "waveform must be a one-dimensional array."
        raise ValueError(msg)
    if samples_per_symbol < 2:
        msg = "samples_per_symbol must be at least 2."
        raise ValueError(msg)
    if waveform_values.size < samples_per_symbol:
        msg = "waveform must contain at least one complete symbol."
        raise ValueError(msg)

    offset = samples_per_symbol // 2 if timing_offset_samples is None else timing_offset_samples
    if offset < 0 or offset >= samples_per_symbol:
        msg = "timing_offset_samples must be in [0, samples_per_symbol)."
        raise ValueError(msg)

    num_symbols = waveform_values.size // samples_per_symbol
    indices = offset + np.arange(num_symbols) * samples_per_symbol
    return waveform_values[indices]


def sampled_eye_levels(samples: np.ndarray, bits: np.ndarray) -> dict[str, float]:
    """Calculate simple sampled NRZ/OOK eye-level metrics.

    samples and bits must have the same length. The sample amplitudes may be
    optical power in W or current in A. isi_penalty_db is an engineering
    indicator calculated as 20*log10(mean_level_separation / vertical_eye_opening);
    it is not a standards-compliance metric.
    """
    sample_values = np.asarray(samples, dtype=float)
    bit_values = np.asarray(bits)

    if sample_values.ndim != 1:
        msg = "samples must be a one-dimensional array."
        raise ValueError(msg)
    if bit_values.ndim != 1:
        msg = "bits must be a one-dimensional array."
        raise ValueError(msg)
    if sample_values.size == 0:
        msg = "samples must contain at least one sample."
        raise ValueError(msg)
    if sample_values.size != bit_values.size:
        msg = "samples length must match bits length."
        raise ValueError(msg)
    if not np.all((bit_values == 0) | (bit_values == 1)):
        msg = "bits must contain only 0 and 1."
        raise ValueError(msg)
    if not np.any(bit_values == 0) or not np.any(bit_values == 1):
        msg = "bits must contain at least one 0 and one 1."
        raise ValueError(msg)

    zero_samples = sample_values[bit_values == 0]
    one_samples = sample_values[bit_values == 1]
    mean_zero_level = float(np.mean(zero_samples))
    mean_one_level = float(np.mean(one_samples))
    min_one_level = float(np.min(one_samples))
    max_zero_level = float(np.max(zero_samples))
    vertical_eye_opening = min_one_level - max_zero_level
    mean_level_separation = mean_one_level - mean_zero_level

    if vertical_eye_opening <= 0:
        isi_penalty_db = float("inf")
    else:
        isi_penalty_db = float(20 * np.log10(mean_level_separation / vertical_eye_opening))

    return {
        "mean_zero_level": mean_zero_level,
        "mean_one_level": mean_one_level,
        "min_one_level": min_one_level,
        "max_zero_level": max_zero_level,
        "vertical_eye_opening": float(vertical_eye_opening),
        "mean_level_separation": float(mean_level_separation),
        "isi_penalty_db": isi_penalty_db,
    }
