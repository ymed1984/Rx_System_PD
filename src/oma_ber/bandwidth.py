"""Simple bandwidth penalty helpers for scalar OMA-to-BER analysis."""

import numpy as np


def first_order_lowpass_mag(
    frequency_hz: float | np.ndarray,
    bandwidth_3db_hz: float,
) -> float | np.ndarray:
    """Calculate first-order low-pass magnitude at frequency_hz.

    frequency_hz and bandwidth_3db_hz are in Hz. The returned magnitude is
    linear voltage/current amplitude gain, not dB.
    """
    if bandwidth_3db_hz <= 0:
        msg = "bandwidth_3db_hz must be positive."
        raise ValueError(msg)

    frequency = np.asarray(frequency_hz, dtype=float)
    if np.any(frequency < 0):
        msg = "frequency_hz must be non-negative."
        raise ValueError(msg)

    magnitude = 1 / np.sqrt(1 + (frequency / bandwidth_3db_hz) ** 2)
    if np.isscalar(frequency_hz):
        return float(magnitude)
    return magnitude


def bandwidth_penalty_db(
    signal_frequency_hz: float,
    bandwidth_3db_hz: float,
) -> float:
    """Calculate OMA penalty in dB from first-order bandwidth attenuation.

    signal_frequency_hz and bandwidth_3db_hz are in Hz. The penalty is
    non-negative and should be subtracted from OMA in dBm.
    """
    magnitude = first_order_lowpass_mag(signal_frequency_hz, bandwidth_3db_hz)
    return float(-20 * np.log10(magnitude))


def apply_oma_penalty_db(
    oma_dbm: float,
    penalty_db: float,
) -> float:
    """Apply a non-negative OMA penalty in dB to OMA in dBm."""
    if penalty_db < 0:
        msg = "penalty_db must be non-negative."
        raise ValueError(msg)
    return oma_dbm - penalty_db
