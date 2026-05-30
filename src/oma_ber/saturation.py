"""Simplified static photodiode saturation helpers."""

import numpy as np


def compressed_responsivity_a_per_w(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
    compression_order: float = 1.0,
) -> float | np.ndarray:
    """Calculate compressed photodiode responsivity in A/W.

    optical_power_w is optical input power in W, responsivity_a_per_w is the
    small-signal responsivity in A/W, and saturation_power_w is the optical
    power scale in W. This default helper uses the simplified tanh model
    I = R * Psat * tanh(P / Psat), and returns R_eff = I / P. At P = 0 W,
    R_eff is defined as the small-signal responsivity.
    """
    _validate_saturation_inputs(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
        compression_order=compression_order,
    )

    optical_power = np.asarray(optical_power_w, dtype=float)
    current = responsivity_a_per_w * saturation_power_w * np.tanh(optical_power / saturation_power_w)
    responsivity = np.full_like(optical_power, responsivity_a_per_w, dtype=float)
    np.divide(current, optical_power, out=responsivity, where=optical_power > 0)
    return _return_scalar_if_scalar_input(responsivity, optical_power_w)


def saturated_photocurrent_a(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
    compression_order: float = 1.0,
) -> float | np.ndarray:
    """Calculate saturated photodiode photocurrent in A.

    optical_power_w is optical input power in W, responsivity_a_per_w is in A/W,
    and saturation_power_w is in W. The default simplified static model is
    I = R * Psat * tanh(P / Psat), which approaches R * Psat at high optical
    power. compression_order is kept for API compatibility but is not used by
    the tanh model.
    """
    _validate_saturation_inputs(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
        compression_order=compression_order,
    )

    optical_power = np.asarray(optical_power_w, dtype=float)
    current = responsivity_a_per_w * saturation_power_w * np.tanh(optical_power / saturation_power_w)
    return _return_scalar_if_scalar_input(current, optical_power_w)


def compressed_responsivity_rational_a_per_w(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
    compression_order: float = 1.0,
) -> float | np.ndarray:
    """Calculate compressed responsivity in A/W with the rational model.

    optical_power_w is optical input power in W, responsivity_a_per_w is the
    small-signal responsivity in A/W, and saturation_power_w is in W. The model
    is R_eff = R / (1 + (P / Psat)**compression_order).
    """
    _validate_saturation_inputs(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
        compression_order=compression_order,
    )

    optical_power = np.asarray(optical_power_w, dtype=float)
    responsivity = responsivity_a_per_w / (1 + (optical_power / saturation_power_w) ** compression_order)
    return _return_scalar_if_scalar_input(responsivity, optical_power_w)


def saturated_photocurrent_rational_a(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
    compression_order: float = 1.0,
) -> float | np.ndarray:
    """Calculate photocurrent in A with the original rational model.

    optical_power_w is optical input power in W, responsivity_a_per_w is in A/W,
    and saturation_power_w is in W. The model is
    I = R * P / (1 + (P / Psat)**compression_order).
    """
    effective_responsivity = compressed_responsivity_rational_a_per_w(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
        compression_order=compression_order,
    )
    current = np.asarray(optical_power_w, dtype=float) * effective_responsivity
    return _return_scalar_if_scalar_input(current, optical_power_w)


def saturated_photocurrent_soft_clip_a(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
    compression_order: float = 1.0,
) -> float | np.ndarray:
    """Calculate photocurrent in A with a soft-clipping saturation model.

    optical_power_w is optical input power in W, responsivity_a_per_w is in A/W,
    and saturation_power_w is in W. The simplified static model is
    I = R * P / (1 + (P / Psat)**order)**(1 / order), which approaches
    R * Psat at high optical power.
    """
    _validate_saturation_inputs(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
        compression_order=compression_order,
    )

    optical_power = np.asarray(optical_power_w, dtype=float)
    current = responsivity_a_per_w * optical_power / (
        1 + (optical_power / saturation_power_w) ** compression_order
    ) ** (1 / compression_order)
    return _return_scalar_if_scalar_input(current, optical_power_w)


def saturated_photocurrent_tanh_a(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
) -> float | np.ndarray:
    """Calculate photocurrent in A with a tanh saturation model.

    optical_power_w is optical input power in W, responsivity_a_per_w is in A/W,
    and saturation_power_w is in W. The simplified static model is
    I = R * Psat * tanh(P / Psat), which approaches R * Psat at high optical
    power.
    """
    _validate_saturation_inputs(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
        compression_order=1.0,
    )

    optical_power = np.asarray(optical_power_w, dtype=float)
    current = responsivity_a_per_w * saturation_power_w * np.tanh(optical_power / saturation_power_w)
    return _return_scalar_if_scalar_input(current, optical_power_w)


def saturated_photocurrent_exponential_a(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
) -> float | np.ndarray:
    """Calculate photocurrent in A with an exponential saturation model.

    optical_power_w is optical input power in W, responsivity_a_per_w is in A/W,
    and saturation_power_w is in W. The simplified static model is
    I = R * Psat * (1 - exp(-P / Psat)), which approaches R * Psat at high
    optical power.
    """
    _validate_saturation_inputs(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
        compression_order=1.0,
    )

    optical_power = np.asarray(optical_power_w, dtype=float)
    current = responsivity_a_per_w * saturation_power_w * (1 - np.exp(-optical_power / saturation_power_w))
    return _return_scalar_if_scalar_input(current, optical_power_w)


def _validate_saturation_inputs(
    optical_power_w: float | np.ndarray,
    responsivity_a_per_w: float,
    saturation_power_w: float,
    compression_order: float,
) -> None:
    if responsivity_a_per_w <= 0:
        msg = "responsivity_a_per_w must be positive."
        raise ValueError(msg)
    if saturation_power_w <= 0:
        msg = "saturation_power_w must be positive."
        raise ValueError(msg)
    if compression_order <= 0:
        msg = "compression_order must be positive."
        raise ValueError(msg)

    optical_power = np.asarray(optical_power_w, dtype=float)
    if np.any(optical_power < 0):
        msg = "optical_power_w must be non-negative."
        raise ValueError(msg)


def _return_scalar_if_scalar_input(value: np.ndarray, original_input: float | np.ndarray) -> float | np.ndarray:
    if np.isscalar(original_input):
        return float(value)
    return value
