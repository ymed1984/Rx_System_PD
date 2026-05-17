"""Gaussian-level BER calculations for receiver current levels."""

from math import sqrt

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import erfc


def _validate_ordered_levels(mu0_a: float, mu1_a: float) -> None:
    if mu1_a <= mu0_a:
        msg = "mu1_a must be larger than mu0_a."
        raise ValueError(msg)


def _validate_positive_noise(sigma0_a: float, sigma1_a: float) -> None:
    if sigma0_a <= 0:
        msg = "sigma0_a must be positive."
        raise ValueError(msg)
    if sigma1_a <= 0:
        msg = "sigma1_a must be positive."
        raise ValueError(msg)


def qfunc(x: float | np.ndarray) -> float | np.ndarray:
    """Calculate the Gaussian Q function for a scalar or array argument."""
    result = 0.5 * erfc(np.asarray(x) / sqrt(2))
    if np.isscalar(x):
        return float(result)
    return result


def ber_for_threshold(
    mu0_a: float,
    mu1_a: float,
    sigma0_a: float,
    sigma1_a: float,
    threshold_a: float,
) -> float:
    """Calculate BER for Gaussian current levels and a threshold in amperes (A).

    For transmitted 0, an error is current above threshold_a:
    Q((threshold_a - mu0_a) / sigma0_a). For transmitted 1, an error is
    current below threshold_a: Q((mu1_a - threshold_a) / sigma1_a).
    """
    _validate_ordered_levels(mu0_a, mu1_a)
    _validate_positive_noise(sigma0_a, sigma1_a)

    p_error_0 = qfunc((threshold_a - mu0_a) / sigma0_a)
    p_error_1 = qfunc((mu1_a - threshold_a) / sigma1_a)

    return 0.5 * (p_error_0 + p_error_1)


def optimum_threshold(mu0_a: float, mu1_a: float, sigma0_a: float, sigma1_a: float) -> float:
    """Calculate the equal-prior optimum threshold current in amperes (A)."""
    _validate_ordered_levels(mu0_a, mu1_a)
    _validate_positive_noise(sigma0_a, sigma1_a)

    if sigma0_a == sigma1_a:
        return (mu0_a + mu1_a) / 2

    result = minimize_scalar(
        lambda threshold_a: ber_for_threshold(
            mu0_a,
            mu1_a,
            sigma0_a,
            sigma1_a,
            threshold_a,
        ),
        bounds=(mu0_a, mu1_a),
        method="bounded",
    )
    if not result.success:
        msg = "failed to find optimum threshold."
        raise ValueError(msg)
    return float(result.x)


def q_from_levels(mu0_a: float, mu1_a: float, sigma0_a: float, sigma1_a: float) -> float:
    """Calculate the receiver Q estimate from current levels and RMS noise in A."""
    _validate_ordered_levels(mu0_a, mu1_a)
    _validate_positive_noise(sigma0_a, sigma1_a)

    return (mu1_a - mu0_a) / (sigma1_a + sigma0_a)


def ber_from_q(q_rx: float) -> float:
    """Calculate the approximate BER from a receiver Q estimate."""
    if q_rx < 0:
        msg = "q_rx must be non-negative."
        raise ValueError(msg)
    return float(0.5 * erfc(q_rx / sqrt(2)))


def ber_from_gaussian_levels(
    mu0_a: float,
    mu1_a: float,
    sigma0_a: float,
    sigma1_a: float,
    optimize_threshold: bool = True,
) -> dict[str, float]:
    """Calculate BER summary values for Gaussian current levels in amperes (A)."""
    _validate_ordered_levels(mu0_a, mu1_a)
    _validate_positive_noise(sigma0_a, sigma1_a)

    threshold_a = (
        optimum_threshold(mu0_a, mu1_a, sigma0_a, sigma1_a)
        if optimize_threshold
        else (mu0_a + mu1_a) / 2
    )
    q_rx = q_from_levels(mu0_a, mu1_a, sigma0_a, sigma1_a)

    return {
        "mu0_a": mu0_a,
        "mu1_a": mu1_a,
        "sigma0_a": sigma0_a,
        "sigma1_a": sigma1_a,
        "threshold_a": threshold_a,
        "q_rx": q_rx,
        "ber": ber_for_threshold(mu0_a, mu1_a, sigma0_a, sigma1_a, threshold_a),
        "ber_from_q": ber_from_q(q_rx),
    }
