"""Scalar RMS noise calculations for photodiode receiver levels."""

from math import sqrt

from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.units import rin_db_per_hz_to_linear

Q_E = 1.602176634e-19


def shot_noise_rms_a(current_a: float, dark_current_a: float, bandwidth_hz: float) -> float:
    """Calculate level-dependent shot noise RMS current in amperes (A).

    current_a and dark_current_a are currents in A. bandwidth_hz is the
    noise bandwidth in Hz.
    """
    if current_a < 0:
        msg = "current_a must be non-negative."
        raise ValueError(msg)
    if dark_current_a < 0:
        msg = "dark_current_a must be non-negative."
        raise ValueError(msg)
    if bandwidth_hz <= 0:
        msg = "bandwidth_hz must be positive."
        raise ValueError(msg)

    return sqrt(2 * Q_E * (current_a + dark_current_a) * bandwidth_hz)


def tia_noise_rms_a(
    input_noise_density_a_per_sqrt_hz: float,
    bandwidth_hz: float,
) -> float:
    """Calculate TIA input-referred RMS current noise in amperes (A).

    input_noise_density_a_per_sqrt_hz is in A/sqrt(Hz). bandwidth_hz is in Hz.
    """
    if input_noise_density_a_per_sqrt_hz < 0:
        msg = "input_noise_density_a_per_sqrt_hz must be non-negative."
        raise ValueError(msg)
    if bandwidth_hz <= 0:
        msg = "bandwidth_hz must be positive."
        raise ValueError(msg)

    return input_noise_density_a_per_sqrt_hz * sqrt(bandwidth_hz)


def rin_noise_rms_a(
    responsivity_a_per_w: float,
    optical_power_w: float,
    rin_db_per_hz: float | None,
    bandwidth_hz: float,
) -> float:
    """Calculate RIN RMS current noise in amperes (A).

    responsivity_a_per_w is in A/W, optical_power_w is in W,
    rin_db_per_hz is in dB/Hz when provided, and bandwidth_hz is in Hz.
    """
    if responsivity_a_per_w <= 0:
        msg = "responsivity_a_per_w must be positive."
        raise ValueError(msg)
    if optical_power_w < 0:
        msg = "optical_power_w must be non-negative."
        raise ValueError(msg)
    if bandwidth_hz <= 0:
        msg = "bandwidth_hz must be positive."
        raise ValueError(msg)
    if rin_db_per_hz is None:
        return 0.0

    rin_linear_per_hz = rin_db_per_hz_to_linear(rin_db_per_hz)
    return responsivity_a_per_w * optical_power_w * sqrt(rin_linear_per_hz * bandwidth_hz)


def total_noise_rms_a(
    optical_power_w: float,
    photocurrent_a: float,
    pd: Photodiode,
    rx: Receiver,
) -> float:
    """Calculate total RMS current noise in amperes (A).

    optical_power_w is the optical level in W and photocurrent_a is the
    corresponding signal current in A.
    """
    sigma_shot_a = shot_noise_rms_a(
        current_a=photocurrent_a,
        dark_current_a=pd.dark_current_a,
        bandwidth_hz=rx.noise_bandwidth_hz,
    )
    sigma_tia_a = tia_noise_rms_a(
        input_noise_density_a_per_sqrt_hz=rx.input_current_noise_density_a_per_sqrt_hz,
        bandwidth_hz=rx.noise_bandwidth_hz,
    )
    sigma_rin_a = rin_noise_rms_a(
        responsivity_a_per_w=pd.responsivity_a_per_w,
        optical_power_w=optical_power_w,
        rin_db_per_hz=rx.rin_db_per_hz,
        bandwidth_hz=rx.noise_bandwidth_hz,
    )

    return sqrt(sigma_shot_a**2 + sigma_tia_a**2 + sigma_rin_a**2)
