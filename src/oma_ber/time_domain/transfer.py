"""Causal discrete transfer functions for PD and TIA waveform responses."""

from dataclasses import dataclass, field
from math import ceil, isfinite, log, pi, tan
from typing import Literal

import numpy as np
from scipy.signal import lfilter

from oma_ber.time_domain._validation import readonly_float_array


ResponseKind = Literal["dimensionless", "transimpedance_ohm"]


@dataclass(frozen=True)
class DiscreteTransferFunction:
    """Stable causal digital transfer function with an explicit gain unit.

    numerator and denominator are z**-1 polynomial coefficients. sample_rate_hz
    is in samples/s. A dimensionless response maps current to current; a
    transimpedance_ohm response maps current in A to voltage in V.
    """

    numerator: np.ndarray = field(repr=False)
    denominator: np.ndarray = field(repr=False)
    sample_rate_hz: float
    response_kind: ResponseKind

    def __post_init__(self) -> None:
        numerator = readonly_float_array(self.numerator, "numerator")
        denominator = readonly_float_array(self.denominator, "denominator")
        if not isfinite(self.sample_rate_hz) or self.sample_rate_hz <= 0:
            msg = "sample_rate_hz must be finite and positive."
            raise ValueError(msg)
        if self.response_kind not in {"dimensionless", "transimpedance_ohm"}:
            msg = "response_kind must be 'dimensionless' or 'transimpedance_ohm'."
            raise ValueError(msg)
        if denominator[0] == 0:
            msg = "denominator[0] must be non-zero."
            raise ValueError(msg)

        normalized_numerator = np.array(numerator / denominator[0], copy=True)
        normalized_denominator = np.array(denominator / denominator[0], copy=True)
        poles = np.roots(normalized_denominator)
        if poles.size and np.any(np.abs(poles) >= 1):
            msg = "denominator must describe a stable causal filter."
            raise ValueError(msg)
        denominator_dc = float(np.sum(normalized_denominator))
        denominator_scale = float(np.sum(np.abs(normalized_denominator)))
        if abs(denominator_dc) <= np.finfo(float).eps * denominator_scale:
            msg = "transfer function must have a finite DC gain."
            raise ValueError(msg)

        normalized_numerator.setflags(write=False)
        normalized_denominator.setflags(write=False)
        object.__setattr__(self, "numerator", normalized_numerator)
        object.__setattr__(self, "denominator", normalized_denominator)

    @property
    def dc_gain(self) -> float:
        """Return the signed DC gain, dimensionless or in ohms (Ω)."""
        return float(np.sum(self.numerator) / np.sum(self.denominator))

    def settling_samples(self, tolerance: float = 1e-9) -> int:
        """Estimate samples needed for zero-state transient decay.

        tolerance is the requested relative pole-envelope remainder. FIR memory
        is included exactly; for stable IIR filters the dominant pole radius is
        used as a conservative first-order estimate.
        """
        if not isfinite(tolerance) or not 0 < tolerance < 1:
            msg = "tolerance must be finite and between 0 and 1."
            raise ValueError(msg)

        finite_memory = max(self.numerator.size, self.denominator.size) - 1
        poles = np.roots(self.denominator)
        if poles.size == 0:
            return finite_memory
        dominant_radius = float(np.max(np.abs(poles)))
        if dominant_radius == 0:
            return finite_memory
        return max(finite_memory, ceil(log(tolerance) / log(dominant_radius)))


def identity_transfer(
    sample_rate_hz: float,
    dc_gain: float = 1.0,
    response_kind: ResponseKind = "dimensionless",
) -> DiscreteTransferFunction:
    """Create a causal flat transfer with explicit DC gain."""
    if not isfinite(dc_gain) or dc_gain == 0:
        msg = "dc_gain must be finite and non-zero."
        raise ValueError(msg)
    return DiscreteTransferFunction(
        numerator=np.array([dc_gain]),
        denominator=np.array([1.0]),
        sample_rate_hz=sample_rate_hz,
        response_kind=response_kind,
    )


def first_order_lowpass_transfer(
    sample_rate_hz: float,
    bandwidth_3db_hz: float,
    dc_gain: float = 1.0,
    response_kind: ResponseKind = "dimensionless",
) -> DiscreteTransferFunction:
    """Create a causal one-pole low-pass with an exact digital 3 dB point.

    The analog one-pole response is mapped with a prewarped bilinear transform.
    sample_rate_hz and bandwidth_3db_hz are in Hz. dc_gain is dimensionless for
    a PD response and in ohms for a transimpedance response.
    """
    if not isfinite(sample_rate_hz) or sample_rate_hz <= 0:
        msg = "sample_rate_hz must be finite and positive."
        raise ValueError(msg)
    if not isfinite(bandwidth_3db_hz) or bandwidth_3db_hz <= 0:
        msg = "bandwidth_3db_hz must be finite and positive."
        raise ValueError(msg)
    if bandwidth_3db_hz >= sample_rate_hz / 2:
        msg = "bandwidth_3db_hz must be below the Nyquist frequency."
        raise ValueError(msg)
    if not isfinite(dc_gain) or dc_gain == 0:
        msg = "dc_gain must be finite and non-zero."
        raise ValueError(msg)

    prewarped = tan(pi * bandwidth_3db_hz / sample_rate_hz)
    b0 = dc_gain * prewarped / (1 + prewarped)
    a1 = (prewarped - 1) / (prewarped + 1)
    return DiscreteTransferFunction(
        numerator=np.array([b0, b0]),
        denominator=np.array([1.0, a1]),
        sample_rate_hz=sample_rate_hz,
        response_kind=response_kind,
    )


def apply_transfer(
    waveform: np.ndarray,
    transfer: DiscreteTransferFunction,
) -> np.ndarray:
    """Apply a causal transfer function to a finite waveform."""
    values = readonly_float_array(waveform, "waveform")
    return np.asarray(
        lfilter(transfer.numerator, transfer.denominator, values), dtype=float
    )


def transfer_impulse_response(
    transfer: DiscreteTransferFunction,
    tolerance: float = 1e-12,
) -> np.ndarray:
    """Return a truncated causal impulse response for a digital transfer.

    The response contains the complete FIR memory. Stable IIR tails are kept
    until the dominant-pole envelope is below tolerance. Coefficients retain
    the transfer gain unit: dimensionless for PD responses and ohms for TIA
    transimpedance responses.
    """
    num_samples = transfer.settling_samples(tolerance) + 1
    impulse = np.zeros(num_samples, dtype=float)
    impulse[0] = 1.0
    response = np.asarray(
        lfilter(transfer.numerator, transfer.denominator, impulse),
        dtype=float,
    )
    response.setflags(write=False)
    return response


def equivalent_noise_bandwidth_hz(
    transfer: DiscreteTransferFunction,
    tolerance: float = 1e-12,
) -> float:
    """Calculate one-sided equivalent noise bandwidth in hertz (Hz).

    For white input noise and sampled causal impulse response h[n],
    Bn = fs/2 * sum(h[n]**2) / abs(H(0))**2. The transfer must have non-zero
    DC gain.
    """
    dc_gain = transfer.dc_gain
    if dc_gain == 0:
        msg = "equivalent noise bandwidth requires a non-zero DC gain."
        raise ValueError(msg)
    impulse_response = transfer_impulse_response(transfer, tolerance)
    return float(
        transfer.sample_rate_hz / 2 * np.sum(impulse_response**2) / abs(dc_gain) ** 2
    )
