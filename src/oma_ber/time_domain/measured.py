"""Measured complex-response import and causal FIR conversion."""

from dataclasses import dataclass, field
from math import isfinite, pi
from pathlib import Path

import numpy as np

from oma_ber.time_domain.transfer import (
    DiscreteTransferFunction,
    ResponseKind,
    equivalent_noise_bandwidth_hz,
)


@dataclass(frozen=True)
class MeasuredFrequencyResponse:
    """Complex amplitude response at an explicit physical reference plane.

    frequency_hz contains non-negative frequencies in Hz. response_complex is
    a linear amplitude response: A/A for ``dimensionless`` or V/A (ohms) for
    ``transimpedance_ohm``. The response must include DC explicitly.
    """

    frequency_hz: np.ndarray = field(repr=False)
    response_complex: np.ndarray = field(repr=False)
    response_kind: ResponseKind
    reference_plane: str
    source_path: Path | None = None

    def __post_init__(self) -> None:
        frequency_hz = np.asarray(self.frequency_hz, dtype=float)
        response_complex = np.asarray(self.response_complex, dtype=complex)
        if frequency_hz.ndim != 1:
            msg = "frequency_hz must be a one-dimensional array."
            raise ValueError(msg)
        if response_complex.ndim != 1:
            msg = "response_complex must be a one-dimensional array."
            raise ValueError(msg)
        if frequency_hz.size < 2:
            msg = "frequency_hz must contain at least two frequencies."
            raise ValueError(msg)
        if frequency_hz.size != response_complex.size:
            msg = "response_complex length must match frequency_hz length."
            raise ValueError(msg)
        if not np.all(np.isfinite(frequency_hz)):
            msg = "frequency_hz must contain only finite values."
            raise ValueError(msg)
        if not np.all(
            np.isfinite(response_complex.real) & np.isfinite(response_complex.imag)
        ):
            msg = "response_complex must contain only finite values."
            raise ValueError(msg)
        if frequency_hz[0] != 0.0:
            msg = "frequency_hz must include DC as its first value (0 Hz)."
            raise ValueError(msg)
        if not np.all(np.diff(frequency_hz) > 0):
            msg = "frequency_hz values must be strictly increasing."
            raise ValueError(msg)
        if self.response_kind not in {"dimensionless", "transimpedance_ohm"}:
            msg = "response_kind must be 'dimensionless' or 'transimpedance_ohm'."
            raise ValueError(msg)
        if (
            not isinstance(self.reference_plane, str)
            or not self.reference_plane.strip()
        ):
            msg = "reference_plane must be a non-empty string."
            raise ValueError(msg)
        if self.source_path is not None and not isinstance(self.source_path, Path):
            msg = "source_path must be a pathlib.Path when provided."
            raise ValueError(msg)

        dc_scale = max(abs(response_complex[0]), 1.0)
        if abs(response_complex[0].imag) > 1e-12 * dc_scale:
            msg = "the DC response must be real for a real impulse response."
            raise ValueError(msg)
        if abs(response_complex[0].real) <= np.finfo(float).tiny:
            msg = "the DC response must be non-zero."
            raise ValueError(msg)
        if np.any(np.abs(response_complex[:-1]) == 0):
            msg = (
                "response magnitude may be zero only at the final frequency; "
                "interior zeros have undefined phase for interpolation."
            )
            raise ValueError(msg)

        copied_frequency = np.array(frequency_hz, copy=True)
        copied_response = np.array(response_complex, copy=True)
        copied_frequency.setflags(write=False)
        copied_response.setflags(write=False)
        object.__setattr__(self, "frequency_hz", copied_frequency)
        object.__setattr__(self, "response_complex", copied_response)
        object.__setattr__(self, "reference_plane", self.reference_plane.strip())


@dataclass(frozen=True)
class MeasuredResponseDiagnostics:
    """Diagnostics for a measured response converted to a causal FIR."""

    measured_dc_gain: float
    fir_dc_gain: float
    dc_gain_relative_error: float
    equivalent_noise_bandwidth_hz: float
    group_delay_s: float
    pre_echo_energy_ratio: float
    discarded_tail_energy_ratio: float
    removed_reference_delay_s: float
    sample_rate_hz: float
    fft_size: int
    fir_length_samples: int

    def __post_init__(self) -> None:
        finite_values = (
            self.measured_dc_gain,
            self.fir_dc_gain,
            self.dc_gain_relative_error,
            self.equivalent_noise_bandwidth_hz,
            self.group_delay_s,
            self.pre_echo_energy_ratio,
            self.discarded_tail_energy_ratio,
            self.removed_reference_delay_s,
            self.sample_rate_hz,
        )
        if not all(isfinite(value) for value in finite_values):
            msg = "measured-response diagnostics must contain finite values."
            raise ValueError(msg)
        if self.measured_dc_gain == 0 or self.fir_dc_gain == 0:
            msg = "measured and FIR DC gains must be non-zero."
            raise ValueError(msg)
        if self.dc_gain_relative_error < 0:
            msg = "dc_gain_relative_error must be non-negative."
            raise ValueError(msg)
        if self.equivalent_noise_bandwidth_hz <= 0:
            msg = "equivalent_noise_bandwidth_hz must be positive."
            raise ValueError(msg)
        if self.pre_echo_energy_ratio < 0 or self.discarded_tail_energy_ratio < 0:
            msg = "response energy ratios must be non-negative."
            raise ValueError(msg)
        if self.removed_reference_delay_s < 0:
            msg = "removed_reference_delay_s must be non-negative."
            raise ValueError(msg)
        if self.sample_rate_hz <= 0:
            msg = "sample_rate_hz must be positive."
            raise ValueError(msg)
        if self.fft_size < 4 or self.fft_size % 2:
            msg = "fft_size must be an even integer of at least 4."
            raise ValueError(msg)
        if not 1 <= self.fir_length_samples <= self.fft_size // 2:
            msg = "fir_length_samples must be in [1, fft_size/2]."
            raise ValueError(msg)


def measured_frequency_response_from_csv(
    path: str | Path,
    response_kind: ResponseKind,
    reference_plane: str,
) -> MeasuredFrequencyResponse:
    """Load a physical complex amplitude response from CSV.

    The CSV columns are ``frequency_hz``, ``magnitude_db``, and ``phase_deg``.
    magnitude_db is always an amplitude quantity using 20*log10: A/A for a
    dimensionless PD response or V/A for a TIA transimpedance response.
    """
    source_path = Path(path)
    data = np.genfromtxt(source_path, delimiter=",", names=True)
    if data.dtype.names is None:
        msg = "CSV must contain a header row."
        raise ValueError(msg)
    required_columns = {"frequency_hz", "magnitude_db", "phase_deg"}
    missing_columns = required_columns - set(data.dtype.names)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        msg = f"CSV is missing required columns: {missing}."
        raise ValueError(msg)

    frequency_hz = np.atleast_1d(np.asarray(data["frequency_hz"], dtype=float))
    magnitude_db = np.atleast_1d(np.asarray(data["magnitude_db"], dtype=float))
    phase_deg = np.atleast_1d(np.asarray(data["phase_deg"], dtype=float))
    if not np.all(np.isfinite(magnitude_db)):
        msg = "magnitude_db must contain only finite values."
        raise ValueError(msg)
    if not np.all(np.isfinite(phase_deg)):
        msg = "phase_deg must contain only finite values."
        raise ValueError(msg)

    response_complex = 10 ** (magnitude_db / 20) * np.exp(1j * np.deg2rad(phase_deg))
    return MeasuredFrequencyResponse(
        frequency_hz=frequency_hz,
        response_complex=response_complex,
        response_kind=response_kind,
        reference_plane=reference_plane,
        source_path=source_path,
    )


def fir_transfer_from_measured_response(
    response: MeasuredFrequencyResponse,
    sample_rate_hz: float,
    fft_size: int,
    fir_length_samples: int,
    removed_reference_delay_s: float = 0.0,
    causality_tolerance: float = 1e-6,
    tail_energy_tolerance: float = 1e-6,
) -> tuple[DiscreteTransferFunction, MeasuredResponseDiagnostics]:
    """Convert a measured complex response to a validated causal FIR.

    The input must explicitly cover DC through Nyquist. Magnitude and unwrapped
    phase are linearly interpolated onto the rFFT grid. Only the caller-supplied
    pure reference delay is removed. Responses with excessive periodic
    negative-time energy or discarded positive-time tail energy are rejected
    rather than silently shifted or truncated.
    """
    if not isinstance(response, MeasuredFrequencyResponse):
        msg = "response must be a MeasuredFrequencyResponse."
        raise ValueError(msg)
    _validate_conversion_inputs(
        sample_rate_hz=sample_rate_hz,
        fft_size=fft_size,
        fir_length_samples=fir_length_samples,
        removed_reference_delay_s=removed_reference_delay_s,
        causality_tolerance=causality_tolerance,
        tail_energy_tolerance=tail_energy_tolerance,
    )
    nyquist_hz = sample_rate_hz / 2
    coverage_tolerance_hz = (
        16
        * np.finfo(float).eps
        * max(
            nyquist_hz,
            response.frequency_hz[-1],
        )
    )
    if response.frequency_hz[-1] < nyquist_hz - coverage_tolerance_hz:
        msg = (
            "measured frequency response must cover the Nyquist frequency; "
            "high-frequency extrapolation is not implicit."
        )
        raise ValueError(msg)

    target_frequency_hz = np.fft.rfftfreq(fft_size, d=1 / sample_rate_hz)
    source_magnitude = np.abs(response.response_complex)
    source_phase_rad = _phase_for_interpolation(response.response_complex)
    deembedded_phase_rad = source_phase_rad + (
        2 * pi * response.frequency_hz * removed_reference_delay_s
    )
    target_magnitude = np.interp(
        target_frequency_hz,
        response.frequency_hz,
        source_magnitude,
    )
    target_phase_rad = np.interp(
        target_frequency_hz,
        response.frequency_hz,
        deembedded_phase_rad,
    )
    target_response = target_magnitude * np.exp(1j * target_phase_rad)
    _validate_real_fft_endpoints(target_response)
    target_response[0] = target_response[0].real
    target_response[-1] = target_response[-1].real

    periodic_impulse = np.asarray(
        np.fft.irfft(target_response, n=fft_size), dtype=float
    )
    total_energy = float(np.sum(periodic_impulse**2))
    if total_energy <= np.finfo(float).tiny:
        msg = "measured response produced zero impulse-response energy."
        raise ValueError(msg)
    positive_limit = fft_size // 2
    pre_echo_energy = float(np.sum(periodic_impulse[positive_limit + 1 :] ** 2))
    discarded_tail_energy = float(
        np.sum(periodic_impulse[fir_length_samples : positive_limit + 1] ** 2)
    )
    pre_echo_energy_ratio = pre_echo_energy / total_energy
    discarded_tail_energy_ratio = discarded_tail_energy / total_energy
    if pre_echo_energy_ratio > causality_tolerance:
        msg = (
            "measured response exceeds causality_tolerance: "
            f"pre_echo_energy_ratio={pre_echo_energy_ratio:.6e}."
        )
        raise ValueError(msg)
    if discarded_tail_energy_ratio > tail_energy_tolerance:
        msg = (
            "measured response exceeds tail_energy_tolerance: "
            f"discarded_tail_energy_ratio={discarded_tail_energy_ratio:.6e}."
        )
        raise ValueError(msg)

    fir = np.array(periodic_impulse[:fir_length_samples], copy=True)
    transfer = DiscreteTransferFunction(
        numerator=fir,
        denominator=np.array([1.0]),
        sample_rate_hz=sample_rate_hz,
        response_kind=response.response_kind,
    )
    measured_dc_gain = float(response.response_complex[0].real)
    fir_dc_gain = transfer.dc_gain
    diagnostics = MeasuredResponseDiagnostics(
        measured_dc_gain=measured_dc_gain,
        fir_dc_gain=fir_dc_gain,
        dc_gain_relative_error=abs(fir_dc_gain - measured_dc_gain)
        / abs(measured_dc_gain),
        equivalent_noise_bandwidth_hz=equivalent_noise_bandwidth_hz(transfer),
        group_delay_s=_low_frequency_group_delay_s(
            target_frequency_hz,
            target_response,
        ),
        pre_echo_energy_ratio=pre_echo_energy_ratio,
        discarded_tail_energy_ratio=discarded_tail_energy_ratio,
        removed_reference_delay_s=removed_reference_delay_s,
        sample_rate_hz=sample_rate_hz,
        fft_size=fft_size,
        fir_length_samples=fir_length_samples,
    )
    return transfer, diagnostics


def _validate_conversion_inputs(
    sample_rate_hz: float,
    fft_size: int,
    fir_length_samples: int,
    removed_reference_delay_s: float,
    causality_tolerance: float,
    tail_energy_tolerance: float,
) -> None:
    if not isfinite(sample_rate_hz) or sample_rate_hz <= 0:
        msg = "sample_rate_hz must be finite and positive."
        raise ValueError(msg)
    if isinstance(fft_size, bool) or not isinstance(fft_size, int):
        msg = "fft_size must be an integer."
        raise ValueError(msg)
    if fft_size < 4 or fft_size % 2:
        msg = "fft_size must be an even integer of at least 4."
        raise ValueError(msg)
    if isinstance(fir_length_samples, bool) or not isinstance(fir_length_samples, int):
        msg = "fir_length_samples must be an integer."
        raise ValueError(msg)
    if not 1 <= fir_length_samples <= fft_size // 2:
        msg = "fir_length_samples must be in [1, fft_size/2]."
        raise ValueError(msg)
    if not isfinite(removed_reference_delay_s) or removed_reference_delay_s < 0:
        msg = "removed_reference_delay_s must be finite and non-negative."
        raise ValueError(msg)
    for value, name in (
        (causality_tolerance, "causality_tolerance"),
        (tail_energy_tolerance, "tail_energy_tolerance"),
    ):
        if not isfinite(value) or not 0 <= value < 1:
            msg = f"{name} must be finite and in [0, 1)."
            raise ValueError(msg)


def _phase_for_interpolation(response_complex: np.ndarray) -> np.ndarray:
    phase_rad = np.angle(response_complex)
    if abs(response_complex[-1]) == 0:
        phase_rad[-1] = phase_rad[-2]
    return np.unwrap(phase_rad)


def _validate_real_fft_endpoints(target_response: np.ndarray) -> None:
    for value, name in (
        (target_response[0], "DC"),
        (target_response[-1], "Nyquist"),
    ):
        scale = max(abs(value), 1.0)
        if abs(value.imag) > 1e-10 * scale:
            msg = (
                f"interpolated {name} response must be real for a real impulse "
                "response."
            )
            raise ValueError(msg)


def _low_frequency_group_delay_s(
    frequency_hz: np.ndarray,
    response_complex: np.ndarray,
) -> float:
    if frequency_hz.size < 2:
        return 0.0
    phase_rad = np.unwrap(np.angle(response_complex))
    return float(
        -(phase_rad[1] - phase_rad[0]) / (2 * pi * (frequency_hz[1] - frequency_hz[0]))
    )
