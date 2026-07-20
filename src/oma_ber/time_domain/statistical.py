"""Gaussian-mixture statistical eye and BER analysis for TIA voltage."""

from dataclasses import dataclass, field
from math import isfinite, sqrt

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import log_ndtr, logsumexp

from oma_ber.time_domain._validation import readonly_float_array
from oma_ber.time_domain.noise import TiaOutputNoiseResult
from oma_ber.time_domain.receiver import ReceiverWaveformResult


@dataclass(frozen=True)
class StatisticalEyePhaseResult:
    """Pattern-mixture voltage statistics and optimized BER at one phase."""

    sampling_phase_ui: float
    decision_delay_s: float
    one_is_high: bool
    threshold_v: float
    ber: float
    effective_q: float
    mean_zero_v: float
    mean_one_v: float
    rms_zero_noise_v: float
    rms_one_noise_v: float
    zero_means_v: np.ndarray = field(repr=False)
    one_means_v: np.ndarray = field(repr=False)
    zero_sigmas_v: np.ndarray = field(repr=False)
    one_sigmas_v: np.ndarray = field(repr=False)

    def __post_init__(self) -> None:
        if not isfinite(self.sampling_phase_ui) or not 0 <= self.sampling_phase_ui < 1:
            msg = "sampling_phase_ui must be finite and in [0, 1)."
            raise ValueError(msg)
        if not isfinite(self.decision_delay_s):
            msg = "decision_delay_s must be finite."
            raise ValueError(msg)
        if not isinstance(self.one_is_high, bool):
            msg = "one_is_high must be a bool."
            raise ValueError(msg)
        if not isfinite(self.threshold_v):
            msg = "threshold_v must be finite."
            raise ValueError(msg)
        if not isfinite(self.ber) or not 0 <= self.ber <= 0.5:
            msg = "ber must be finite and in [0, 0.5]."
            raise ValueError(msg)
        if not isfinite(self.effective_q) or self.effective_q < 0:
            msg = "effective_q must be finite and non-negative."
            raise ValueError(msg)
        scalar_voltage_values = (
            self.mean_zero_v,
            self.mean_one_v,
            self.rms_zero_noise_v,
            self.rms_one_noise_v,
        )
        if not all(isfinite(value) for value in scalar_voltage_values):
            msg = "statistical eye voltage summaries must be finite."
            raise ValueError(msg)
        if self.rms_zero_noise_v <= 0 or self.rms_one_noise_v <= 0:
            msg = "statistical eye RMS noise values must be positive."
            raise ValueError(msg)

        pairs = (
            ("zero_means_v", "zero_sigmas_v"),
            ("one_means_v", "one_sigmas_v"),
        )
        for mean_name, sigma_name in pairs:
            means = readonly_float_array(getattr(self, mean_name), mean_name)
            sigmas = readonly_float_array(getattr(self, sigma_name), sigma_name)
            if means.size != sigmas.size:
                msg = f"{mean_name} and {sigma_name} must have the same length."
                raise ValueError(msg)
            if np.any(sigmas <= 0):
                msg = f"{sigma_name} must be positive."
                raise ValueError(msg)
            object.__setattr__(self, mean_name, means)
            object.__setattr__(self, sigma_name, sigmas)


@dataclass(frozen=True)
class StatisticalEyeAnalysis:
    """Phase-dependent Gaussian-mixture BER and statistical eye density."""

    phase_results: tuple[StatisticalEyePhaseResult, ...]
    optimum_result: StatisticalEyePhaseResult
    sampling_phases_ui: np.ndarray = field(repr=False)
    amplitude_axis_v: np.ndarray = field(repr=False)
    density_per_v: np.ndarray = field(repr=False)

    def __post_init__(self) -> None:
        phases = readonly_float_array(self.sampling_phases_ui, "sampling_phases_ui")
        amplitude = readonly_float_array(self.amplitude_axis_v, "amplitude_axis_v")
        density = np.asarray(self.density_per_v, dtype=float)
        if (
            np.any((phases < 0) | (phases >= 1))
            or np.unique(phases).size != phases.size
        ):
            msg = "sampling_phases_ui must contain unique values in [0, 1)."
            raise ValueError(msg)
        if not np.all(np.diff(amplitude) > 0):
            msg = "amplitude_axis_v must be strictly increasing."
            raise ValueError(msg)
        if len(self.phase_results) != phases.size:
            msg = "phase_results length must match sampling_phases_ui."
            raise ValueError(msg)
        if density.shape != (phases.size, amplitude.size):
            msg = "density_per_v shape must be (num_phases, num_amplitude_points)."
            raise ValueError(msg)
        if not np.all(np.isfinite(density)) or np.any(density < 0):
            msg = "density_per_v must contain finite non-negative values."
            raise ValueError(msg)
        if not any(result is self.optimum_result for result in self.phase_results):
            msg = "optimum_result must be one of phase_results."
            raise ValueError(msg)
        copied_density = np.array(density, copy=True)
        copied_density.setflags(write=False)
        object.__setattr__(self, "sampling_phases_ui", phases)
        object.__setattr__(self, "amplitude_axis_v", amplitude)
        object.__setattr__(self, "density_per_v", copied_density)


def gaussian_mixture_ber_for_threshold(
    zero_means_v: np.ndarray,
    one_means_v: np.ndarray,
    zero_sigmas_v: np.ndarray,
    one_sigmas_v: np.ndarray,
    threshold_v: float,
    one_is_high: bool,
    zero_weights: np.ndarray | None = None,
    one_weights: np.ndarray | None = None,
) -> float:
    """Calculate equal-prior BER for pattern mixtures at a voltage threshold.

    Every pattern-conditioned component is Gaussian. Bit 0 and bit 1 receive
    equal prior probability 0.5, independent of the finite PRBS balance.
    Components within each bit value are uniform unless normalized positive
    zero_weights and one_weights are provided.
    """
    zero_means, one_means, zero_sigmas, one_sigmas = _validated_mixture_arrays(
        zero_means_v,
        one_means_v,
        zero_sigmas_v,
        one_sigmas_v,
    )
    if not isfinite(threshold_v):
        msg = "threshold_v must be finite."
        raise ValueError(msg)
    if not isinstance(one_is_high, bool):
        msg = "one_is_high must be a bool."
        raise ValueError(msg)
    validated_zero_weights = _validated_mixture_weights(
        zero_weights,
        zero_means.size,
        "zero_weights",
    )
    validated_one_weights = _validated_mixture_weights(
        one_weights,
        one_means.size,
        "one_weights",
    )

    polarity = 1.0 if one_is_high else -1.0
    threshold = polarity * threshold_v
    transformed_zero = polarity * zero_means
    transformed_one = polarity * one_means
    log_ber = _transformed_mixture_log_ber_array(
        np.array([threshold]),
        transformed_zero,
        transformed_one,
        zero_sigmas,
        one_sigmas,
        validated_zero_weights,
        validated_one_weights,
    )[0]
    return float(np.exp(log_ber))


def optimize_gaussian_mixture_threshold(
    zero_means_v: np.ndarray,
    one_means_v: np.ndarray,
    zero_sigmas_v: np.ndarray,
    one_sigmas_v: np.ndarray,
    one_is_high: bool,
    grid_points: int = 1025,
    zero_weights: np.ndarray | None = None,
    one_weights: np.ndarray | None = None,
) -> tuple[float, float]:
    """Find a robust global-grid-plus-local optimum voltage threshold."""
    zero_means, one_means, zero_sigmas, one_sigmas = _validated_mixture_arrays(
        zero_means_v,
        one_means_v,
        zero_sigmas_v,
        one_sigmas_v,
    )
    if not isinstance(one_is_high, bool):
        msg = "one_is_high must be a bool."
        raise ValueError(msg)
    if isinstance(grid_points, bool) or not isinstance(grid_points, int):
        msg = "grid_points must be an integer."
        raise ValueError(msg)
    if grid_points < 33:
        msg = "grid_points must be at least 33."
        raise ValueError(msg)
    validated_zero_weights = _validated_mixture_weights(
        zero_weights,
        zero_means.size,
        "zero_weights",
    )
    validated_one_weights = _validated_mixture_weights(
        one_weights,
        one_means.size,
        "one_weights",
    )

    polarity = 1.0 if one_is_high else -1.0
    transformed_zero = polarity * zero_means
    transformed_one = polarity * one_means
    maximum_sigma = float(max(np.max(zero_sigmas), np.max(one_sigmas)))
    lower = float(
        min(np.min(transformed_zero), np.min(transformed_one)) - 8 * maximum_sigma
    )
    upper = float(
        max(np.max(transformed_zero), np.max(transformed_one)) + 8 * maximum_sigma
    )
    threshold_grid = np.linspace(lower, upper, grid_points)
    log_ber_grid = _transformed_mixture_log_ber_array(
        threshold_grid,
        transformed_zero,
        transformed_one,
        zero_sigmas,
        one_sigmas,
        validated_zero_weights,
        validated_one_weights,
    )
    best_index = int(np.argmin(log_ber_grid))

    left_index = max(0, best_index - 1)
    right_index = min(grid_points - 1, best_index + 1)
    local_lower = float(threshold_grid[left_index])
    local_upper = float(threshold_grid[right_index])
    if local_lower == local_upper:
        optimum_transformed = float(threshold_grid[best_index])
    else:
        scale = max(
            abs(local_lower), abs(local_upper), maximum_sigma, np.finfo(float).tiny
        )
        optimum = minimize_scalar(
            lambda threshold: float(
                _transformed_mixture_log_ber_array(
                    np.array([threshold]),
                    transformed_zero,
                    transformed_one,
                    zero_sigmas,
                    one_sigmas,
                    validated_zero_weights,
                    validated_one_weights,
                )[0]
            ),
            bounds=(local_lower, local_upper),
            method="bounded",
            options={"xatol": 64 * np.finfo(float).eps * scale},
        )
        optimum_transformed = float(optimum.x)

    threshold_v = polarity * optimum_transformed
    ber = gaussian_mixture_ber_for_threshold(
        zero_means,
        one_means,
        zero_sigmas,
        one_sigmas,
        threshold_v,
        one_is_high,
        zero_weights=validated_zero_weights,
        one_weights=validated_one_weights,
    )
    return threshold_v, ber


def analyze_statistical_tia_eye(
    waveforms: ReceiverWaveformResult,
    noise: TiaOutputNoiseResult,
    decision_delay_s: float = 0.0,
    sampling_phases_ui: np.ndarray | None = None,
    amplitude_points: int = 401,
    density_sigma_span: float = 7.0,
) -> StatisticalEyeAnalysis:
    """Calculate statistical TIA-voltage eye density and BER versus phase.

    The deterministic pattern means come from the causal Phase 1 TIA waveform.
    Pattern-dependent variances and adjacent-sample covariance come from
    calculate_tia_output_noise(). Fractional sampling propagates both variance
    and lag-one covariance through the linear interpolation operation.
    """
    voltage_v = waveforms.tia_output_voltage_v
    if voltage_v is None:
        msg = "statistical TIA eye requires a TIA voltage waveform."
        raise ValueError(msg)
    if noise.total_variance_v2.size != voltage_v.size:
        msg = "noise length must match the TIA voltage waveform length."
        raise ValueError(msg)
    if not isfinite(decision_delay_s):
        msg = "decision_delay_s must be finite."
        raise ValueError(msg)
    if isinstance(amplitude_points, bool) or not isinstance(amplitude_points, int):
        msg = "amplitude_points must be an integer."
        raise ValueError(msg)
    if amplitude_points < 101:
        msg = "amplitude_points must be at least 101."
        raise ValueError(msg)
    if not isfinite(density_sigma_span) or density_sigma_span <= 0:
        msg = "density_sigma_span must be finite and positive."
        raise ValueError(msg)

    if sampling_phases_ui is None:
        phases = np.arange(waveforms.time_grid.samples_per_symbol, dtype=float) / (
            waveforms.time_grid.samples_per_symbol
        )
    else:
        phases = readonly_float_array(sampling_phases_ui, "sampling_phases_ui")
        if np.any((phases < 0) | (phases >= 1)):
            msg = "sampling_phases_ui values must be in [0, 1)."
            raise ValueError(msg)
        if np.unique(phases).size != phases.size:
            msg = "sampling_phases_ui values must be unique."
            raise ValueError(msg)

    phase_results = tuple(
        _analyze_phase(
            waveforms=waveforms,
            noise=noise,
            sampling_phase_ui=float(phase),
            decision_delay_s=decision_delay_s,
        )
        for phase in phases
    )
    optimum_result = min(
        phase_results,
        key=lambda result: (
            result.ber,
            -result.effective_q,
            abs(result.sampling_phase_ui - 0.5),
        ),
    )
    amplitude_axis_v = _amplitude_axis(
        phase_results,
        amplitude_points,
        density_sigma_span,
    )
    density_per_v = np.vstack(
        [_phase_density(result, amplitude_axis_v) for result in phase_results]
    )
    return StatisticalEyeAnalysis(
        phase_results=phase_results,
        optimum_result=optimum_result,
        sampling_phases_ui=phases,
        amplitude_axis_v=amplitude_axis_v,
        density_per_v=density_per_v,
    )


def _analyze_phase(
    waveforms: ReceiverWaveformResult,
    noise: TiaOutputNoiseResult,
    sampling_phase_ui: float,
    decision_delay_s: float,
) -> StatisticalEyePhaseResult:
    voltage_v = waveforms.tia_output_voltage_v
    if voltage_v is None:
        msg = "statistical TIA eye requires a TIA voltage waveform."
        raise ValueError(msg)
    positions, sampled_bits = _sample_positions(
        num_symbols=waveforms.bits.size,
        num_samples=voltage_v.size,
        samples_per_symbol=waveforms.time_grid.samples_per_symbol,
        sample_rate_hz=waveforms.time_grid.sample_rate_hz,
        sampling_phase_ui=sampling_phase_ui,
        decision_delay_s=decision_delay_s,
        bits=waveforms.bits,
    )
    voltage_v = np.asarray(voltage_v)
    sampled_means_v = np.interp(
        positions,
        np.arange(voltage_v.size, dtype=float),
        voltage_v,
    )
    sampled_variances_v2 = _interpolated_variance(
        noise.total_variance_v2,
        noise.total_lag1_covariance_v2,
        positions,
    )
    sampled_sigmas_v = np.sqrt(sampled_variances_v2)
    if np.any(sampled_sigmas_v <= 0):
        msg = (
            "statistical Gaussian eye requires positive noise at every sampled "
            "pattern component."
        )
        raise ValueError(msg)

    zero_means_v = sampled_means_v[sampled_bits == 0]
    one_means_v = sampled_means_v[sampled_bits == 1]
    zero_sigmas_v = sampled_sigmas_v[sampled_bits == 0]
    one_sigmas_v = sampled_sigmas_v[sampled_bits == 1]
    if zero_means_v.size == 0 or one_means_v.size == 0:
        msg = "each phase must contain at least one sampled 0 and 1 bit."
        raise ValueError(msg)

    mean_zero_v = float(np.mean(zero_means_v))
    mean_one_v = float(np.mean(one_means_v))
    one_is_high = mean_one_v >= mean_zero_v
    threshold_v, ber = optimize_gaussian_mixture_threshold(
        zero_means_v,
        one_means_v,
        zero_sigmas_v,
        one_sigmas_v,
        one_is_high,
    )
    effective_sigma_zero_v = sqrt(
        float(np.mean(zero_sigmas_v**2 + (zero_means_v - mean_zero_v) ** 2))
    )
    effective_sigma_one_v = sqrt(
        float(np.mean(one_sigmas_v**2 + (one_means_v - mean_one_v) ** 2))
    )
    effective_q = abs(mean_one_v - mean_zero_v) / (
        effective_sigma_zero_v + effective_sigma_one_v
    )
    return StatisticalEyePhaseResult(
        sampling_phase_ui=sampling_phase_ui,
        decision_delay_s=decision_delay_s,
        one_is_high=one_is_high,
        threshold_v=threshold_v,
        ber=ber,
        effective_q=effective_q,
        mean_zero_v=mean_zero_v,
        mean_one_v=mean_one_v,
        rms_zero_noise_v=sqrt(float(np.mean(zero_sigmas_v**2))),
        rms_one_noise_v=sqrt(float(np.mean(one_sigmas_v**2))),
        zero_means_v=zero_means_v,
        one_means_v=one_means_v,
        zero_sigmas_v=zero_sigmas_v,
        one_sigmas_v=one_sigmas_v,
    )


def _sample_positions(
    num_symbols: int,
    num_samples: int,
    samples_per_symbol: int,
    sample_rate_hz: float,
    sampling_phase_ui: float,
    decision_delay_s: float,
    bits: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    bit_indices = np.arange(num_symbols)
    positions = (
        bit_indices * samples_per_symbol
        + sampling_phase_ui * samples_per_symbol
        + decision_delay_s * sample_rate_hz
    )
    valid = (positions >= 0) & (positions <= num_samples - 1)
    if not np.any(valid):
        msg = "decision timing does not overlap the statistical waveform."
        raise ValueError(msg)
    return positions[valid], np.asarray(bits[bit_indices[valid]], dtype=np.int_)


def _interpolated_variance(
    variance_v2: np.ndarray,
    lag1_covariance_v2: np.ndarray,
    positions: np.ndarray,
) -> np.ndarray:
    lower_indices = np.floor(positions).astype(int)
    fractions = positions - lower_indices
    upper_indices = np.minimum(lower_indices + 1, variance_v2.size - 1)
    interpolated = (
        (1 - fractions) ** 2 * variance_v2[lower_indices]
        + fractions**2 * variance_v2[upper_indices]
        + 2 * fractions * (1 - fractions) * lag1_covariance_v2[lower_indices]
    )
    scale = max(float(np.max(variance_v2)), np.finfo(float).tiny)
    if np.any(interpolated < -128 * np.finfo(float).eps * scale):
        msg = "interpolated noise variance became physically negative."
        raise ValueError(msg)
    return np.maximum(interpolated, 0.0)


def _validated_mixture_arrays(
    zero_means_v: np.ndarray,
    one_means_v: np.ndarray,
    zero_sigmas_v: np.ndarray,
    one_sigmas_v: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    zero_means = readonly_float_array(zero_means_v, "zero_means_v")
    one_means = readonly_float_array(one_means_v, "one_means_v")
    zero_sigmas = readonly_float_array(zero_sigmas_v, "zero_sigmas_v")
    one_sigmas = readonly_float_array(one_sigmas_v, "one_sigmas_v")
    if zero_means.size != zero_sigmas.size:
        msg = "zero_means_v and zero_sigmas_v must have the same length."
        raise ValueError(msg)
    if one_means.size != one_sigmas.size:
        msg = "one_means_v and one_sigmas_v must have the same length."
        raise ValueError(msg)
    if np.any(zero_sigmas <= 0) or np.any(one_sigmas <= 0):
        msg = "Gaussian mixture sigmas must be positive."
        raise ValueError(msg)
    return zero_means, one_means, zero_sigmas, one_sigmas


def _transformed_mixture_log_ber_array(
    thresholds: np.ndarray,
    transformed_zero_means: np.ndarray,
    transformed_one_means: np.ndarray,
    zero_sigmas: np.ndarray,
    one_sigmas: np.ndarray,
    zero_weights: np.ndarray,
    one_weights: np.ndarray,
) -> np.ndarray:
    """Calculate natural-log mixture BER without Gaussian-tail underflow."""
    threshold_values = np.asarray(thresholds, dtype=float)
    if threshold_values.ndim != 1 or threshold_values.size == 0:
        msg = "thresholds must be a non-empty one-dimensional array."
        raise ValueError(msg)
    if not np.all(np.isfinite(threshold_values)):
        msg = "thresholds must contain only finite values."
        raise ValueError(msg)

    maximum_components = max(zero_sigmas.size, one_sigmas.size)
    maximum_temporary_elements = 2_000_000
    chunk_size = max(
        1,
        min(
            threshold_values.size,
            maximum_temporary_elements // maximum_components,
        ),
    )
    log_zero_weights = np.log(zero_weights)
    log_one_weights = np.log(one_weights)
    output = np.empty(threshold_values.size, dtype=float)
    for start in range(0, threshold_values.size, chunk_size):
        stop = min(start + chunk_size, threshold_values.size)
        threshold_column = threshold_values[start:stop, np.newaxis]
        log_error_zero = logsumexp(
            log_ndtr(
                -(threshold_column - transformed_zero_means) / zero_sigmas,
            )
            + log_zero_weights,
            axis=1,
        )
        log_error_one = logsumexp(
            log_ndtr(
                -(transformed_one_means - threshold_column) / one_sigmas,
            )
            + log_one_weights,
            axis=1,
        )
        output[start:stop] = logsumexp(
            np.vstack((log_error_zero, log_error_one)), axis=0
        ) - np.log(2.0)
    return output


def _validated_mixture_weights(
    weights: np.ndarray | None,
    expected_size: int,
    name: str,
) -> np.ndarray:
    if weights is None:
        uniform = np.full(expected_size, 1.0 / expected_size)
        uniform.setflags(write=False)
        return uniform
    validated = readonly_float_array(weights, name)
    if validated.size != expected_size:
        msg = f"{name} length must match its mixture component count."
        raise ValueError(msg)
    if np.any(validated <= 0):
        msg = f"{name} values must be positive."
        raise ValueError(msg)
    normalized = np.array(validated / np.sum(validated), copy=True)
    normalized.setflags(write=False)
    return normalized


def _amplitude_axis(
    phase_results: tuple[StatisticalEyePhaseResult, ...],
    amplitude_points: int,
    density_sigma_span: float,
) -> np.ndarray:
    lower = min(
        float(np.min(result.zero_means_v - density_sigma_span * result.zero_sigmas_v))
        for result in phase_results
    )
    lower = min(
        lower,
        *(
            float(np.min(result.one_means_v - density_sigma_span * result.one_sigmas_v))
            for result in phase_results
        ),
    )
    upper = max(
        float(np.max(result.zero_means_v + density_sigma_span * result.zero_sigmas_v))
        for result in phase_results
    )
    upper = max(
        upper,
        *(
            float(np.max(result.one_means_v + density_sigma_span * result.one_sigmas_v))
            for result in phase_results
        ),
    )
    return np.linspace(lower, upper, amplitude_points)


def _phase_density(
    result: StatisticalEyePhaseResult,
    amplitude_axis_v: np.ndarray,
) -> np.ndarray:
    zero_density = _gaussian_mixture_density(
        amplitude_axis_v,
        result.zero_means_v,
        result.zero_sigmas_v,
    )
    one_density = _gaussian_mixture_density(
        amplitude_axis_v,
        result.one_means_v,
        result.one_sigmas_v,
    )
    return 0.5 * (zero_density + one_density)


def _gaussian_mixture_density(
    amplitude_axis_v: np.ndarray,
    means_v: np.ndarray,
    sigmas_v: np.ndarray,
) -> np.ndarray:
    normalized = (amplitude_axis_v[:, np.newaxis] - means_v[np.newaxis, :]) / sigmas_v[
        np.newaxis, :
    ]
    component_density = np.exp(-0.5 * normalized**2) / (
        sqrt(2 * np.pi) * sigmas_v[np.newaxis, :]
    )
    return np.mean(component_density, axis=1)
