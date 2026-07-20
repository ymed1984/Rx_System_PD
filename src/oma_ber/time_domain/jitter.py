"""Residual timing-jitter averaging and BER bathtub analysis."""

from dataclasses import dataclass, field
from math import isfinite, pi, sqrt

import numpy as np
from numpy.polynomial.hermite import hermgauss

from oma_ber.time_domain._validation import readonly_float_array
from oma_ber.time_domain.noise import TiaOutputNoiseResult
from oma_ber.time_domain.receiver import ReceiverWaveformResult
from oma_ber.time_domain.statistical import optimize_gaussian_mixture_threshold


@dataclass(frozen=True)
class ResidualTimingJitter:
    """Timing error remaining between the data eye and decision clock.

    random_jitter_rms_s is the standard deviation of zero-mean Gaussian random
    jitter in seconds. sinusoidal_jitter_peak_s is the peak time deviation of a
    sinusoid whose phase is uniformly distributed. These are residual sampler
    errors, not pre-CDR input jitter specifications.
    """

    random_jitter_rms_s: float = 0.0
    sinusoidal_jitter_peak_s: float = 0.0

    def __post_init__(self) -> None:
        for value, name in (
            (self.random_jitter_rms_s, "random_jitter_rms_s"),
            (self.sinusoidal_jitter_peak_s, "sinusoidal_jitter_peak_s"),
        ):
            if not isfinite(value) or value < 0:
                msg = f"{name} must be finite and non-negative."
                raise ValueError(msg)

    @property
    def total_rms_s(self) -> float:
        """Return RMS time error in seconds for independent RJ and sinusoidal jitter."""
        return sqrt(self.random_jitter_rms_s**2 + self.sinusoidal_jitter_peak_s**2 / 2)


@dataclass(frozen=True)
class JitteredEyePhaseResult:
    """Fixed-threshold jitter-averaged result at one nominal phase."""

    nominal_sampling_phase_ui: float
    decision_delay_s: float
    one_is_high: bool
    threshold_v: float
    ber: float
    effective_q: float
    mean_zero_v: float
    mean_one_v: float
    rms_zero_noise_v: float
    rms_one_noise_v: float
    quadrature_nodes: int

    def __post_init__(self) -> None:
        if (
            not isfinite(self.nominal_sampling_phase_ui)
            or not 0 <= self.nominal_sampling_phase_ui < 1
        ):
            msg = "nominal_sampling_phase_ui must be finite and in [0, 1)."
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
        scalar_values = (
            self.effective_q,
            self.mean_zero_v,
            self.mean_one_v,
            self.rms_zero_noise_v,
            self.rms_one_noise_v,
        )
        if not all(isfinite(value) for value in scalar_values):
            msg = "jittered eye voltage summaries must be finite."
            raise ValueError(msg)
        if self.effective_q < 0:
            msg = "effective_q must be non-negative."
            raise ValueError(msg)
        if self.rms_zero_noise_v <= 0 or self.rms_one_noise_v <= 0:
            msg = "jittered eye RMS noise values must be positive."
            raise ValueError(msg)
        if isinstance(self.quadrature_nodes, bool) or not isinstance(
            self.quadrature_nodes, int
        ):
            msg = "quadrature_nodes must be an integer."
            raise ValueError(msg)
        if self.quadrature_nodes < 1:
            msg = "quadrature_nodes must be positive."
            raise ValueError(msg)


@dataclass(frozen=True)
class JitteredEyeAnalysis:
    """BER bathtub and optimum fixed-threshold residual-jitter result."""

    phase_results: tuple[JitteredEyePhaseResult, ...]
    optimum_result: JitteredEyePhaseResult
    jitter: ResidualTimingJitter
    nominal_phases_ui: np.ndarray = field(repr=False)
    optimized_thresholds_v: np.ndarray = field(repr=False)
    jitter_averaged_ber: np.ndarray = field(repr=False)
    random_jitter_rms_ui: float
    sinusoidal_jitter_peak_ui: float

    def __post_init__(self) -> None:
        phases = readonly_float_array(self.nominal_phases_ui, "nominal_phases_ui")
        thresholds = readonly_float_array(
            self.optimized_thresholds_v,
            "optimized_thresholds_v",
        )
        ber = readonly_float_array(self.jitter_averaged_ber, "jitter_averaged_ber")
        if np.any((phases < 0) | (phases >= 1)):
            msg = "nominal_phases_ui values must be in [0, 1)."
            raise ValueError(msg)
        if np.unique(phases).size != phases.size:
            msg = "nominal_phases_ui values must be unique."
            raise ValueError(msg)
        if thresholds.size != phases.size or ber.size != phases.size:
            msg = "jittered bathtub arrays must have matching lengths."
            raise ValueError(msg)
        if np.any((ber < 0) | (ber > 0.5)):
            msg = "jitter_averaged_ber values must be in [0, 0.5]."
            raise ValueError(msg)
        if len(self.phase_results) != phases.size:
            msg = "phase_results length must match nominal_phases_ui."
            raise ValueError(msg)
        if not any(result is self.optimum_result for result in self.phase_results):
            msg = "optimum_result must be one of phase_results."
            raise ValueError(msg)
        if not isinstance(self.jitter, ResidualTimingJitter):
            msg = "jitter must be a ResidualTimingJitter."
            raise ValueError(msg)
        if (
            not isfinite(self.random_jitter_rms_ui)
            or self.random_jitter_rms_ui < 0
            or not isfinite(self.sinusoidal_jitter_peak_ui)
            or self.sinusoidal_jitter_peak_ui < 0
        ):
            msg = "jitter values in UI must be finite and non-negative."
            raise ValueError(msg)
        object.__setattr__(self, "nominal_phases_ui", phases)
        object.__setattr__(self, "optimized_thresholds_v", thresholds)
        object.__setattr__(self, "jitter_averaged_ber", ber)


def analyze_jittered_tia_eye(
    waveforms: ReceiverWaveformResult,
    noise: TiaOutputNoiseResult,
    jitter: ResidualTimingJitter,
    nominal_phases_ui: np.ndarray | None = None,
    decision_delay_s: float = 0.0,
    random_quadrature_points: int = 17,
    sinusoidal_phase_points: int = 32,
    threshold_grid_points: int = 257,
) -> JitteredEyeAnalysis:
    """Calculate a fixed-threshold residual-jitter BER bathtub.

    A single voltage threshold is optimized for each nominal phase after
    averaging over all timing-error quadrature nodes. The threshold is not
    re-optimized for each instantaneous jitter value. Finite signal and noise
    waveforms are treated as periodic PRBS continuations while the decision bit
    label remains fixed when a timing offset crosses a UI boundary.
    """
    if waveforms.tia_output_voltage_v is None:
        msg = "jittered TIA eye requires a TIA voltage waveform."
        raise ValueError(msg)
    if noise.total_variance_v2.size != waveforms.tia_output_voltage_v.size:
        msg = "noise length must match the TIA voltage waveform length."
        raise ValueError(msg)
    if not isinstance(jitter, ResidualTimingJitter):
        msg = "jitter must be a ResidualTimingJitter."
        raise ValueError(msg)
    if not isfinite(decision_delay_s):
        msg = "decision_delay_s must be finite."
        raise ValueError(msg)
    _validate_quadrature_points(
        random_quadrature_points,
        sinusoidal_phase_points,
        threshold_grid_points,
    )
    if nominal_phases_ui is None:
        phases = np.arange(waveforms.time_grid.samples_per_symbol, dtype=float) / (
            waveforms.time_grid.samples_per_symbol
        )
    else:
        phases = readonly_float_array(nominal_phases_ui, "nominal_phases_ui")
        if np.any((phases < 0) | (phases >= 1)):
            msg = "nominal_phases_ui values must be in [0, 1)."
            raise ValueError(msg)
        if np.unique(phases).size != phases.size:
            msg = "nominal_phases_ui values must be unique."
            raise ValueError(msg)

    timing_offsets_s, timing_weights = _timing_quadrature(
        jitter,
        random_quadrature_points,
        sinusoidal_phase_points,
    )
    phase_results = tuple(
        _analyze_jittered_phase(
            waveforms=waveforms,
            noise=noise,
            nominal_phase_ui=float(phase),
            decision_delay_s=decision_delay_s,
            timing_offsets_s=timing_offsets_s,
            timing_weights=timing_weights,
            threshold_grid_points=threshold_grid_points,
        )
        for phase in phases
    )
    optimum_result = min(
        phase_results,
        key=lambda result: (
            result.ber,
            -result.effective_q,
            abs(result.nominal_sampling_phase_ui - 0.5),
        ),
    )
    unit_interval_s = waveforms.time_grid.unit_interval_s
    return JitteredEyeAnalysis(
        phase_results=phase_results,
        optimum_result=optimum_result,
        jitter=jitter,
        nominal_phases_ui=phases,
        optimized_thresholds_v=np.array(
            [result.threshold_v for result in phase_results]
        ),
        jitter_averaged_ber=np.array([result.ber for result in phase_results]),
        random_jitter_rms_ui=jitter.random_jitter_rms_s / unit_interval_s,
        sinusoidal_jitter_peak_ui=jitter.sinusoidal_jitter_peak_s / unit_interval_s,
    )


def _analyze_jittered_phase(
    waveforms: ReceiverWaveformResult,
    noise: TiaOutputNoiseResult,
    nominal_phase_ui: float,
    decision_delay_s: float,
    timing_offsets_s: np.ndarray,
    timing_weights: np.ndarray,
    threshold_grid_points: int,
) -> JitteredEyePhaseResult:
    zero_means_parts: list[np.ndarray] = []
    one_means_parts: list[np.ndarray] = []
    zero_sigmas_parts: list[np.ndarray] = []
    one_sigmas_parts: list[np.ndarray] = []
    zero_weight_parts: list[np.ndarray] = []
    one_weight_parts: list[np.ndarray] = []

    for timing_offset_s, timing_weight in zip(
        timing_offsets_s,
        timing_weights,
        strict=True,
    ):
        sampling_phase_ui = nominal_phase_ui + (
            timing_offset_s / waveforms.time_grid.unit_interval_s
        )
        zero_means, one_means, zero_sigmas, one_sigmas = _periodic_pattern_components(
            waveforms,
            noise,
            sampling_phase_ui,
            decision_delay_s,
        )
        zero_means_parts.append(zero_means)
        one_means_parts.append(one_means)
        zero_sigmas_parts.append(zero_sigmas)
        one_sigmas_parts.append(one_sigmas)
        zero_weight_parts.append(
            np.full(zero_means.size, timing_weight / zero_means.size)
        )
        one_weight_parts.append(np.full(one_means.size, timing_weight / one_means.size))

    zero_means_v = np.concatenate(zero_means_parts)
    one_means_v = np.concatenate(one_means_parts)
    zero_sigmas_v = np.concatenate(zero_sigmas_parts)
    one_sigmas_v = np.concatenate(one_sigmas_parts)
    zero_weights = np.concatenate(zero_weight_parts)
    one_weights = np.concatenate(one_weight_parts)
    mean_zero_v = float(np.dot(zero_weights, zero_means_v))
    mean_one_v = float(np.dot(one_weights, one_means_v))
    one_is_high = mean_one_v >= mean_zero_v
    threshold_v, ber = optimize_gaussian_mixture_threshold(
        zero_means_v,
        one_means_v,
        zero_sigmas_v,
        one_sigmas_v,
        one_is_high,
        grid_points=threshold_grid_points,
        zero_weights=zero_weights,
        one_weights=one_weights,
    )
    effective_sigma_zero_v = sqrt(
        float(
            np.dot(
                zero_weights,
                zero_sigmas_v**2 + (zero_means_v - mean_zero_v) ** 2,
            )
        )
    )
    effective_sigma_one_v = sqrt(
        float(
            np.dot(
                one_weights,
                one_sigmas_v**2 + (one_means_v - mean_one_v) ** 2,
            )
        )
    )
    effective_q = abs(mean_one_v - mean_zero_v) / (
        effective_sigma_zero_v + effective_sigma_one_v
    )
    return JitteredEyePhaseResult(
        nominal_sampling_phase_ui=nominal_phase_ui,
        decision_delay_s=decision_delay_s,
        one_is_high=one_is_high,
        threshold_v=threshold_v,
        ber=ber,
        effective_q=effective_q,
        mean_zero_v=mean_zero_v,
        mean_one_v=mean_one_v,
        rms_zero_noise_v=sqrt(float(np.dot(zero_weights, zero_sigmas_v**2))),
        rms_one_noise_v=sqrt(float(np.dot(one_weights, one_sigmas_v**2))),
        quadrature_nodes=timing_offsets_s.size,
    )


def _periodic_pattern_components(
    waveforms: ReceiverWaveformResult,
    noise: TiaOutputNoiseResult,
    sampling_phase_ui: float,
    decision_delay_s: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample periodic waveforms while preserving each intended bit label."""
    if not isfinite(sampling_phase_ui):
        msg = "sampling_phase_ui must be finite."
        raise ValueError(msg)
    voltage_v = waveforms.tia_output_voltage_v
    if voltage_v is None:
        msg = "periodic mixture sampling requires a TIA voltage waveform."
        raise ValueError(msg)
    num_samples = voltage_v.size
    samples_per_symbol = waveforms.time_grid.samples_per_symbol
    unwrapped_positions = (
        np.arange(waveforms.bits.size, dtype=float) * samples_per_symbol
        + sampling_phase_ui * samples_per_symbol
        + decision_delay_s * waveforms.time_grid.sample_rate_hz
    )
    lower_unwrapped = np.floor(unwrapped_positions)
    fractions = unwrapped_positions - lower_unwrapped
    lower_indices = lower_unwrapped.astype(np.int64) % num_samples
    upper_indices = (lower_indices + 1) % num_samples
    sampled_means_v = (1 - fractions) * voltage_v[
        lower_indices
    ] + fractions * voltage_v[upper_indices]
    sampled_variance_v2 = (
        (1 - fractions) ** 2 * noise.total_variance_v2[lower_indices]
        + fractions**2 * noise.total_variance_v2[upper_indices]
        + 2
        * fractions
        * (1 - fractions)
        * noise.total_lag1_covariance_v2[lower_indices]
    )
    variance_scale = max(float(np.max(noise.total_variance_v2)), np.finfo(float).tiny)
    if np.any(sampled_variance_v2 < -128 * np.finfo(float).eps * variance_scale):
        msg = "periodic interpolated noise variance became physically negative."
        raise ValueError(msg)
    sampled_sigmas_v = np.sqrt(np.maximum(sampled_variance_v2, 0.0))
    if np.any(sampled_sigmas_v <= 0):
        msg = "jittered Gaussian eye requires positive sampled noise."
        raise ValueError(msg)
    zero_mask = waveforms.bits == 0
    one_mask = ~zero_mask
    if not np.any(zero_mask) or not np.any(one_mask):
        msg = "jittered eye requires at least one 0 bit and one 1 bit."
        raise ValueError(msg)
    return (
        sampled_means_v[zero_mask],
        sampled_means_v[one_mask],
        sampled_sigmas_v[zero_mask],
        sampled_sigmas_v[one_mask],
    )


def _timing_quadrature(
    jitter: ResidualTimingJitter,
    random_quadrature_points: int,
    sinusoidal_phase_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    if jitter.random_jitter_rms_s == 0:
        random_offsets_s = np.array([0.0])
        random_weights = np.array([1.0])
    else:
        nodes, weights = hermgauss(random_quadrature_points)
        random_offsets_s = sqrt(2) * jitter.random_jitter_rms_s * nodes
        random_weights = weights / sqrt(pi)

    if jitter.sinusoidal_jitter_peak_s == 0:
        sinusoidal_offsets_s = np.array([0.0])
        sinusoidal_weights = np.array([1.0])
    else:
        phases_rad = (
            2
            * pi
            * (np.arange(sinusoidal_phase_points, dtype=float) + 0.5)
            / sinusoidal_phase_points
        )
        sinusoidal_offsets_s = jitter.sinusoidal_jitter_peak_s * np.sin(phases_rad)
        sinusoidal_weights = np.full(
            sinusoidal_phase_points,
            1.0 / sinusoidal_phase_points,
        )

    offsets_s = (
        random_offsets_s[:, np.newaxis] + sinusoidal_offsets_s[np.newaxis, :]
    ).ravel()
    weights = (
        random_weights[:, np.newaxis] * sinusoidal_weights[np.newaxis, :]
    ).ravel()
    normalized_weights = np.asarray(weights / np.sum(weights), dtype=float)
    offsets_s.setflags(write=False)
    normalized_weights.setflags(write=False)
    return offsets_s, normalized_weights


def _validate_quadrature_points(
    random_quadrature_points: int,
    sinusoidal_phase_points: int,
    threshold_grid_points: int,
) -> None:
    for value, name, minimum, maximum in (
        (random_quadrature_points, "random_quadrature_points", 3, 64),
        (sinusoidal_phase_points, "sinusoidal_phase_points", 4, 256),
        (threshold_grid_points, "threshold_grid_points", 33, 4097),
    ):
        if isinstance(value, bool) or not isinstance(value, int):
            msg = f"{name} must be an integer."
            raise ValueError(msg)
        if not minimum <= value <= maximum:
            msg = f"{name} must be in [{minimum}, {maximum}]."
            raise ValueError(msg)
