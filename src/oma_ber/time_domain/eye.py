"""Clock-referenced deterministic eye analysis for current or voltage waveforms."""

from dataclasses import dataclass, field
from math import isfinite, log10
from typing import Literal

import numpy as np

from oma_ber.time_domain._validation import readonly_float_array, validated_bits
from oma_ber.time_domain.timebase import TimeGrid


AmplitudeUnit = Literal["A", "V"]


@dataclass(frozen=True)
class EyeLevelMetrics:
    """Deterministic sampled-eye metrics at one clock phase."""

    sampling_phase_ui: float
    decision_delay_s: float
    mean_zero_level: float
    mean_one_level: float
    min_one_level: float
    max_zero_level: float
    lower_inner_level: float
    upper_inner_level: float
    one_is_high: bool
    vertical_eye_opening: float
    mean_level_separation: float
    isi_penalty_db: float
    zero_sample_count: int
    one_sample_count: int


@dataclass(frozen=True)
class DeterministicEyeAnalysis:
    """Deterministic eye traces and phase-dependent level metrics."""

    amplitude_unit: AmplitudeUnit
    time_grid: TimeGrid
    decision_delay_s: float
    metrics_by_phase: tuple[EyeLevelMetrics, ...]
    optimum_metrics: EyeLevelMetrics
    eye_time_ui: np.ndarray = field(repr=False)
    eye_traces: np.ndarray = field(repr=False)

    def __post_init__(self) -> None:
        if self.amplitude_unit not in {"A", "V"}:
            msg = "amplitude_unit must be 'A' or 'V'."
            raise ValueError(msg)
        axis = readonly_float_array(self.eye_time_ui, "eye_time_ui")
        traces = np.asarray(self.eye_traces, dtype=float)
        if traces.ndim != 2 or traces.shape[0] == 0:
            msg = "eye_traces must be a non-empty two-dimensional array."
            raise ValueError(msg)
        if traces.shape[1] != axis.size:
            msg = "eye_traces width must match eye_time_ui length."
            raise ValueError(msg)
        if not np.all(np.isfinite(traces)):
            msg = "eye_traces must contain only finite values."
            raise ValueError(msg)
        copied_traces = np.array(traces, copy=True)
        copied_traces.setflags(write=False)
        object.__setattr__(self, "eye_time_ui", axis)
        object.__setattr__(self, "eye_traces", copied_traces)


def sample_eye_at_phase(
    waveform: np.ndarray,
    bits: np.ndarray,
    time_grid: TimeGrid,
    sampling_phase_ui: float,
    decision_delay_s: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample a waveform at a fractional clock phase for matching bits.

    sampling_phase_ui is in [0, 1) UI. decision_delay_s is an explicit signed
    channel/clock alignment delay in seconds. Linear interpolation is used at
    fractional sample positions. Only bit/sample pairs inside the waveform are
    returned.
    """
    values = readonly_float_array(waveform, "waveform")
    bit_values = validated_bits(bits)
    expected_samples = bit_values.size * time_grid.samples_per_symbol
    if values.size != expected_samples:
        msg = "waveform length must equal bits * samples_per_symbol."
        raise ValueError(msg)
    if not isfinite(sampling_phase_ui) or not 0 <= sampling_phase_ui < 1:
        msg = "sampling_phase_ui must be finite and in [0, 1)."
        raise ValueError(msg)
    if not isfinite(decision_delay_s):
        msg = "decision_delay_s must be finite."
        raise ValueError(msg)

    samples_per_symbol = time_grid.samples_per_symbol
    bit_indices = np.arange(bit_values.size)
    sample_positions = (
        bit_indices * samples_per_symbol
        + sampling_phase_ui * samples_per_symbol
        + decision_delay_s * time_grid.sample_rate_hz
    )
    valid = (sample_positions >= 0) & (sample_positions <= values.size - 1)
    if not np.any(valid):
        msg = "decision timing does not overlap the waveform."
        raise ValueError(msg)

    sampled_values = np.interp(
        sample_positions[valid],
        np.arange(values.size, dtype=float),
        values,
    )
    sampled_bits = np.asarray(bit_values[bit_indices[valid]], dtype=np.int_)
    return sampled_values, sampled_bits


def analyze_deterministic_eye(
    waveform: np.ndarray,
    bits: np.ndarray,
    time_grid: TimeGrid,
    amplitude_unit: AmplitudeUnit,
    decision_delay_s: float = 0.0,
    sampling_phases_ui: np.ndarray | None = None,
) -> DeterministicEyeAnalysis:
    """Analyze a deterministic 2 UI eye and select maximum-opening phase.

    The selected phase maximizes the separation between the inner 0/1 levels,
    with either normal or inverted receiver polarity. This is an ISI-only
    diagnostic, not a CDR or jitter-aware optimum. Ties are resolved toward
    0.5 UI. Random-noise extrema are intentionally not accepted here.
    """
    values = readonly_float_array(waveform, "waveform")
    bit_values = validated_bits(bits)
    if amplitude_unit not in {"A", "V"}:
        msg = "amplitude_unit must be 'A' or 'V'."
        raise ValueError(msg)
    if not isfinite(decision_delay_s):
        msg = "decision_delay_s must be finite."
        raise ValueError(msg)

    if sampling_phases_ui is None:
        phases = np.arange(time_grid.samples_per_symbol, dtype=float) / (
            time_grid.samples_per_symbol
        )
    else:
        phases = readonly_float_array(sampling_phases_ui, "sampling_phases_ui")
        if np.any((phases < 0) | (phases >= 1)):
            msg = "sampling_phases_ui values must be in [0, 1)."
            raise ValueError(msg)
        if np.unique(phases).size != phases.size:
            msg = "sampling_phases_ui values must be unique."
            raise ValueError(msg)

    metrics = tuple(
        _metrics_at_phase(
            waveform=values,
            bits=bit_values,
            time_grid=time_grid,
            sampling_phase_ui=float(phase),
            decision_delay_s=decision_delay_s,
        )
        for phase in phases
    )
    optimum = _choose_optimum_phase(metrics)
    eye_time_ui, eye_traces = _centered_eye_traces(
        waveform=values,
        num_symbols=bit_values.size,
        time_grid=time_grid,
        sampling_phase_ui=optimum.sampling_phase_ui,
        decision_delay_s=decision_delay_s,
    )
    return DeterministicEyeAnalysis(
        amplitude_unit=amplitude_unit,
        time_grid=time_grid,
        decision_delay_s=decision_delay_s,
        metrics_by_phase=metrics,
        optimum_metrics=optimum,
        eye_time_ui=eye_time_ui,
        eye_traces=eye_traces,
    )


def _metrics_at_phase(
    waveform: np.ndarray,
    bits: np.ndarray,
    time_grid: TimeGrid,
    sampling_phase_ui: float,
    decision_delay_s: float,
) -> EyeLevelMetrics:
    samples, sampled_bits = sample_eye_at_phase(
        waveform=waveform,
        bits=bits,
        time_grid=time_grid,
        sampling_phase_ui=sampling_phase_ui,
        decision_delay_s=decision_delay_s,
    )
    zero_samples = samples[sampled_bits == 0]
    one_samples = samples[sampled_bits == 1]
    if zero_samples.size == 0 or one_samples.size == 0:
        msg = "each sampling phase must contain at least one 0 and one 1 sample."
        raise ValueError(msg)

    mean_zero = float(np.mean(zero_samples))
    mean_one = float(np.mean(one_samples))
    min_one = float(np.min(one_samples))
    max_zero = float(np.max(zero_samples))
    one_is_high = mean_one >= mean_zero
    if one_is_high:
        lower_inner_level = max_zero
        upper_inner_level = min_one
    else:
        lower_inner_level = float(np.max(one_samples))
        upper_inner_level = float(np.min(zero_samples))
    opening = upper_inner_level - lower_inner_level
    separation = abs(mean_one - mean_zero)
    isi_penalty_db = float("inf") if opening <= 0 else 20 * log10(separation / opening)
    return EyeLevelMetrics(
        sampling_phase_ui=sampling_phase_ui,
        decision_delay_s=decision_delay_s,
        mean_zero_level=mean_zero,
        mean_one_level=mean_one,
        min_one_level=min_one,
        max_zero_level=max_zero,
        lower_inner_level=lower_inner_level,
        upper_inner_level=upper_inner_level,
        one_is_high=one_is_high,
        vertical_eye_opening=opening,
        mean_level_separation=separation,
        isi_penalty_db=isi_penalty_db,
        zero_sample_count=int(zero_samples.size),
        one_sample_count=int(one_samples.size),
    )


def _choose_optimum_phase(metrics: tuple[EyeLevelMetrics, ...]) -> EyeLevelMetrics:
    maximum_opening = max(metric.vertical_eye_opening for metric in metrics)
    scale = max(
        abs(maximum_opening),
        *(abs(metric.mean_zero_level) for metric in metrics),
        *(abs(metric.mean_one_level) for metric in metrics),
        np.finfo(float).tiny,
    )
    tolerance = 64 * np.finfo(float).eps * scale
    candidates = [
        metric
        for metric in metrics
        if abs(metric.vertical_eye_opening - maximum_opening) <= tolerance
    ]
    return min(candidates, key=lambda metric: abs(metric.sampling_phase_ui - 0.5))


def _centered_eye_traces(
    waveform: np.ndarray,
    num_symbols: int,
    time_grid: TimeGrid,
    sampling_phase_ui: float,
    decision_delay_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    samples_per_symbol = time_grid.samples_per_symbol
    offsets = np.arange(-samples_per_symbol, samples_per_symbol + 1, dtype=float)
    centers = (
        np.arange(num_symbols, dtype=float) * samples_per_symbol
        + sampling_phase_ui * samples_per_symbol
        + decision_delay_s * time_grid.sample_rate_hz
    )
    valid_centers = centers[
        (centers + offsets[0] >= 0) & (centers + offsets[-1] <= waveform.size - 1)
    ]
    if valid_centers.size == 0:
        msg = "waveform does not contain a complete 2 UI trace at the decision timing."
        raise ValueError(msg)

    sample_axis = np.arange(waveform.size, dtype=float)
    traces = np.vstack(
        [np.interp(center + offsets, sample_axis, waveform) for center in valid_centers]
    )
    return offsets / samples_per_symbol, traces
