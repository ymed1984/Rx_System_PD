"""Matplotlib helpers for visualizing OMA-to-BER results."""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterSciNotation

from oma_ber.eye import eye_traces, eye_unit_interval_axis
from oma_ber.time_domain.eye import DeterministicEyeAnalysis
from oma_ber.time_domain.statistical import StatisticalEyeAnalysis


def _require_keys(sweep_results: list[dict], required_keys: set[str]) -> None:
    for index, result in enumerate(sweep_results):
        missing_keys = required_keys - result.keys()
        if missing_keys:
            keys = ", ".join(sorted(missing_keys))
            msg = f"sweep_results[{index}] is missing required key(s): {keys}."
            raise ValueError(msg)


def _get_axes(ax: Axes | None) -> Axes:
    if ax is not None:
        return ax
    _, new_ax = plt.subplots()
    return new_ax


def plot_oma_sweep(
    sweep_results: list[dict],
    ax: Axes | None = None,
) -> Axes:
    """Plot BER versus OMA in dBm for sweep result dictionaries.

    The BER axis uses a display-only lower limit of 1e-12 so very small BER
    values do not make the practical receiver region unreadable. Data values
    are not clipped or modified.
    """
    _require_keys(sweep_results, {"oma_dbm", "ber"})
    oma_dbm_values = np.array(
        [result["oma_dbm"] for result in sweep_results], dtype=float
    )
    ber_values = np.array([result["ber"] for result in sweep_results], dtype=float)

    if np.any(ber_values <= 0):
        msg = "BER values must be positive for logarithmic plotting."
        raise ValueError(msg)

    ax = _get_axes(ax)
    ax.plot(oma_dbm_values, ber_values, marker="o")
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(LogFormatterSciNotation())
    ax.set_ylim(bottom=1e-12)
    ax.set_xlabel("OMA [dBm]")
    ax.set_ylabel("BER (axis floor 1e-12)")
    ax.grid(True, which="both", alpha=0.35)
    return ax


def plot_q_vs_oma(
    sweep_results: list[dict],
    ax: Axes | None = None,
) -> Axes:
    """Plot receiver Q estimate versus OMA in dBm for sweep result dictionaries."""
    _require_keys(sweep_results, {"oma_dbm", "q_rx"})
    oma_dbm_values = np.array(
        [result["oma_dbm"] for result in sweep_results], dtype=float
    )
    q_values = np.array([result["q_rx"] for result in sweep_results], dtype=float)

    ax = _get_axes(ax)
    ax.plot(oma_dbm_values, q_values, marker="o")
    ax.set_xlabel("OMA [dBm]")
    ax.set_ylabel("Q")
    ax.grid(True, alpha=0.35)
    return ax


def plot_required_oma_vs_responsivity(
    responsivity_a_per_w_values: np.ndarray,
    required_oma_dbm_values: np.ndarray,
    ax: Axes | None = None,
) -> Axes:
    """Plot required OMA in dBm versus photodiode responsivity in A/W."""
    if len(responsivity_a_per_w_values) != len(required_oma_dbm_values):
        msg = "responsivity_a_per_w_values and required_oma_dbm_values must have the same length."
        raise ValueError(msg)

    ax = _get_axes(ax)
    ax.plot(responsivity_a_per_w_values, required_oma_dbm_values, marker="o")
    ax.set_xlabel("Responsivity [A/W]")
    ax.set_ylabel("Required OMA [dBm]")
    ax.grid(True, alpha=0.35)
    return ax


def plot_eye_diagram(
    waveform: np.ndarray,
    samples_per_symbol: int,
    symbols_per_trace: int = 2,
    sample_rate_hz: float | None = None,
    max_traces: int | None = 200,
    eye_metrics: dict[str, float] | None = None,
    ax: Axes | None = None,
    ylabel: str = "Amplitude",
) -> Axes:
    """Plot an eye diagram from a deterministic waveform.

    waveform uses arbitrary but consistent amplitude units, such as optical
    power in W or current in A. samples_per_symbol is in samples/symbol.
    If sample_rate_hz is provided, the x-axis is shown in ps; otherwise it is
    shown in unit intervals (UI). eye_metrics may be a sampled_eye_levels()
    result dictionary; when provided, min/max eye-opening guides are drawn.
    """
    if sample_rate_hz is not None and sample_rate_hz <= 0:
        msg = "sample_rate_hz must be positive when provided."
        raise ValueError(msg)
    if max_traces is not None and max_traces <= 0:
        msg = "max_traces must be positive when provided."
        raise ValueError(msg)

    traces = eye_traces(
        waveform=waveform,
        samples_per_symbol=samples_per_symbol,
        symbols_per_trace=symbols_per_trace,
    )
    if max_traces is not None:
        traces = traces[:max_traces]

    if sample_rate_hz is None:
        x_values = eye_unit_interval_axis(samples_per_symbol, symbols_per_trace)
        xlabel = "Time [UI]"
    else:
        x_values = np.arange(traces.shape[1], dtype=float) / sample_rate_hz * 1e12
        xlabel = "Time [ps]"

    ax = _get_axes(ax)
    for trace in traces:
        ax.plot(x_values, trace, color="C0", alpha=0.22, linewidth=0.8)

    if eye_metrics is not None:
        _require_eye_metric_keys(eye_metrics)
        ax.axhline(
            eye_metrics["min_one_level"], color="C2", linestyle="--", linewidth=1.0
        )
        ax.axhline(
            eye_metrics["max_zero_level"], color="C3", linestyle="--", linewidth=1.0
        )
        ax.axhspan(
            eye_metrics["max_zero_level"],
            eye_metrics["min_one_level"],
            color="C2",
            alpha=0.08,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(float(x_values[0]), float(x_values[-1]))
    ax.grid(True, alpha=0.3)
    return ax


def plot_deterministic_eye_analysis(
    analysis: DeterministicEyeAnalysis,
    max_traces: int | None = 200,
    ax: Axes | None = None,
) -> Axes:
    """Plot a clock-centered deterministic current or voltage eye analysis.

    analysis uses the optimum deterministic sampling phase selected by
    analyze_deterministic_eye(). The x-axis is centered on the decision instant
    at 0 UI. Current eyes are labeled in A and voltage eyes in V. If max_traces
    limits the display, traces are selected uniformly across the full pattern.
    """
    if not isinstance(analysis, DeterministicEyeAnalysis):
        msg = "analysis must be a DeterministicEyeAnalysis."
        raise ValueError(msg)
    if max_traces is not None and max_traces <= 0:
        msg = "max_traces must be positive when provided."
        raise ValueError(msg)

    traces = analysis.eye_traces
    if max_traces is not None and traces.shape[0] > max_traces:
        trace_indices = np.linspace(0, traces.shape[0] - 1, max_traces, dtype=int)
        traces = traces[trace_indices]

    ax = _get_axes(ax)
    for trace in traces:
        ax.plot(
            analysis.eye_time_ui,
            trace,
            color="C0",
            alpha=0.22,
            linewidth=0.8,
        )

    optimum = analysis.optimum_metrics
    ax.axvline(0.0, color="black", linestyle=":", linewidth=1.0)
    ax.axhline(optimum.upper_inner_level, color="C2", linestyle="--", linewidth=1.0)
    ax.axhline(optimum.lower_inner_level, color="C3", linestyle="--", linewidth=1.0)
    if optimum.vertical_eye_opening > 0:
        ax.axhspan(
            optimum.lower_inner_level,
            optimum.upper_inner_level,
            color="C2",
            alpha=0.08,
        )

    ax.set_xlabel("Time relative to decision [UI]")
    ax.set_ylabel(f"Receiver amplitude [{analysis.amplitude_unit}]")
    ax.set_xlim(float(analysis.eye_time_ui[0]), float(analysis.eye_time_ui[-1]))
    ax.grid(True, alpha=0.3)
    return ax


def plot_statistical_eye_density(
    analysis: StatisticalEyeAnalysis,
    density_floor_ratio: float = 1e-8,
    ax: Axes | None = None,
) -> Axes:
    """Plot a two-UI Gaussian-mixture TIA-voltage eye-density map.

    Time is centered on the BER-optimum decision phase. Color uses logarithmic
    probability density so low-probability tails remain visible. The horizontal
    line is the optimized physical TIA voltage threshold in V.
    """
    if not isinstance(analysis, StatisticalEyeAnalysis):
        msg = "analysis must be a StatisticalEyeAnalysis."
        raise ValueError(msg)
    if not 0 < density_floor_ratio < 1:
        msg = "density_floor_ratio must be between 0 and 1."
        raise ValueError(msg)

    order = np.argsort(analysis.sampling_phases_ui)
    phases = analysis.sampling_phases_ui[order]
    density = analysis.density_per_v[order]
    extended_phases = np.concatenate((phases - 1, phases, phases + 1))
    extended_density = np.vstack((density, density, density))
    optimum_phase = analysis.optimum_result.sampling_phase_ui
    selected = (extended_phases >= optimum_phase - 1) & (
        extended_phases < optimum_phase + 1
    )
    relative_phase_ui = extended_phases[selected] - optimum_phase
    selected_density = extended_density[selected]

    maximum_density = float(np.max(selected_density))
    display_density = np.maximum(
        selected_density,
        maximum_density * density_floor_ratio,
    )
    ax = _get_axes(ax)
    ax.pcolormesh(
        relative_phase_ui,
        analysis.amplitude_axis_v,
        display_density.T,
        shading="auto",
        cmap="viridis",
        norm=LogNorm(
            vmin=maximum_density * density_floor_ratio,
            vmax=maximum_density,
        ),
    )
    ax.axvline(0.0, color="white", linestyle=":", linewidth=1.0)
    ax.axhline(
        analysis.optimum_result.threshold_v,
        color="white",
        linestyle="--",
        linewidth=1.0,
    )
    ax.set_xlabel("Time relative to BER-optimum decision [UI]")
    ax.set_ylabel("TIA output voltage [V]")
    ax.set_xlim(-1.0, 1.0)
    return ax


def plot_statistical_ber_vs_phase(
    analysis: StatisticalEyeAnalysis,
    ax: Axes | None = None,
) -> Axes:
    """Plot optimized Gaussian-mixture BER versus sampling phase in UI."""
    if not isinstance(analysis, StatisticalEyeAnalysis):
        msg = "analysis must be a StatisticalEyeAnalysis."
        raise ValueError(msg)

    phases = np.array(
        [result.sampling_phase_ui for result in analysis.phase_results],
        dtype=float,
    )
    ber_values = np.array(
        [max(result.ber, np.finfo(float).tiny) for result in analysis.phase_results],
        dtype=float,
    )
    optimum = analysis.optimum_result
    ax = _get_axes(ax)
    ax.semilogy(phases, ber_values, marker="o")
    ax.semilogy(
        [optimum.sampling_phase_ui],
        [max(optimum.ber, np.finfo(float).tiny)],
        marker="*",
        markersize=11,
        color="C3",
    )
    ax.set_xlabel("Sampling phase [UI]")
    ax.set_ylabel("Optimized Gaussian-mixture BER")
    ax.set_xlim(0.0, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    return ax


def _require_eye_metric_keys(eye_metrics: dict[str, float]) -> None:
    required_keys = {"min_one_level", "max_zero_level"}
    missing_keys = required_keys - eye_metrics.keys()
    if missing_keys:
        keys = ", ".join(sorted(missing_keys))
        msg = f"eye_metrics is missing required key(s): {keys}."
        raise ValueError(msg)
