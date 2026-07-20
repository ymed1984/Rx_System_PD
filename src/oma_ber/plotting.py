"""Matplotlib helpers for visualizing OMA-to-BER results."""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterSciNotation

from oma_ber.eye import eye_traces, eye_unit_interval_axis
from oma_ber.time_domain.eye import DeterministicEyeAnalysis
from oma_ber.time_domain.jitter import JitteredEyeAnalysis
from oma_ber.time_domain.statistical import StatisticalEyeAnalysis
from oma_ber.time_domain.transfer import DiscreteTransferFunction


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


def plot_transfer_magnitude_comparison(
    transfers: dict[str, DiscreteTransferFunction],
    num_points: int = 1025,
    maximum_frequency_hz: float | None = None,
    normalize_to_dc: bool = True,
    ax: Axes | None = None,
) -> Axes:
    """Compare discrete PD or TIA amplitude responses versus frequency.

    All transfers must have the same sample rate and response kind. When
    normalize_to_dc is true, each curve is divided by its own absolute DC gain.
    Otherwise, transimpedance magnitude is displayed as 20*log10(|Z|/1 ohm)
    and a dimensionless response as 20*log10(|H|).
    """
    if not transfers:
        msg = "transfers must contain at least one response."
        raise ValueError(msg)
    if isinstance(num_points, bool) or not isinstance(num_points, int):
        msg = "num_points must be an integer."
        raise ValueError(msg)
    if num_points < 2:
        msg = "num_points must be at least 2."
        raise ValueError(msg)
    first_transfer = next(iter(transfers.values()))
    if not isinstance(first_transfer, DiscreteTransferFunction):
        msg = "all transfer values must be DiscreteTransferFunction instances."
        raise ValueError(msg)
    sample_rate_hz = first_transfer.sample_rate_hz
    response_kind = first_transfer.response_kind
    for name, transfer in transfers.items():
        if not isinstance(name, str) or not name:
            msg = "transfer labels must be non-empty strings."
            raise ValueError(msg)
        if not isinstance(transfer, DiscreteTransferFunction):
            msg = "all transfer values must be DiscreteTransferFunction instances."
            raise ValueError(msg)
        if not np.isclose(transfer.sample_rate_hz, sample_rate_hz, rtol=1e-12):
            msg = "all transfers must have the same sample_rate_hz."
            raise ValueError(msg)
        if transfer.response_kind != response_kind:
            msg = "all transfers must have the same response_kind."
            raise ValueError(msg)
    if not isinstance(normalize_to_dc, bool):
        msg = "normalize_to_dc must be a bool."
        raise ValueError(msg)

    nyquist_hz = sample_rate_hz / 2
    resolved_maximum_hz = (
        nyquist_hz if maximum_frequency_hz is None else maximum_frequency_hz
    )
    if (
        not np.isfinite(resolved_maximum_hz)
        or not 0 < resolved_maximum_hz <= nyquist_hz
    ):
        msg = "maximum_frequency_hz must be finite and in (0, Nyquist]."
        raise ValueError(msg)
    frequency_hz = np.linspace(0.0, resolved_maximum_hz, num_points)
    z_inverse = np.exp(-1j * 2 * np.pi * frequency_hz / sample_rate_hz)

    ax = _get_axes(ax)
    for label, transfer in transfers.items():
        numerator = np.polynomial.polynomial.polyval(z_inverse, transfer.numerator)
        denominator = np.polynomial.polynomial.polyval(
            z_inverse,
            transfer.denominator,
        )
        magnitude = np.abs(numerator / denominator)
        if normalize_to_dc:
            magnitude = magnitude / abs(transfer.dc_gain)
        magnitude_db = 20 * np.log10(np.maximum(magnitude, np.finfo(float).tiny))
        ax.plot(frequency_hz * 1e-9, magnitude_db, label=label)

    ax.set_xlabel("Frequency [GHz]")
    if normalize_to_dc:
        ax.set_ylabel("Magnitude relative to DC [dB]")
    elif response_kind == "transimpedance_ohm":
        ax.set_ylabel("Transimpedance magnitude [dBΩ re 1 Ω]")
    else:
        ax.set_ylabel("Magnitude [dB]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    return ax


def plot_statistical_ber_comparison(
    analyses: dict[str, StatisticalEyeAnalysis],
    ax: Axes | None = None,
) -> Axes:
    """Compare optimized Gaussian-mixture BER versus phase for RX cases."""
    if not analyses:
        msg = "analyses must contain at least one statistical eye result."
        raise ValueError(msg)

    ax = _get_axes(ax)
    for label, analysis in analyses.items():
        if not isinstance(label, str) or not label:
            msg = "analysis labels must be non-empty strings."
            raise ValueError(msg)
        if not isinstance(analysis, StatisticalEyeAnalysis):
            msg = "all analysis values must be StatisticalEyeAnalysis instances."
            raise ValueError(msg)
        phases = np.array(
            [result.sampling_phase_ui for result in analysis.phase_results],
            dtype=float,
        )
        ber_values = np.array(
            [
                max(result.ber, np.finfo(float).tiny)
                for result in analysis.phase_results
            ],
            dtype=float,
        )
        ax.semilogy(phases, ber_values, marker="o", label=label)

    ax.set_xlabel("Sampling phase [UI]")
    ax.set_ylabel("Optimized Gaussian-mixture BER")
    ax.set_xlim(0.0, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    return ax


def plot_jittered_ber_bathtub(
    analysis: JitteredEyeAnalysis,
    target_ber: float | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Plot fixed-threshold residual-jitter BER versus nominal phase."""
    if not isinstance(analysis, JitteredEyeAnalysis):
        msg = "analysis must be a JitteredEyeAnalysis."
        raise ValueError(msg)
    if target_ber is not None and (
        not np.isfinite(target_ber) or not 0 < target_ber <= 0.5
    ):
        msg = "target_ber must be finite and in (0, 0.5] when provided."
        raise ValueError(msg)

    display_ber = np.maximum(analysis.jitter_averaged_ber, np.finfo(float).tiny)
    optimum = analysis.optimum_result
    ax = _get_axes(ax)
    ax.semilogy(
        analysis.nominal_phases_ui,
        display_ber,
        marker="o",
        label="residual-jitter BER",
    )
    ax.semilogy(
        [optimum.nominal_sampling_phase_ui],
        [max(optimum.ber, np.finfo(float).tiny)],
        marker="*",
        markersize=11,
        color="C3",
        linestyle="none",
        label="optimum",
    )
    if target_ber is not None:
        ax.axhline(
            target_ber,
            color="black",
            linestyle="--",
            linewidth=1.0,
            label="target BER",
        )
    ax.set_xlabel("Nominal sampling phase [UI]")
    ax.set_ylabel("Jitter-averaged BER")
    ax.set_xlim(0.0, 1.0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    return ax


def _require_eye_metric_keys(eye_metrics: dict[str, float]) -> None:
    required_keys = {"min_one_level", "max_zero_level"}
    missing_keys = required_keys - eye_metrics.keys()
    if missing_keys:
        keys = ", ".join(sorted(missing_keys))
        msg = f"eye_metrics is missing required key(s): {keys}."
        raise ValueError(msg)
