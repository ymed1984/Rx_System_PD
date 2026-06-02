"""Matplotlib helpers for visualizing OMA-to-BER results."""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.ticker import LogFormatterSciNotation

from oma_ber.eye import eye_traces, eye_unit_interval_axis


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
    oma_dbm_values = np.array([result["oma_dbm"] for result in sweep_results], dtype=float)
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
    oma_dbm_values = np.array([result["oma_dbm"] for result in sweep_results], dtype=float)
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
        ax.axhline(eye_metrics["min_one_level"], color="C2", linestyle="--", linewidth=1.0)
        ax.axhline(eye_metrics["max_zero_level"], color="C3", linestyle="--", linewidth=1.0)
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


def _require_eye_metric_keys(eye_metrics: dict[str, float]) -> None:
    required_keys = {"min_one_level", "max_zero_level"}
    missing_keys = required_keys - eye_metrics.keys()
    if missing_keys:
        keys = ", ".join(sorted(missing_keys))
        msg = f"eye_metrics is missing required key(s): {keys}."
        raise ValueError(msg)
