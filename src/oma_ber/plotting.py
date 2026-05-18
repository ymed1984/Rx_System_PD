"""Matplotlib helpers for visualizing OMA-to-BER results."""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.ticker import LogFormatterSciNotation


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
