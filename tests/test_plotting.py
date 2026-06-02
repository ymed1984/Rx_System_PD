import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import numpy as np
import pytest

from oma_ber.plotting import (
    plot_eye_diagram,
    plot_oma_sweep,
    plot_q_vs_oma,
    plot_required_oma_vs_responsivity,
)


def sample_sweep_results() -> list[dict]:
    return [
        {"oma_dbm": -20.0, "ber": 1e-3, "q_rx": 3.0},
        {"oma_dbm": -18.0, "ber": 1e-5, "q_rx": 4.3},
        {"oma_dbm": -16.0, "ber": 1e-8, "q_rx": 5.6},
    ]


def test_plot_oma_sweep_returns_axes_and_labels_units() -> None:
    ax = plot_oma_sweep(sample_sweep_results())

    assert isinstance(ax, Axes)
    assert "OMA" in ax.get_xlabel()
    assert "dBm" in ax.get_xlabel()
    assert "BER" in ax.get_ylabel()
    assert "1e-12" in ax.get_ylabel()
    assert ax.get_yscale() == "log"
    assert ax.get_ylim()[0] == pytest.approx(1e-12)
    plt.close(ax.figure)


def test_plot_q_vs_oma_returns_axes_and_labels_units() -> None:
    ax = plot_q_vs_oma(sample_sweep_results())

    assert isinstance(ax, Axes)
    assert "OMA" in ax.get_xlabel()
    assert "dBm" in ax.get_xlabel()
    assert ax.get_ylabel() == "Q"
    plt.close(ax.figure)


def test_plot_required_oma_vs_responsivity_returns_axes_and_labels_units() -> None:
    ax = plot_required_oma_vs_responsivity(
        responsivity_a_per_w_values=np.array([0.6, 0.8, 1.0]),
        required_oma_dbm_values=np.array([-14.0, -16.0, -17.5]),
    )

    assert isinstance(ax, Axes)
    assert "Responsivity" in ax.get_xlabel()
    assert "A/W" in ax.get_xlabel()
    assert "Required OMA" in ax.get_ylabel()
    assert "dBm" in ax.get_ylabel()
    plt.close(ax.figure)


def test_plot_eye_diagram_returns_axes_and_labels_ui_axis() -> None:
    waveform = np.array([0.0, 0.0, 1.0, 1.0] * 8)

    ax = plot_eye_diagram(waveform, samples_per_symbol=4)

    assert isinstance(ax, Axes)
    assert "UI" in ax.get_xlabel()
    assert ax.get_ylabel() == "Amplitude"
    assert len(ax.lines) > 0
    plt.close(ax.figure)


def test_plot_eye_diagram_accepts_sample_rate_and_metrics() -> None:
    waveform = np.array([0.0, 0.0, 1.0, 1.0] * 8)
    metrics = {"min_one_level": 0.8, "max_zero_level": 0.2}

    ax = plot_eye_diagram(
        waveform,
        samples_per_symbol=4,
        sample_rate_hz=100e9,
        eye_metrics=metrics,
        ylabel="Optical power [W]",
    )

    assert "ps" in ax.get_xlabel()
    assert ax.get_ylabel() == "Optical power [W]"
    assert len(ax.lines) > 2
    plt.close(ax.figure)


def test_plot_helpers_use_provided_axes() -> None:
    _, ax = plt.subplots()

    returned_ax = plot_oma_sweep(sample_sweep_results(), ax=ax)

    assert returned_ax is ax
    plt.close(ax.figure)


def test_plot_oma_sweep_rejects_missing_required_key() -> None:
    with pytest.raises(ValueError, match="missing required key"):
        plot_oma_sweep([{"oma_dbm": -20.0, "q_rx": 3.0}])


def test_plot_q_vs_oma_rejects_missing_required_key() -> None:
    with pytest.raises(ValueError, match="missing required key"):
        plot_q_vs_oma([{"oma_dbm": -20.0, "ber": 1e-3}])


def test_plot_oma_sweep_rejects_non_positive_ber_for_log_axis() -> None:
    with pytest.raises(ValueError, match="BER values must be positive"):
        plot_oma_sweep([{"oma_dbm": -20.0, "ber": 0.0}])


def test_plot_required_oma_vs_responsivity_rejects_length_mismatch() -> None:
    with pytest.raises(ValueError, match="same length"):
        plot_required_oma_vs_responsivity(
            responsivity_a_per_w_values=np.array([0.6, 0.8]),
            required_oma_dbm_values=np.array([-14.0]),
        )


def test_plot_eye_diagram_rejects_invalid_sample_rate() -> None:
    with pytest.raises(ValueError, match="sample_rate_hz must be positive"):
        plot_eye_diagram(np.ones(8), samples_per_symbol=4, sample_rate_hz=0.0)


def test_plot_eye_diagram_rejects_invalid_max_traces() -> None:
    with pytest.raises(ValueError, match="max_traces must be positive"):
        plot_eye_diagram(np.ones(8), samples_per_symbol=4, max_traces=0)


def test_plot_eye_diagram_rejects_missing_eye_metric_key() -> None:
    with pytest.raises(ValueError, match="eye_metrics is missing"):
        plot_eye_diagram(np.ones(8), samples_per_symbol=4, eye_metrics={"min_one_level": 0.8})
