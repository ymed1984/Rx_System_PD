import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import numpy as np
import pytest

from oma_ber.plotting import (
    plot_deterministic_eye_analysis,
    plot_eye_diagram,
    plot_jittered_ber_bathtub,
    plot_oma_sweep,
    plot_q_vs_oma,
    plot_required_oma_vs_responsivity,
    plot_statistical_ber_comparison,
    plot_statistical_ber_vs_phase,
    plot_statistical_eye_density,
    plot_transfer_magnitude_comparison,
)
from oma_ber.time_domain import (
    TimeDomainNoiseModel,
    TimeGrid,
    ResidualTimingJitter,
    analyze_jittered_tia_eye,
    analyze_deterministic_eye,
    analyze_statistical_tia_eye,
    calculate_tia_output_noise,
    identity_transfer,
    simulate_pd_tia_waveform,
)
from oma_ber import Photodiode


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


def test_plot_deterministic_eye_analysis_marks_clock_and_physical_unit() -> None:
    time_grid = TimeGrid(25e9, 8)
    bits = np.array([0, 1, 1, 0, 1, 0])
    waveform_a = np.repeat(np.where(bits == 0, 1e-6, 5e-6), 8)
    analysis = analyze_deterministic_eye(waveform_a, bits, time_grid, "A")

    ax = plot_deterministic_eye_analysis(analysis, max_traces=2)

    assert "decision" in ax.get_xlabel()
    assert "[A]" in ax.get_ylabel()
    assert ax.get_xlim() == pytest.approx((-1.0, 1.0))
    assert len(ax.lines) >= 5
    plt.close(ax.figure)


def test_statistical_eye_plot_helpers_show_density_and_ber_phase() -> None:
    time_grid = TimeGrid(10e9, 8)
    bits = np.array([0, 1, 1, 0, 1, 0] * 8)
    pd = Photodiode(0.8, 1e-9)
    tia = identity_transfer(
        time_grid.sample_rate_hz,
        dc_gain=1e3,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits,
        2e-6,
        8e-6,
        time_grid,
        pd,
        tia_transimpedance_response=tia,
    )
    noise = calculate_tia_output_noise(waveforms, TimeDomainNoiseModel(10e-12))
    analysis = analyze_statistical_tia_eye(waveforms, noise)

    density_ax = plot_statistical_eye_density(analysis)
    ber_ax = plot_statistical_ber_vs_phase(analysis)

    assert "BER-optimum" in density_ax.get_xlabel()
    assert "[V]" in density_ax.get_ylabel()
    assert len(density_ax.collections) > 0
    assert ber_ax.get_yscale() == "log"
    assert "BER" in ber_ax.get_ylabel()
    plt.close(density_ax.figure)
    plt.close(ber_ax.figure)


def test_phase3_plot_helpers_compare_transfer_and_statistical_ber() -> None:
    time_grid = TimeGrid(10e9, 8)
    bits = np.array([0, 1, 1, 0, 1, 0] * 8)
    pd = Photodiode(0.8, 1e-9)
    flat_pd = identity_transfer(time_grid.sample_rate_hz)
    tia = identity_transfer(
        time_grid.sample_rate_hz,
        dc_gain=1e3,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits,
        2e-6,
        8e-6,
        time_grid,
        pd,
        pd_current_response=flat_pd,
        tia_transimpedance_response=tia,
    )
    noise = calculate_tia_output_noise(waveforms, TimeDomainNoiseModel(10e-12))
    analysis = analyze_statistical_tia_eye(waveforms, noise)

    response_ax = plot_transfer_magnitude_comparison(
        {"baseline": flat_pd, "candidate": identity_transfer(time_grid.sample_rate_hz)},
        maximum_frequency_hz=30e9,
    )
    ber_ax = plot_statistical_ber_comparison(
        {"baseline": analysis, "candidate": analysis}
    )

    assert "Frequency [GHz]" == response_ax.get_xlabel()
    assert "relative to DC" in response_ax.get_ylabel()
    assert len(response_ax.lines) == 2
    assert ber_ax.get_yscale() == "log"
    assert len(ber_ax.lines) == 2
    plt.close(response_ax.figure)
    plt.close(ber_ax.figure)


def test_phase4_plot_helper_shows_jitter_bathtub() -> None:
    time_grid = TimeGrid(10e9, 8)
    bits = np.array([0, 1, 1, 0, 1, 0] * 8)
    pd = Photodiode(0.8, 1e-9)
    tia = identity_transfer(
        time_grid.sample_rate_hz,
        dc_gain=1e3,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits,
        2e-6,
        8e-6,
        time_grid,
        pd,
        tia_transimpedance_response=tia,
    )
    noise = calculate_tia_output_noise(waveforms, TimeDomainNoiseModel(10e-12))
    analysis = analyze_jittered_tia_eye(
        waveforms,
        noise,
        ResidualTimingJitter(),
        threshold_grid_points=129,
    )

    ax = plot_jittered_ber_bathtub(analysis, target_ber=0.2)

    assert ax.get_yscale() == "log"
    assert "Nominal" in ax.get_xlabel()
    assert "Jitter-averaged" in ax.get_ylabel()
    assert len(ax.lines) == 3
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
        plot_eye_diagram(
            np.ones(8), samples_per_symbol=4, eye_metrics={"min_one_level": 0.8}
        )
