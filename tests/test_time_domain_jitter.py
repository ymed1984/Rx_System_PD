from collections.abc import Callable
from dataclasses import replace

import numpy as np
import pytest

from oma_ber import Photodiode
from oma_ber.time_domain import (
    ResidualTimingJitter,
    TimeDomainNoiseModel,
    TimeGrid,
    analyze_jittered_tia_eye,
    analyze_statistical_tia_eye,
    calculate_tia_output_noise,
    first_order_lowpass_transfer,
    identity_transfer,
    prbs_bits,
    simulate_pd_tia_waveform,
)
from oma_ber.time_domain.jitter import _periodic_pattern_components


def test_zero_jitter_matches_phase2_statistical_eye() -> None:
    waveforms, noise = _bandlimited_receiver_case()
    phase2 = analyze_statistical_tia_eye(waveforms, noise)

    phase4 = analyze_jittered_tia_eye(
        waveforms,
        noise,
        ResidualTimingJitter(),
        threshold_grid_points=1025,
    )

    np.testing.assert_allclose(
        phase4.jitter_averaged_ber,
        np.array([result.ber for result in phase2.phase_results]),
        rtol=1e-10,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        phase4.optimized_thresholds_v,
        np.array([result.threshold_v for result in phase2.phase_results]),
        rtol=1e-10,
        atol=1e-15,
    )
    assert phase4.optimum_result.nominal_sampling_phase_ui == pytest.approx(
        phase2.optimum_result.sampling_phase_ui
    )


def test_small_random_jitter_does_not_change_flat_memoryless_mid_ui_sample() -> None:
    time_grid = TimeGrid(10e9, 16)
    bits = prbs_bits(7, 127)
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
    noise = calculate_tia_output_noise(waveforms, TimeDomainNoiseModel(8e-12))
    stationary_interpolated_noise = replace(
        noise,
        total_lag1_covariance_v2=noise.total_variance_v2,
    )

    no_jitter = analyze_jittered_tia_eye(
        waveforms,
        stationary_interpolated_noise,
        ResidualTimingJitter(),
        nominal_phases_ui=np.array([0.5]),
        threshold_grid_points=257,
    )
    random_jitter = analyze_jittered_tia_eye(
        waveforms,
        stationary_interpolated_noise,
        ResidualTimingJitter(random_jitter_rms_s=0.01 * time_grid.unit_interval_s),
        nominal_phases_ui=np.array([0.5]),
        random_quadrature_points=9,
        threshold_grid_points=257,
    )

    assert random_jitter.optimum_result.ber == pytest.approx(
        no_jitter.optimum_result.ber,
        rel=1e-10,
    )


def test_larger_random_jitter_worsens_best_bandlimited_ber() -> None:
    waveforms, noise = _bandlimited_receiver_case()

    low_jitter = analyze_jittered_tia_eye(
        waveforms,
        noise,
        ResidualTimingJitter(random_jitter_rms_s=0.01e-12),
        random_quadrature_points=9,
        threshold_grid_points=129,
    )
    high_jitter = analyze_jittered_tia_eye(
        waveforms,
        noise,
        ResidualTimingJitter(random_jitter_rms_s=3e-12),
        random_quadrature_points=9,
        threshold_grid_points=129,
    )

    assert high_jitter.optimum_result.ber > low_jitter.optimum_result.ber


def test_bandlimited_jitter_bathtub_worsens_at_both_ui_edges() -> None:
    waveforms, noise = _bandlimited_receiver_case()
    unit_interval_s = waveforms.time_grid.unit_interval_s

    analysis = analyze_jittered_tia_eye(
        waveforms,
        noise,
        ResidualTimingJitter(random_jitter_rms_s=2e-12),
        nominal_phases_ui=np.array([0.0, 0.5, 0.875]),
        decision_delay_s=0.375 * unit_interval_s,
        random_quadrature_points=9,
        threshold_grid_points=129,
    )

    assert analysis.optimum_result.nominal_sampling_phase_ui == pytest.approx(0.5)
    assert analysis.jitter_averaged_ber[0] > analysis.jitter_averaged_ber[1]
    assert analysis.jitter_averaged_ber[2] > analysis.jitter_averaged_ber[1]


def test_ui_boundary_crossing_keeps_intended_bit_labels() -> None:
    time_grid = TimeGrid(10e9, 8)
    bits = np.array([0, 1, 1, 0])
    pd = Photodiode(1.0, 1e-9)
    tia = identity_transfer(
        time_grid.sample_rate_hz,
        dc_gain=1.0,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits,
        1.0,
        3.0,
        time_grid,
        pd,
        tia_transimpedance_response=tia,
    )
    noise = calculate_tia_output_noise(waveforms, TimeDomainNoiseModel(1e-12))

    zero_means, one_means, _, _ = _periodic_pattern_components(
        waveforms,
        noise,
        sampling_phase_ui=1.0,
        decision_delay_s=0.0,
    )

    np.testing.assert_allclose(np.sort(zero_means), np.array([1.0, 3.0]) + 1e-9)
    np.testing.assert_allclose(np.sort(one_means), np.array([1.0, 3.0]) + 1e-9)


def test_random_jitter_quadrature_converges() -> None:
    waveforms, noise = _bandlimited_receiver_case()
    jitter = ResidualTimingJitter(random_jitter_rms_s=0.8e-12)

    nine_point = analyze_jittered_tia_eye(
        waveforms,
        noise,
        jitter,
        nominal_phases_ui=np.array([0.75]),
        random_quadrature_points=9,
        threshold_grid_points=129,
    )
    seventeen_point = analyze_jittered_tia_eye(
        waveforms,
        noise,
        jitter,
        nominal_phases_ui=np.array([0.75]),
        random_quadrature_points=17,
        threshold_grid_points=129,
    )

    assert nine_point.optimum_result.ber == pytest.approx(
        seventeen_point.optimum_result.ber,
        rel=1.5e-2,
    )


def test_sinusoidal_jitter_phase_quadrature_converges() -> None:
    waveforms, noise = _bandlimited_receiver_case()
    jitter = ResidualTimingJitter(sinusoidal_jitter_peak_s=1.5e-12)

    thirty_two_point = analyze_jittered_tia_eye(
        waveforms,
        noise,
        jitter,
        nominal_phases_ui=np.array([0.75]),
        sinusoidal_phase_points=32,
        threshold_grid_points=129,
    )
    reference = analyze_jittered_tia_eye(
        waveforms,
        noise,
        jitter,
        nominal_phases_ui=np.array([0.75]),
        sinusoidal_phase_points=128,
        threshold_grid_points=129,
    )

    assert thirty_two_point.optimum_result.ber == pytest.approx(
        reference.optimum_result.ber,
        rel=2e-3,
    )


def test_jittered_eye_supports_inverting_tia_and_reports_ui_units() -> None:
    time_grid = TimeGrid(10e9, 8)
    bits = prbs_bits(7, 63)
    pd = Photodiode(0.8, 1e-9)
    tia = identity_transfer(
        time_grid.sample_rate_hz,
        dc_gain=-1e3,
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
    noise = calculate_tia_output_noise(waveforms, TimeDomainNoiseModel(8e-12))
    jitter = ResidualTimingJitter(
        random_jitter_rms_s=1e-12,
        sinusoidal_jitter_peak_s=2e-12,
    )

    analysis = analyze_jittered_tia_eye(
        waveforms,
        noise,
        jitter,
        nominal_phases_ui=np.array([0.5]),
        random_quadrature_points=5,
        sinusoidal_phase_points=8,
        threshold_grid_points=129,
    )

    assert not analysis.optimum_result.one_is_high
    assert analysis.optimum_result.threshold_v < 0
    assert analysis.random_jitter_rms_ui == pytest.approx(0.01)
    assert analysis.sinusoidal_jitter_peak_ui == pytest.approx(0.02)
    assert analysis.optimum_result.quadrature_nodes == 40
    assert jitter.total_rms_s == pytest.approx(np.sqrt(3) * 1e-12)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (
            lambda: ResidualTimingJitter(random_jitter_rms_s=-1e-12),
            "random_jitter_rms_s must be finite and non-negative",
        ),
        (
            lambda: analyze_jittered_tia_eye(
                *_bandlimited_receiver_case(),
                ResidualTimingJitter(),
                random_quadrature_points=2,
            ),
            "random_quadrature_points must be in",
        ),
    ],
)
def test_jitter_model_rejects_invalid_inputs(
    call: Callable[[], object],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        call()


def _bandlimited_receiver_case():
    time_grid = TimeGrid(25e9, 8)
    bits = prbs_bits(7, 63)
    pd = Photodiode(0.8, 2e-9, shunt_resistance_ohm=1e8)
    pd_response = first_order_lowpass_transfer(time_grid.sample_rate_hz, 14e9)
    tia_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        18e9,
        dc_gain=1.5e3,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits,
        5e-6,
        20e-6,
        time_grid,
        pd,
        pd_response,
        tia_response,
    )
    noise = calculate_tia_output_noise(
        waveforms,
        TimeDomainNoiseModel(10e-12, rin_db_per_hz=-150.0),
    )
    return waveforms, noise
