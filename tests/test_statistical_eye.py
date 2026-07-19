from collections.abc import Callable

import numpy as np
import pytest

from oma_ber import Photodiode, ber_from_gaussian_levels
from oma_ber.time_domain import (
    TimeDomainNoiseModel,
    TimeGrid,
    analyze_statistical_tia_eye,
    calculate_tia_output_noise,
    first_order_lowpass_transfer,
    gaussian_mixture_ber_for_threshold,
    identity_transfer,
    optimize_gaussian_mixture_threshold,
    prbs_bits,
    simulate_pd_tia_waveform,
)


def test_single_component_gaussian_mixture_matches_scalar_ber() -> None:
    zero_means_v = np.array([0.0])
    one_means_v = np.array([1.0])
    sigma_v = np.array([0.1])

    threshold_v, ber = optimize_gaussian_mixture_threshold(
        zero_means_v,
        one_means_v,
        sigma_v,
        sigma_v,
        one_is_high=True,
    )
    scalar = ber_from_gaussian_levels(0.0, 1.0, 0.1, 0.1)

    assert threshold_v == pytest.approx(scalar["threshold_a"])
    assert ber == pytest.approx(scalar["ber"])


def test_gaussian_mixture_supports_inverted_voltage_polarity() -> None:
    zero_means_v = np.array([1.0])
    one_means_v = np.array([0.0])
    sigma_v = np.array([0.1])

    threshold_v, ber = optimize_gaussian_mixture_threshold(
        zero_means_v,
        one_means_v,
        sigma_v,
        sigma_v,
        one_is_high=False,
    )

    assert threshold_v == pytest.approx(0.5)
    assert ber == pytest.approx(2.866515718791947e-7)
    assert gaussian_mixture_ber_for_threshold(
        zero_means_v,
        one_means_v,
        sigma_v,
        sigma_v,
        threshold_v,
        one_is_high=False,
    ) == pytest.approx(ber)


def test_gaussian_mixture_threshold_remains_stable_when_ber_underflows() -> None:
    sigma_v = np.array([1e-3])

    threshold_v, ber = optimize_gaussian_mixture_threshold(
        np.array([0.0]),
        np.array([1.0]),
        sigma_v,
        sigma_v,
        one_is_high=True,
    )

    assert threshold_v == pytest.approx(0.5, abs=1e-9)
    assert ber == 0.0


def test_flat_memoryless_statistical_eye_matches_scalar_level_model() -> None:
    time_grid = TimeGrid(10e9, 8)
    bits = prbs_bits(7, 127)
    pd = Photodiode(0.8, 1e-9)
    tia_gain_ohm = 1e3
    tia = identity_transfer(
        time_grid.sample_rate_hz,
        dc_gain=tia_gain_ohm,
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
    noise = calculate_tia_output_noise(
        waveforms,
        TimeDomainNoiseModel(8e-12),
    )

    analysis = analyze_statistical_tia_eye(waveforms, noise)
    optimum = analysis.optimum_result
    expected = ber_from_gaussian_levels(
        mu0_a=optimum.mean_zero_v,
        mu1_a=optimum.mean_one_v,
        sigma0_a=optimum.rms_zero_noise_v,
        sigma1_a=optimum.rms_one_noise_v,
    )

    assert optimum.sampling_phase_ui == pytest.approx(0.5)
    assert optimum.threshold_v == pytest.approx(expected["threshold_a"], rel=5e-9)
    assert optimum.ber == pytest.approx(expected["ber"], rel=1e-9)


def test_statistical_eye_density_integrates_to_unity_at_each_phase() -> None:
    waveforms, noise = _bandlimited_receiver_case()

    analysis = analyze_statistical_tia_eye(
        waveforms,
        noise,
        amplitude_points=801,
    )
    integrals = np.trapezoid(
        analysis.density_per_v,
        analysis.amplitude_axis_v,
        axis=1,
    )

    np.testing.assert_allclose(integrals, 1.0, rtol=2e-3, atol=2e-3)
    assert analysis.density_per_v.shape == (
        waveforms.time_grid.samples_per_symbol,
        801,
    )


def test_statistical_eye_selects_minimum_ber_phase() -> None:
    waveforms, noise = _bandlimited_receiver_case()

    analysis = analyze_statistical_tia_eye(waveforms, noise)

    assert analysis.optimum_result.ber == min(
        result.ber for result in analysis.phase_results
    )
    assert 0 <= analysis.optimum_result.sampling_phase_ui < 1
    assert np.isfinite(analysis.optimum_result.threshold_v)
    assert analysis.optimum_result.effective_q > 0


def test_fractional_noise_sampling_uses_lag_one_covariance() -> None:
    waveforms, noise = _bandlimited_receiver_case()
    half_sample_phase_ui = 0.5 / waveforms.time_grid.samples_per_symbol

    analysis = analyze_statistical_tia_eye(
        waveforms,
        noise,
        sampling_phases_ui=np.array([half_sample_phase_ui]),
    )
    result = analysis.optimum_result
    expected_first_variance = (
        0.25 * noise.total_variance_v2[0]
        + 0.25 * noise.total_variance_v2[1]
        + 0.5 * noise.total_lag1_covariance_v2[0]
    )
    first_bit = waveforms.bits[0]
    expected_sigma_v = np.sqrt(expected_first_variance)
    sampled_sigmas = result.zero_sigmas_v if first_bit == 0 else result.one_sigmas_v

    assert sampled_sigmas[0] == pytest.approx(expected_sigma_v)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (
            lambda: optimize_gaussian_mixture_threshold(
                np.array([0.0]),
                np.array([1.0]),
                np.array([0.0]),
                np.array([0.1]),
                True,
            ),
            "sigmas must be positive",
        ),
        (
            lambda: gaussian_mixture_ber_for_threshold(
                np.array([0.0]),
                np.array([1.0]),
                np.array([0.1]),
                np.array([0.1]),
                float("nan"),
                True,
            ),
            "threshold_v must be finite",
        ),
    ],
)
def test_statistical_eye_rejects_invalid_inputs(
    call: Callable[[], object],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        call()


def _bandlimited_receiver_case():
    time_grid = TimeGrid(25e9, 16)
    bits = prbs_bits(7, 127)
    pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=2e-9,
        shunt_resistance_ohm=1e8,
    )
    pd_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        14e9,
    )
    tia_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        18e9,
        dc_gain=1.5e3,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits,
        2e-5,
        8e-5,
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
