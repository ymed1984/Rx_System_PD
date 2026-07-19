from collections.abc import Callable
from math import pi

import numpy as np
import pytest

from oma_ber import K_B, Q_E, Photodiode
from oma_ber.time_domain import (
    DiscreteTransferFunction,
    TimeDomainNoiseModel,
    TimeGrid,
    calculate_tia_output_noise,
    equivalent_noise_bandwidth_hz,
    first_order_lowpass_transfer,
    identity_transfer,
    simulate_pd_tia_waveform,
)
from oma_ber.units import rin_db_per_hz_to_linear


def test_identity_transfer_noise_bandwidth_is_nyquist() -> None:
    response = identity_transfer(sample_rate_hz=200e9)

    assert equivalent_noise_bandwidth_hz(response) == pytest.approx(100e9)


def test_high_oversampling_one_pole_noise_bandwidth_converges_to_analog_value() -> None:
    sample_rate_hz = 20e12
    bandwidth_3db_hz = 10e9
    response = first_order_lowpass_transfer(sample_rate_hz, bandwidth_3db_hz)

    calculated_hz = equivalent_noise_bandwidth_hz(response)

    assert calculated_hz == pytest.approx(pi / 2 * bandwidth_3db_hz, rel=3e-3)


def test_constant_level_noise_components_match_psd_times_enbw() -> None:
    time_grid = TimeGrid(25e9, 32)
    bits = np.ones(64, dtype=int)
    optical_power_w = 20e-6
    responsivity_a_per_w = 0.8
    dark_current_a = 2e-9
    shunt_resistance_ohm = 2e6
    temperature_k = 310.0
    tia_density_a_per_sqrt_hz = 9e-12
    rin_db_per_hz = -150.0
    tia_gain_ohm = 1.2e3
    pd = Photodiode(
        responsivity_a_per_w=responsivity_a_per_w,
        dark_current_a=dark_current_a,
        shunt_resistance_ohm=shunt_resistance_ohm,
        temperature_k=temperature_k,
        dark_current_fano_factor=1.3,
    )
    tia_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        18e9,
        dc_gain=tia_gain_ohm,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits=bits,
        p0_w=1e-6,
        p1_w=optical_power_w,
        time_grid=time_grid,
        pd=pd,
        tia_transimpedance_response=tia_response,
    )
    noise = calculate_tia_output_noise(
        waveforms,
        TimeDomainNoiseModel(
            tia_input_current_noise_density_a_per_sqrt_hz=tia_density_a_per_sqrt_hz,
            rin_db_per_hz=rin_db_per_hz,
        ),
    )
    bandwidth_hz = noise.tia_noise_bandwidth_hz
    gain_squared = tia_gain_ohm**2
    photocurrent_a = responsivity_a_per_w * optical_power_w

    expected_photo_v2 = 2 * Q_E * photocurrent_a * bandwidth_hz * gain_squared
    expected_dark_v2 = 2 * Q_E * 1.3 * dark_current_a * bandwidth_hz * gain_squared
    expected_tia_v2 = tia_density_a_per_sqrt_hz**2 * bandwidth_hz * gain_squared
    expected_rin_v2 = (
        photocurrent_a**2
        * rin_db_per_hz_to_linear(rin_db_per_hz)
        * bandwidth_hz
        * gain_squared
    )
    expected_thermal_v2 = (
        4 * K_B * temperature_k / shunt_resistance_ohm * bandwidth_hz * gain_squared
    )

    np.testing.assert_allclose(noise.photocurrent_shot_variance_v2, expected_photo_v2)
    np.testing.assert_allclose(noise.dark_current_shot_variance_v2, expected_dark_v2)
    np.testing.assert_allclose(noise.tia_variance_v2, expected_tia_v2)
    np.testing.assert_allclose(noise.rin_variance_v2, expected_rin_v2)
    np.testing.assert_allclose(noise.shunt_thermal_variance_v2, expected_thermal_v2)
    np.testing.assert_allclose(
        noise.total_variance_v2,
        noise.photocurrent_shot_variance_v2
        + noise.dark_current_shot_variance_v2
        + noise.tia_variance_v2
        + noise.rin_variance_v2
        + noise.shunt_thermal_variance_v2,
    )


def test_filtered_noise_has_lag_one_covariance() -> None:
    time_grid = TimeGrid(25e9, 16)
    bits = np.ones(32, dtype=int)
    pd = Photodiode(0.8, 1e-9)
    tia = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        15e9,
        dc_gain=1e3,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits,
        1e-6,
        10e-6,
        time_grid,
        pd,
        tia_transimpedance_response=tia,
    )

    noise = calculate_tia_output_noise(
        waveforms,
        TimeDomainNoiseModel(10e-12),
    )

    assert np.all(noise.total_lag1_covariance_v2 > 0)
    assert np.all(noise.total_lag1_covariance_v2 < noise.total_variance_v2)


def test_static_saturation_reduces_rin_small_signal_conversion() -> None:
    time_grid = TimeGrid(10e9, 32)
    bits = np.ones(32, dtype=int)
    tia = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        10e9,
        dc_gain=1e3,
        response_kind="transimpedance_ohm",
    )
    linear_pd = Photodiode(0.8, 0.0)
    saturated_pd = Photodiode(0.8, 0.0, saturation_power_w=1e-3)

    linear_waveforms = simulate_pd_tia_waveform(
        bits,
        1e-6,
        2e-3,
        time_grid,
        linear_pd,
        tia_transimpedance_response=tia,
    )
    saturated_waveforms = simulate_pd_tia_waveform(
        bits,
        1e-6,
        2e-3,
        time_grid,
        saturated_pd,
        tia_transimpedance_response=tia,
    )
    model = TimeDomainNoiseModel(0.0, rin_db_per_hz=-140.0)

    linear_noise = calculate_tia_output_noise(linear_waveforms, model)
    saturated_noise = calculate_tia_output_noise(saturated_waveforms, model)

    assert np.all(saturated_noise.rin_variance_v2 < linear_noise.rin_variance_v2)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (
            lambda: TimeDomainNoiseModel(-1e-12),
            "must be finite and non-negative",
        ),
        (
            lambda: TimeDomainNoiseModel(1e-12, rin_db_per_hz=float("nan")),
            "rin_db_per_hz must be finite",
        ),
        (
            lambda: equivalent_noise_bandwidth_hz(
                # H(z)=1-z^-1 has zero DC gain.
                DiscreteTransferFunction(
                    np.array([1.0, -1.0]),
                    np.array([1.0]),
                    100e9,
                    "dimensionless",
                ),
            ),
            "non-zero DC gain",
        ),
    ],
)
def test_time_domain_noise_rejects_invalid_inputs(
    call: Callable[[], object],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        call()
