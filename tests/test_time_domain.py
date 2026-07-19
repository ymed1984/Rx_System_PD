from collections.abc import Callable
from math import pi

import numpy as np
import pytest
from scipy.signal import freqz

from oma_ber import Photodiode
from oma_ber.time_domain import (
    TimeGrid,
    analyze_deterministic_eye,
    first_order_lowpass_transfer,
    identity_transfer,
    nrz_optical_power_waveform,
    photodetect_power_waveform,
    prbs_bits,
    sample_eye_at_phase,
    simulate_pd_tia_waveform,
)


def test_time_grid_derives_exact_sample_rate_without_rounding() -> None:
    time_grid = TimeGrid(symbol_rate_baud=25e9, samples_per_symbol=16)

    assert time_grid.unit_interval_s == pytest.approx(40e-12)
    assert time_grid.sample_rate_hz == pytest.approx(400e9)
    assert time_grid.sample_interval_s == pytest.approx(2.5e-12)
    np.testing.assert_allclose(
        time_grid.sample_times_s(3),
        np.array([0.0, 2.5e-12, 5e-12]),
    )


def test_prbs7_has_maximal_period_and_expected_balance() -> None:
    bits = prbs_bits(order=7, num_bits=2 * (2**7 - 1))

    np.testing.assert_array_equal(bits[:127], bits[127:])
    assert np.count_nonzero(bits[:127]) == 64
    assert np.count_nonzero(bits[:127] == 0) == 63


@pytest.mark.parametrize(("order", "second_period"), [(9, 511), (15, 32767)])
def test_supported_prbs_polynomials_repeat_after_maximal_period(
    order: int,
    second_period: int,
) -> None:
    state_length = order
    bits = prbs_bits(order=order, num_bits=second_period + state_length)

    np.testing.assert_array_equal(bits[:state_length], bits[second_period:])


def test_optical_waveform_and_photodetection_keep_physical_units() -> None:
    time_grid = TimeGrid(10e9, 4)
    optical = nrz_optical_power_waveform(
        bits=np.array([0, 1]),
        p0_w=1e-6,
        p1_w=5e-6,
        time_grid=time_grid,
    )
    pd = Photodiode(responsivity_a_per_w=0.8, dark_current_a=2e-9)

    currents = photodetect_power_waveform(optical, pd)

    np.testing.assert_allclose(optical.power_w[:4], 1e-6)
    np.testing.assert_allclose(optical.power_w[4:], 5e-6)
    np.testing.assert_allclose(currents.signal_current_a[:4], 0.8e-6)
    np.testing.assert_allclose(currents.signal_current_a[4:], 4e-6)
    np.testing.assert_allclose(
        currents.total_current_a,
        currents.signal_current_a + 2e-9,
    )


def test_photodetection_applies_static_pd_saturation_before_bandwidth() -> None:
    time_grid = TimeGrid(10e9, 4)
    optical = nrz_optical_power_waveform(
        bits=np.array([0, 1]),
        p0_w=1e-6,
        p1_w=2e-3,
        time_grid=time_grid,
    )
    pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=0.0,
        saturation_power_w=1e-3,
    )

    currents = photodetect_power_waveform(optical, pd)

    assert currents.signal_current_a[-1] < 0.8 * 2e-3
    assert currents.signal_current_a[-1] < 0.8 * 1e-3


def test_first_order_transfer_has_requested_dc_gain_and_digital_3db_point() -> None:
    sample_rate_hz = 400e9
    bandwidth_3db_hz = 20e9
    response = first_order_lowpass_transfer(
        sample_rate_hz=sample_rate_hz,
        bandwidth_3db_hz=bandwidth_3db_hz,
        dc_gain=1.5e3,
        response_kind="transimpedance_ohm",
    )

    frequency_hz, transfer = freqz(
        response.numerator,
        response.denominator,
        worN=np.array([0.0, bandwidth_3db_hz]),
        fs=sample_rate_hz,
    )

    np.testing.assert_allclose(frequency_hz, np.array([0.0, bandwidth_3db_hz]))
    assert abs(transfer[0]) == pytest.approx(1.5e3)
    assert abs(transfer[1]) == pytest.approx(1.5e3 / np.sqrt(2))


def test_flat_pd_tia_chain_matches_dc_physics() -> None:
    time_grid = TimeGrid(25e9, 8)
    bits = np.array([0, 1, 0, 1])
    pd = Photodiode(responsivity_a_per_w=0.75, dark_current_a=3e-9)
    tia_gain_ohm = 2e3
    tia = identity_transfer(
        sample_rate_hz=time_grid.sample_rate_hz,
        dc_gain=tia_gain_ohm,
        response_kind="transimpedance_ohm",
    )

    result = simulate_pd_tia_waveform(
        bits=bits,
        p0_w=2e-6,
        p1_w=8e-6,
        time_grid=time_grid,
        pd=pd,
        tia_transimpedance_response=tia,
    )

    expected_current_a = np.repeat(
        np.array([0.75 * 2e-6 + 3e-9, 0.75 * 8e-6 + 3e-9] * 2),
        time_grid.samples_per_symbol,
    )
    assert result.warmup_symbols == 0
    np.testing.assert_allclose(result.pd_output_current_a, expected_current_a)
    assert result.tia_output_voltage_v is not None
    np.testing.assert_allclose(
        result.tia_output_voltage_v,
        expected_current_a * tia_gain_ohm,
    )


def test_automatic_warmup_removes_zero_state_filter_transient() -> None:
    time_grid = TimeGrid(25e9, 16)
    bits = prbs_bits(7, 127)
    pd = Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9)
    pd_response = first_order_lowpass_transfer(time_grid.sample_rate_hz, 14e9)
    tia_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        18e9,
        dc_gain=1.5e3,
        response_kind="transimpedance_ohm",
    )

    automatic = simulate_pd_tia_waveform(
        bits,
        2e-5,
        8e-5,
        time_grid,
        pd,
        pd_response,
        tia_response,
    )
    longer = simulate_pd_tia_waveform(
        bits,
        2e-5,
        8e-5,
        time_grid,
        pd,
        pd_response,
        tia_response,
        warmup_symbols=automatic.warmup_symbols + 10,
    )

    assert automatic.warmup_symbols > 0
    np.testing.assert_allclose(
        automatic.pd_output_current_a,
        longer.pd_output_current_a,
        rtol=2e-8,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        automatic.tia_output_voltage_v,
        longer.tia_output_voltage_v,
        rtol=2e-8,
        atol=1e-12,
    )


def test_ideal_eye_selects_center_phase_and_exact_opening() -> None:
    time_grid = TimeGrid(25e9, 8)
    bits = np.array([0, 1, 1, 0, 1, 0])
    waveform = np.repeat(np.where(bits == 0, 1e-6, 5e-6), 8)

    analysis = analyze_deterministic_eye(waveform, bits, time_grid, "A")

    assert analysis.optimum_metrics.sampling_phase_ui == pytest.approx(0.5)
    assert analysis.optimum_metrics.vertical_eye_opening == pytest.approx(4e-6)
    assert analysis.eye_time_ui[0] == pytest.approx(-1.0)
    assert analysis.eye_time_ui[-1] == pytest.approx(1.0)
    assert analysis.eye_traces.shape[1] == 2 * time_grid.samples_per_symbol + 1


def test_eye_analysis_supports_inverting_tia_polarity() -> None:
    time_grid = TimeGrid(25e9, 8)
    bits = np.array([0, 1, 1, 0, 1, 0])
    voltage_v = np.repeat(np.where(bits == 0, 0.2, -0.8), 8)

    analysis = analyze_deterministic_eye(voltage_v, bits, time_grid, "V")

    assert analysis.optimum_metrics.one_is_high is False
    assert analysis.optimum_metrics.vertical_eye_opening == pytest.approx(1.0)
    assert analysis.optimum_metrics.lower_inner_level == pytest.approx(-0.8)
    assert analysis.optimum_metrics.upper_inner_level == pytest.approx(0.2)


def test_fractional_sampling_and_explicit_delay_preserve_bit_alignment() -> None:
    time_grid = TimeGrid(10e9, 4)
    bits = np.array([0, 1, 0, 1])
    # The waveform is delayed by exactly one UI relative to the bit labels.
    undelayed = np.repeat(bits.astype(float), 4)
    delayed = np.concatenate((np.zeros(4), undelayed[:-4]))

    samples, sampled_bits = sample_eye_at_phase(
        delayed,
        bits,
        time_grid,
        sampling_phase_ui=0.5,
        decision_delay_s=time_grid.unit_interval_s,
    )

    np.testing.assert_allclose(samples, sampled_bits)


def test_one_pole_noise_bandwidth_reference_value() -> None:
    # This test records the continuous-time ENBW reference that Phase 2 must
    # reproduce when PSD integration is added.
    bandwidth_3db_hz = 14e9
    expected_noise_bandwidth_hz = pi / 2 * bandwidth_3db_hz

    assert expected_noise_bandwidth_hz == pytest.approx(21.991148575e9)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: TimeGrid(0.0, 8), "symbol_rate_baud must be finite and positive"),
        (lambda: TimeGrid(25e9, 1), "samples_per_symbol must be at least 2"),
        (lambda: prbs_bits(8, 10), "order must be one of"),
        (lambda: prbs_bits(7, 10, [0] * 7), "all-zero state"),
        (
            lambda: first_order_lowpass_transfer(100e9, 50e9),
            "below the Nyquist frequency",
        ),
        (
            lambda: simulate_pd_tia_waveform(
                np.array([0, 1]),
                1e-6,
                2e-6,
                TimeGrid(25e9, 8),
                Photodiode(0.8, 0.0),
                pd_current_response=identity_transfer(
                    200e9,
                    response_kind="transimpedance_ohm",
                ),
            ),
            "pd_current_response must have",
        ),
        (
            lambda: analyze_deterministic_eye(
                np.ones(16),
                np.array([0, 1]),
                TimeGrid(25e9, 8),
                "W",  # type: ignore[arg-type]
            ),
            "amplitude_unit must be",
        ),
    ],
)
def test_time_domain_helpers_reject_invalid_inputs(
    call: Callable[[], object],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        call()
