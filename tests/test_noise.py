from collections.abc import Callable

import pytest

from oma_ber.noise import (
    K_B,
    Q_E,
    rin_noise_rms_a,
    shot_noise_rms_a,
    thermal_noise_rms_a,
    tia_noise_rms_a,
    total_noise_rms_a,
)
from oma_ber.photodiode import Photodiode, photocurrent_a
from oma_ber.receiver import Receiver


def test_photodiode_accepts_valid_parameters() -> None:
    pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=1e-9,
        bandwidth_3db_hz=40e9,
        capacitance_f=100e-15,
        saturation_power_w=10e-3,
        return_loss_db=20,
        bias_v=2,
        shunt_resistance_ohm=10e3,
        temperature_k=300,
    )

    assert pd.responsivity_a_per_w == pytest.approx(0.8)


def test_photocurrent_is_linear_without_saturation_power() -> None:
    pd = Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9)

    assert photocurrent_a(1e-3, pd) == pytest.approx(0.8e-3)


def test_photocurrent_uses_tanh_saturation_when_saturation_power_is_provided() -> None:
    pd = Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9, saturation_power_w=1e-3)

    low_power_current_a = photocurrent_a(1e-6, pd)
    high_power_current_a = photocurrent_a(100e-3, pd)

    assert low_power_current_a == pytest.approx(0.8e-6, rel=1e-3)
    assert high_power_current_a == pytest.approx(0.8e-3, rel=1e-2)


def test_photocurrent_rejects_negative_optical_power() -> None:
    pd = Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9)

    with pytest.raises(ValueError, match="optical_power_w must be non-negative"):
        photocurrent_a(-1e-6, pd)


@pytest.mark.parametrize(
    ("field_name", "field_value", "match"),
    [
        ("responsivity_a_per_w", 0.0, "responsivity_a_per_w must be positive"),
        ("dark_current_a", -1e-9, "dark_current_a must be non-negative"),
        ("bandwidth_3db_hz", 0.0, "bandwidth_3db_hz must be positive"),
        ("capacitance_f", 0.0, "capacitance_f must be positive"),
        ("saturation_power_w", 0.0, "saturation_power_w must be positive"),
        ("return_loss_db", 0.0, "return_loss_db must be positive"),
        ("bias_v", 0.0, "bias_v must be positive"),
        ("shunt_resistance_ohm", 0.0, "shunt_resistance_ohm must be positive"),
        ("temperature_k", 0.0, "temperature_k must be positive"),
    ],
)
def test_photodiode_rejects_invalid_physical_inputs(
    field_name: str,
    field_value: float,
    match: str,
) -> None:
    kwargs = {
        "responsivity_a_per_w": 0.8,
        "dark_current_a": 1e-9,
        field_name: field_value,
    }

    with pytest.raises(ValueError, match=match):
        Photodiode(**kwargs)


def test_receiver_accepts_valid_parameters() -> None:
    rx = Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )

    assert rx.noise_bandwidth_hz == pytest.approx(25e9)


@pytest.mark.parametrize(
    ("field_name", "field_value", "match"),
    [
        ("noise_bandwidth_hz", 0.0, "noise_bandwidth_hz must be positive"),
        (
            "input_current_noise_density_a_per_sqrt_hz",
            -1e-12,
            "input_current_noise_density_a_per_sqrt_hz must be non-negative",
        ),
    ],
)
def test_receiver_rejects_invalid_physical_inputs(
    field_name: str,
    field_value: float,
    match: str,
) -> None:
    kwargs = {
        "noise_bandwidth_hz": 25e9,
        "input_current_noise_density_a_per_sqrt_hz": 10e-12,
        field_name: field_value,
    }

    with pytest.raises(ValueError, match=match):
        Receiver(**kwargs)


def test_shot_noise_matches_equation() -> None:
    current_a = 1e-3
    dark_current_a = 1e-9
    bandwidth_hz = 25e9

    expected = (2 * Q_E * (current_a + dark_current_a) * bandwidth_hz) ** 0.5

    assert shot_noise_rms_a(current_a, dark_current_a, bandwidth_hz) == pytest.approx(expected)


def test_shot_noise_increases_with_current_a() -> None:
    low_noise_a = shot_noise_rms_a(current_a=1e-6, dark_current_a=1e-9, bandwidth_hz=25e9)
    high_noise_a = shot_noise_rms_a(current_a=1e-3, dark_current_a=1e-9, bandwidth_hz=25e9)

    assert high_noise_a > low_noise_a


def test_shot_noise_increases_with_bandwidth_hz() -> None:
    narrow_noise_a = shot_noise_rms_a(current_a=1e-3, dark_current_a=1e-9, bandwidth_hz=1e9)
    wide_noise_a = shot_noise_rms_a(current_a=1e-3, dark_current_a=1e-9, bandwidth_hz=25e9)

    assert wide_noise_a > narrow_noise_a


def test_tia_noise_is_zero_if_density_is_zero() -> None:
    assert tia_noise_rms_a(0.0, bandwidth_hz=25e9) == 0.0


def test_tia_noise_matches_equation() -> None:
    assert tia_noise_rms_a(10e-12, bandwidth_hz=25e9) == pytest.approx(
        10e-12 * (25e9) ** 0.5,
    )


def test_thermal_noise_matches_equation() -> None:
    shunt_resistance_ohm = 10e3
    temperature_k = 300
    bandwidth_hz = 25e9
    expected = (4 * K_B * temperature_k * bandwidth_hz / shunt_resistance_ohm) ** 0.5

    assert thermal_noise_rms_a(shunt_resistance_ohm, temperature_k, bandwidth_hz) == pytest.approx(expected)


def test_thermal_noise_increases_with_temperature_k() -> None:
    cold_noise_a = thermal_noise_rms_a(10e3, temperature_k=250, bandwidth_hz=25e9)
    hot_noise_a = thermal_noise_rms_a(10e3, temperature_k=350, bandwidth_hz=25e9)

    assert hot_noise_a > cold_noise_a


def test_thermal_noise_decreases_with_shunt_resistance_ohm() -> None:
    low_resistance_noise_a = thermal_noise_rms_a(5e3, temperature_k=300, bandwidth_hz=25e9)
    high_resistance_noise_a = thermal_noise_rms_a(20e3, temperature_k=300, bandwidth_hz=25e9)

    assert high_resistance_noise_a < low_resistance_noise_a


def test_rin_noise_is_zero_if_rin_db_per_hz_is_none() -> None:
    assert rin_noise_rms_a(0.8, optical_power_w=1e-3, rin_db_per_hz=None, bandwidth_hz=25e9) == 0.0


def test_rin_noise_matches_equation() -> None:
    expected = 0.8 * 1e-3 * (1e-15 * 25e9) ** 0.5

    assert rin_noise_rms_a(0.8, optical_power_w=1e-3, rin_db_per_hz=-150, bandwidth_hz=25e9) == pytest.approx(
        expected,
    )


def test_total_noise_is_at_least_each_contribution() -> None:
    pd = Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9, shunt_resistance_ohm=10e3)
    rx = Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )
    optical_power_w = 1e-3
    photocurrent_a = pd.responsivity_a_per_w * optical_power_w

    total_a = total_noise_rms_a(optical_power_w, photocurrent_a, pd, rx)
    shot_a = shot_noise_rms_a(photocurrent_a, pd.dark_current_a, rx.noise_bandwidth_hz)
    tia_a = tia_noise_rms_a(rx.input_current_noise_density_a_per_sqrt_hz, rx.noise_bandwidth_hz)
    rin_a = rin_noise_rms_a(
        pd.responsivity_a_per_w,
        optical_power_w,
        rx.rin_db_per_hz,
        rx.noise_bandwidth_hz,
    )
    thermal_a = thermal_noise_rms_a(
        pd.shunt_resistance_ohm,
        pd.temperature_k,
        rx.noise_bandwidth_hz,
    )

    assert total_a >= shot_a
    assert total_a >= tia_a
    assert total_a >= rin_a
    assert total_a >= thermal_a


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: shot_noise_rms_a(-1e-9, 0.0, 1e9), "current_a must be non-negative"),
        (lambda: shot_noise_rms_a(0.0, -1e-9, 1e9), "dark_current_a must be non-negative"),
        (lambda: shot_noise_rms_a(0.0, 0.0, 0.0), "bandwidth_hz must be positive"),
        (lambda: tia_noise_rms_a(-1e-12, 1e9), "input_noise_density_a_per_sqrt_hz must be non-negative"),
        (lambda: tia_noise_rms_a(0.0, 0.0), "bandwidth_hz must be positive"),
        (lambda: thermal_noise_rms_a(0.0, 300, 1e9), "shunt_resistance_ohm must be positive"),
        (lambda: thermal_noise_rms_a(10e3, 0.0, 1e9), "temperature_k must be positive"),
        (lambda: thermal_noise_rms_a(10e3, 300, 0.0), "bandwidth_hz must be positive"),
        (lambda: rin_noise_rms_a(0.0, 1e-3, -150, 1e9), "responsivity_a_per_w must be positive"),
        (lambda: rin_noise_rms_a(0.8, -1e-3, -150, 1e9), "optical_power_w must be non-negative"),
        (lambda: rin_noise_rms_a(0.8, 1e-3, -150, 0.0), "bandwidth_hz must be positive"),
    ],
)
def test_noise_functions_reject_invalid_physical_inputs(call: Callable[[], float], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
