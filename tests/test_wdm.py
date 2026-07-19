from collections.abc import Callable

import pytest

from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.transmitter import ExternalLaser, OOKModulator
from oma_ber.units import db_to_linear
from oma_ber.wdm import (
    OpticalFiberLink,
    WdmChannel,
    analyze_wdm_channel,
    analyze_wdm_link,
)


def example_channel(name: str = "ch_1310", wavelength_nm: float = 1310.0) -> WdmChannel:
    return WdmChannel(
        name=name,
        wavelength_nm=wavelength_nm,
        laser=ExternalLaser(output_power_dbm=0.0),
        modulator=OOKModulator(insertion_loss_db=3.0, er_db=6.0),
        laser_to_modulator_coupling_loss_db=1.0,
        awg_split_loss_db=2.0,
        tx_fiber_coupling_loss_db=1.0,
        rx_fiber_coupling_loss_db=1.2,
    )


def example_fiber() -> OpticalFiberLink:
    return OpticalFiberLink(
        length_km=10.0,
        attenuation_db_per_km=0.35,
        additional_loss_db=0.5,
    )


def example_pd() -> Photodiode:
    return Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9)


def example_rx() -> Receiver:
    return Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
    )


def test_fiber_loss_uses_distance_attenuation_and_fixed_loss() -> None:
    assert example_fiber().propagation_loss_db == pytest.approx(4.0)


def test_wdm_channel_follows_requested_reference_planes() -> None:
    result = analyze_wdm_channel(
        channel=example_channel(),
        fiber_link=example_fiber(),
        margin_db=0.3,
        pd=example_pd(),
        rx=example_rx(),
    )

    assert [point.name for point in result.level_diagram] == [
        "External Laser",
        "Laser-to-MOD coupling",
        "MOD output",
        "AWG split loss",
        "TX fiber coupling",
        "Optical fiber (10 km)",
        "RX fiber coupling",
        "Design margin",
        "PD input",
    ]
    assert [point.section for point in result.level_diagram] == [
        "TX",
        "TX",
        "TX",
        "TX",
        "TX",
        "OpticalFiberLink",
        "RX",
        "RX",
        "RX",
    ]
    assert result.level_diagram[0].cw_power_dbm == pytest.approx(0.0)
    assert result.level_diagram[1].cw_power_dbm == pytest.approx(-1.0)
    assert result.level_diagram[2].p1_w is not None
    assert result.pd_input.is_modulated


def test_wdm_channel_cumulative_loss_matches_pd_input_one_level() -> None:
    channel = example_channel()
    result = analyze_wdm_channel(
        channel=channel,
        fiber_link=example_fiber(),
        margin_db=0.3,
        pd=example_pd(),
        rx=example_rx(),
    )

    expected_cumulative_loss_db = 1.0 + 3.0 + 2.0 + 1.0 + 4.0 + 1.2 + 0.3
    expected_p1_w = channel.laser.output_power_w / db_to_linear(expected_cumulative_loss_db)
    assert result.pd_input.cumulative_loss_db == pytest.approx(expected_cumulative_loss_db)
    assert result.pd_input.p1_w == pytest.approx(expected_p1_w)
    assert result.pd_input.er_db == pytest.approx(channel.modulator.er_db)
    assert result.receiver_result["p1_w"] == pytest.approx(expected_p1_w)


def test_multiple_wavelengths_are_reported_by_unique_channel_name() -> None:
    wavelength_dependent_channel = WdmChannel(
        name="ch_1330",
        wavelength_nm=1330.0,
        laser=ExternalLaser(output_power_dbm=0.0),
        modulator=OOKModulator(insertion_loss_db=3.0, er_db=6.0),
        laser_to_modulator_coupling_loss_db=1.0,
        awg_split_loss_db=2.0,
        tx_fiber_coupling_loss_db=1.0,
        rx_fiber_coupling_loss_db=1.2,
        fiber_attenuation_db_per_km=0.30,
    )
    results = analyze_wdm_link(
        channels=[
            example_channel("ch_1310", 1310.0),
            wavelength_dependent_channel,
        ],
        fiber_link=example_fiber(),
        pd=example_pd(),
        rx=example_rx(),
    )

    assert list(results) == ["ch_1310", "ch_1330"]
    assert results["ch_1310"].pd_input.wavelength_nm == pytest.approx(1310.0)
    assert results["ch_1330"].pd_input.wavelength_nm == pytest.approx(1330.0)
    assert results["ch_1310"].fiber_propagation_loss_db == pytest.approx(4.0)
    assert results["ch_1330"].fiber_propagation_loss_db == pytest.approx(3.5)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (
            lambda: OpticalFiberLink(-1.0, 0.35),
            "length_km must be finite and non-negative",
        ),
        (
            lambda: OpticalFiberLink(10.0, -0.35),
            "attenuation_db_per_km must be finite and non-negative",
        ),
        (
            lambda: WdmChannel(
                name="",
                wavelength_nm=1310.0,
                laser=ExternalLaser(0.0),
                modulator=OOKModulator(3.0, 6.0),
                laser_to_modulator_coupling_loss_db=1.0,
                awg_split_loss_db=2.0,
                tx_fiber_coupling_loss_db=1.0,
                rx_fiber_coupling_loss_db=1.0,
            ),
            "name must not be empty",
        ),
        (
            lambda: WdmChannel(
                name="bad",
                wavelength_nm=0.0,
                laser=ExternalLaser(0.0),
                modulator=OOKModulator(3.0, 6.0),
                laser_to_modulator_coupling_loss_db=1.0,
                awg_split_loss_db=2.0,
                tx_fiber_coupling_loss_db=1.0,
                rx_fiber_coupling_loss_db=1.0,
            ),
            "wavelength_nm must be finite and positive",
        ),
        (
            lambda: analyze_wdm_link([], example_fiber(), example_pd(), example_rx()),
            "channels must contain at least one WdmChannel",
        ),
        (
            lambda: analyze_wdm_link(
                [example_channel("same"), example_channel("same", 1330.0)],
                example_fiber(),
                example_pd(),
                example_rx(),
            ),
            "WDM channel names must be unique",
        ),
    ],
)
def test_wdm_models_reject_invalid_inputs(call: Callable[[], object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
