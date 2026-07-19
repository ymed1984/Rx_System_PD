from collections.abc import Callable

import pytest

from oma_ber.link_budget import OpticalLoss, apply_optical_loss, build_level_diagram
from oma_ber.modulation import nrz_levels_from_powers
from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.sweep import calculate_ber_from_optical_levels
from oma_ber.system import analyze_laser_to_receiver
from oma_ber.transmitter import ExternalLaser, OOKModulator


def example_pd() -> Photodiode:
    return Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9)


def example_rx() -> Receiver:
    return Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )


def test_optical_loss_preserves_er_and_attenuates_all_levels() -> None:
    levels = nrz_levels_from_powers(p0_w=1e-4, p1_w=4e-4)

    output = apply_optical_loss(levels, loss_db=3.0)

    transmission = 10 ** (-3.0 / 10)
    assert output.p0_w == pytest.approx(levels.p0_w * transmission)
    assert output.p1_w == pytest.approx(levels.p1_w * transmission)
    assert output.oma_w == pytest.approx(levels.oma_w * transmission)
    assert output.er_linear == pytest.approx(levels.er_linear)


def test_level_diagram_reports_named_and_cumulative_losses() -> None:
    levels = nrz_levels_from_powers(p0_w=1e-4, p1_w=4e-4)

    points = build_level_diagram(
        initial_levels=levels,
        losses=[OpticalLoss("Fiber coupling", 1.0), OpticalLoss("Waveguide", 2.0)],
        margin_db=0.5,
    )

    assert [point.name for point in points] == [
        "Modulator output",
        "Fiber coupling",
        "Waveguide",
        "Design margin",
    ]
    assert [point.cumulative_loss_db for point in points] == pytest.approx([0.0, 1.0, 3.0, 3.5])
    assert points[-1].oma_dbm == pytest.approx(points[0].oma_dbm - 3.5)


def test_integrated_analysis_connects_tx_link_and_receiver_boundaries() -> None:
    result = analyze_laser_to_receiver(
        laser=ExternalLaser(output_power_dbm=3.0),
        modulator=OOKModulator(insertion_loss_db=3.0, er_db=6.0),
        losses=[OpticalLoss("Fiber coupling", 1.5), OpticalLoss("Waveguide", 2.0)],
        margin_db=0.5,
        pd=example_pd(),
        rx=example_rx(),
    )

    assert result.level_diagram[0].name == "Modulator output"
    assert result.pd_input.name == "PD input"
    assert result.pd_input.cumulative_loss_db == pytest.approx(4.0)
    assert result.pd_input.oma_dbm == pytest.approx(
        result.modulator_output.calculated_oma_dbm - 4.0,
    )
    direct = calculate_ber_from_optical_levels(
        result.pd_input.p0_w,
        result.pd_input.p1_w,
        pd=example_pd(),
        rx=example_rx(),
    )
    assert result.receiver_result["ber"] == pytest.approx(direct["ber"])


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: OpticalLoss("", 1.0), "name must not be empty"),
        (lambda: OpticalLoss("Bad", -1.0), "loss_db must be finite and non-negative"),
        (
            lambda: apply_optical_loss(nrz_levels_from_powers(1e-4, 2e-4), -1.0),
            "loss_db must be finite and non-negative",
        ),
        (
            lambda: build_level_diagram(
                nrz_levels_from_powers(1e-4, 2e-4),
                [],
                margin_db=-1.0,
            ),
            "margin_db must be finite and non-negative",
        ),
    ],
)
def test_link_models_reject_invalid_inputs(call: Callable[[], object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
