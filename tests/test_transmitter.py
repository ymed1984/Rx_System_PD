from collections.abc import Callable

import pytest

from oma_ber.transmitter import ExternalLaser, OOKModulator, modulator_output_from_laser
from oma_ber.units import db_to_linear, dbm_to_watt


def test_external_laser_reports_power_in_watts() -> None:
    laser = ExternalLaser(output_power_dbm=0.0)

    assert laser.output_power_w == pytest.approx(1e-3)


def test_modulator_uses_one_level_insertion_loss() -> None:
    laser = ExternalLaser(output_power_dbm=0.0)
    modulator = OOKModulator(insertion_loss_db=3.0, er_db=6.0)

    output = modulator_output_from_laser(laser, modulator)

    expected_p1_w = dbm_to_watt(0.0) / db_to_linear(3.0)
    expected_p0_w = expected_p1_w / db_to_linear(6.0)
    assert output.levels.p1_w == pytest.approx(expected_p1_w)
    assert output.levels.p0_w == pytest.approx(expected_p0_w)
    assert output.levels.oma_w == pytest.approx(expected_p1_w - expected_p0_w)


def test_specified_oma_is_reported_as_consistency_error() -> None:
    baseline = modulator_output_from_laser(
        ExternalLaser(output_power_dbm=0.0),
        OOKModulator(insertion_loss_db=3.0, er_db=6.0),
    )
    specified_oma_dbm = baseline.calculated_oma_dbm - 0.4
    compared = modulator_output_from_laser(
        ExternalLaser(output_power_dbm=0.0),
        OOKModulator(
            insertion_loss_db=3.0,
            er_db=6.0,
            specified_oma_dbm=specified_oma_dbm,
        ),
    )

    assert compared.specified_oma_dbm == pytest.approx(specified_oma_dbm)
    assert compared.oma_error_db == pytest.approx(0.4)
    assert compared.levels == baseline.levels


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: ExternalLaser(float("nan")), "output_power_dbm must be finite"),
        (
            lambda: OOKModulator(insertion_loss_db=-1.0, er_db=6.0),
            "insertion_loss_db must be finite and non-negative",
        ),
        (
            lambda: OOKModulator(insertion_loss_db=1.0, er_db=0.0),
            "er_db must be finite and positive",
        ),
        (
            lambda: OOKModulator(1.0, 6.0, specified_oma_dbm=float("inf")),
            "specified_oma_dbm must be finite",
        ),
    ],
)
def test_transmitter_models_reject_invalid_inputs(
    call: Callable[[], object],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        call()
