import numpy as np
import pytest

from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.sweep import calculate_ber_from_oma, required_oma_dbm, sweep_oma
from oma_ber.units import dbm_to_watt


def example_pd() -> Photodiode:
    return Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=1e-9,
        bandwidth_3db_hz=40e9,
    )


def example_rx() -> Receiver:
    return Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )


def test_calculate_ber_from_oma_returns_required_keys() -> None:
    result = calculate_ber_from_oma(
        oma_dbm=-10,
        er_db=6,
        pd=example_pd(),
        rx=example_rx(),
    )

    expected_keys = {
        "oma_dbm",
        "oma_w",
        "er_db",
        "p0_w",
        "p1_w",
        "pavg_w",
        "i0_a",
        "i1_a",
        "delta_i_a",
        "sigma0_a",
        "sigma1_a",
        "threshold_a",
        "q_rx",
        "ber",
        "ber_from_q",
    }

    assert expected_keys <= result.keys()


def test_calculate_ber_from_oma_preserves_units_and_current_delta() -> None:
    pd = example_pd()
    result = calculate_ber_from_oma(
        oma_dbm=-10,
        er_db=6,
        pd=pd,
        rx=example_rx(),
    )

    assert result["oma_w"] == pytest.approx(dbm_to_watt(-10))
    assert result["i1_a"] - result["i0_a"] == pytest.approx(result["delta_i_a"])
    assert result["delta_i_a"] == pytest.approx(pd.responsivity_a_per_w * result["oma_w"])
    assert result["threshold_a"] > result["i0_a"]
    assert result["threshold_a"] < result["i1_a"]


def test_sweep_oma_returns_one_result_per_oma_value() -> None:
    oma_values_dbm = np.array([-20.0, -15.0, -10.0])

    results = sweep_oma(
        oma_dbm_values=oma_values_dbm,
        er_db=6,
        pd=example_pd(),
        rx=example_rx(),
    )

    assert [result["oma_dbm"] for result in results] == pytest.approx(oma_values_dbm)


def test_ber_decreases_as_oma_increases() -> None:
    results = sweep_oma(
        oma_dbm_values=np.array([-25.0, -20.0, -15.0, -10.0]),
        er_db=6,
        pd=example_pd(),
        rx=example_rx(),
    )
    ber_values = [result["ber"] for result in results]

    assert ber_values == sorted(ber_values, reverse=True)
    assert len(set(ber_values)) == len(ber_values)


def test_required_oma_dbm_reaches_target_ber() -> None:
    pd = example_pd()
    rx = example_rx()
    target_ber = 1e-6

    oma_dbm = required_oma_dbm(
        target_ber=target_ber,
        er_db=6,
        pd=pd,
        rx=rx,
        oma_dbm_min=-40,
        oma_dbm_max=0,
    )
    result = calculate_ber_from_oma(oma_dbm, er_db=6, pd=pd, rx=rx)

    assert result["ber"] == pytest.approx(target_ber, rel=1e-4)


def test_required_oma_dbm_returns_minimum_if_min_already_meets_target() -> None:
    pd = example_pd()
    rx = example_rx()

    assert required_oma_dbm(
        target_ber=0.49,
        er_db=6,
        pd=pd,
        rx=rx,
        oma_dbm_min=-20,
        oma_dbm_max=0,
    ) == pytest.approx(-20)


def test_required_oma_dbm_fails_clearly_when_target_is_out_of_range() -> None:
    with pytest.raises(ValueError, match="cannot be reached"):
        required_oma_dbm(
            target_ber=1e-30,
            er_db=6,
            pd=example_pd(),
            rx=example_rx(),
            oma_dbm_min=-40,
            oma_dbm_max=-35,
        )


@pytest.mark.parametrize(
    ("target_ber", "match"),
    [
        (0.0, "target_ber must be between 0 and 0.5"),
        (0.5, "target_ber must be between 0 and 0.5"),
    ],
)
def test_required_oma_dbm_rejects_invalid_target_ber(target_ber: float, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        required_oma_dbm(
            target_ber=target_ber,
            er_db=6,
            pd=example_pd(),
            rx=example_rx(),
        )


def test_required_oma_dbm_rejects_invalid_search_range() -> None:
    with pytest.raises(ValueError, match="oma_dbm_min must be smaller"):
        required_oma_dbm(
            target_ber=1e-6,
            er_db=6,
            pd=example_pd(),
            rx=example_rx(),
            oma_dbm_min=0,
            oma_dbm_max=-10,
        )
