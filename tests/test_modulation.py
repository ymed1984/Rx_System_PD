import pytest

from oma_ber.modulation import NRZLevels, nrz_levels_from_oma_er
from oma_ber.units import db_to_linear


def test_nrz_levels_from_oma_er_returns_dataclass() -> None:
    levels = nrz_levels_from_oma_er(oma_w=1e-3, er_db=6)

    assert isinstance(levels, NRZLevels)
    assert levels.oma_w == pytest.approx(1e-3)
    assert levels.er_db == pytest.approx(6)


def test_nrz_levels_preserve_oma_w() -> None:
    oma_w = 2e-3
    levels = nrz_levels_from_oma_er(oma_w=oma_w, er_db=4)

    assert levels.p1_w - levels.p0_w == pytest.approx(oma_w)


def test_nrz_levels_preserve_extinction_ratio() -> None:
    er_db = 6
    levels = nrz_levels_from_oma_er(oma_w=1e-3, er_db=er_db)

    assert levels.p1_w / levels.p0_w == pytest.approx(db_to_linear(er_db))
    assert levels.er_linear == pytest.approx(db_to_linear(er_db))


def test_nrz_levels_calculate_average_power_w() -> None:
    levels = nrz_levels_from_oma_er(oma_w=1e-3, er_db=3)

    assert levels.pavg_w == pytest.approx((levels.p0_w + levels.p1_w) / 2)


def test_nrz_levels_reject_non_positive_oma_w() -> None:
    with pytest.raises(ValueError, match="oma_w must be positive"):
        nrz_levels_from_oma_er(oma_w=0, er_db=3)


def test_nrz_levels_reject_er_not_larger_than_one() -> None:
    with pytest.raises(ValueError, match="extinction ratio larger than 1"):
        nrz_levels_from_oma_er(oma_w=1e-3, er_db=0)
