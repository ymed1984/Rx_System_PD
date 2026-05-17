import pytest

from oma_ber.units import (
    db_to_linear,
    dbm_to_watt,
    linear_to_db,
    rin_db_per_hz_to_linear,
    watt_to_dbm,
)


def test_dbm_to_watt_reference_values() -> None:
    assert dbm_to_watt(0) == pytest.approx(1e-3)
    assert dbm_to_watt(-10) == pytest.approx(1e-4)


def test_db_to_linear_reference_value() -> None:
    assert db_to_linear(3) == pytest.approx(2, rel=0.004)


def test_dbm_watt_round_trip() -> None:
    for dbm in [-30.0, -10.0, 0.0, 3.0]:
        assert watt_to_dbm(dbm_to_watt(dbm)) == pytest.approx(dbm)


def test_db_linear_round_trip() -> None:
    for db in [-20.0, -3.0, 0.0, 10.0]:
        assert linear_to_db(db_to_linear(db)) == pytest.approx(db)


def test_rin_db_per_hz_to_linear() -> None:
    assert rin_db_per_hz_to_linear(-150) == pytest.approx(1e-15)


def test_watt_to_dbm_rejects_non_positive_watt() -> None:
    with pytest.raises(ValueError, match="watt must be positive"):
        watt_to_dbm(0)


def test_linear_to_db_rejects_non_positive_value() -> None:
    with pytest.raises(ValueError, match="value must be positive"):
        linear_to_db(0)
