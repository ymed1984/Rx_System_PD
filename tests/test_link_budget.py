import pytest

from oma_ber.link_budget import pd_input_oma_dbm


def test_pd_input_oma_dbm_subtracts_losses_and_margin() -> None:
    assert pd_input_oma_dbm(
        tx_oma_dbm=-2.0,
        losses_db=[1.0, 2.5, 0.5],
        margin_db=1.0,
    ) == pytest.approx(-7.0)


def test_pd_input_oma_dbm_accepts_tuple_losses() -> None:
    assert pd_input_oma_dbm(
        tx_oma_dbm=0.0,
        losses_db=(1.0, 2.0),
    ) == pytest.approx(-3.0)


def test_pd_input_oma_dbm_accepts_empty_losses() -> None:
    assert pd_input_oma_dbm(
        tx_oma_dbm=-5.0,
        losses_db=[],
        margin_db=1.5,
    ) == pytest.approx(-6.5)


def test_pd_input_oma_dbm_rejects_non_sequence_losses() -> None:
    with pytest.raises(ValueError, match="losses_db must be a list or tuple"):
        pd_input_oma_dbm(tx_oma_dbm=0.0, losses_db=1.0)  # type: ignore[arg-type]


def test_pd_input_oma_dbm_rejects_negative_loss() -> None:
    with pytest.raises(ValueError, match="losses_db values must be non-negative"):
        pd_input_oma_dbm(tx_oma_dbm=0.0, losses_db=[1.0, -0.5])


def test_pd_input_oma_dbm_rejects_negative_margin() -> None:
    with pytest.raises(ValueError, match="margin_db must be non-negative"):
        pd_input_oma_dbm(tx_oma_dbm=0.0, losses_db=[1.0], margin_db=-0.1)
