"""Optical link-budget and NRZ/OOK level-diagram helpers."""

from dataclasses import dataclass
from math import isfinite
from typing import Sequence

from oma_ber.modulation import NRZLevels, nrz_levels_from_powers
from oma_ber.units import db_to_linear, watt_to_dbm


@dataclass(frozen=True)
class OpticalLoss:
    """Named passive optical power loss in decibels (dB)."""

    name: str
    loss_db: float

    def __post_init__(self) -> None:
        if not self.name.strip():
            msg = "name must not be empty."
            raise ValueError(msg)
        if not isfinite(self.loss_db) or self.loss_db < 0:
            msg = "loss_db must be finite and non-negative."
            raise ValueError(msg)


@dataclass(frozen=True)
class OpticalLevelPoint:
    """NRZ/OOK optical levels at one link reference plane.

    Optical powers and OMA are in W, OMA is also reported in dBm, and loss
    values are passive optical power losses in dB.
    """

    name: str
    incremental_loss_db: float
    cumulative_loss_db: float
    p0_w: float
    p1_w: float
    pavg_w: float
    oma_w: float
    oma_dbm: float
    er_linear: float
    er_db: float


def _level_point(
    name: str,
    levels: NRZLevels,
    incremental_loss_db: float,
    cumulative_loss_db: float,
) -> OpticalLevelPoint:
    return OpticalLevelPoint(
        name=name,
        incremental_loss_db=incremental_loss_db,
        cumulative_loss_db=cumulative_loss_db,
        p0_w=levels.p0_w,
        p1_w=levels.p1_w,
        pavg_w=levels.pavg_w,
        oma_w=levels.oma_w,
        oma_dbm=watt_to_dbm(levels.oma_w),
        er_linear=levels.er_linear,
        er_db=levels.er_db,
    )


def apply_optical_loss(levels: NRZLevels, loss_db: float) -> NRZLevels:
    """Apply a non-negative, level-independent optical power loss in dB."""
    if not isfinite(loss_db) or loss_db < 0:
        msg = "loss_db must be finite and non-negative."
        raise ValueError(msg)
    transmission = 1 / db_to_linear(loss_db)
    return nrz_levels_from_powers(
        p0_w=levels.p0_w * transmission,
        p1_w=levels.p1_w * transmission,
    )


def build_level_diagram(
    initial_levels: NRZLevels,
    losses: Sequence[OpticalLoss],
    initial_name: str = "Modulator output",
    margin_db: float = 0.0,
) -> list[OpticalLevelPoint]:
    """Build a passive optical level diagram from modulator output to PD input.

    margin_db is represented as a final design-margin point when positive. The
    returned first point is the unattenuated initial reference plane.
    """
    if not initial_name.strip():
        msg = "initial_name must not be empty."
        raise ValueError(msg)
    if not isfinite(margin_db) or margin_db < 0:
        msg = "margin_db must be finite and non-negative."
        raise ValueError(msg)
    if not all(isinstance(loss, OpticalLoss) for loss in losses):
        msg = "losses must contain only OpticalLoss values."
        raise ValueError(msg)

    points = [_level_point(initial_name, initial_levels, 0.0, 0.0)]
    current_levels = initial_levels
    cumulative_loss_db = 0.0
    all_losses = list(losses)
    if margin_db > 0:
        all_losses.append(OpticalLoss(name="Design margin", loss_db=margin_db))

    for loss in all_losses:
        current_levels = apply_optical_loss(current_levels, loss.loss_db)
        cumulative_loss_db += loss.loss_db
        points.append(
            _level_point(
                name=loss.name,
                levels=current_levels,
                incremental_loss_db=loss.loss_db,
                cumulative_loss_db=cumulative_loss_db,
            ),
        )
    return points


def pd_input_oma_dbm(
    tx_oma_dbm: float,
    losses_db: list[float] | tuple[float, ...],
    margin_db: float = 0.0,
) -> float:
    """Calculate photodiode-input OMA in dBm from TX OMA and dB losses.

    tx_oma_dbm is the transmitter OMA in dBm. losses_db and margin_db are
    power penalties in dB.
    """
    if not isinstance(losses_db, list | tuple):
        msg = "losses_db must be a list or tuple of dB losses."
        raise ValueError(msg)
    if any(loss_db < 0 for loss_db in losses_db):
        msg = "losses_db values must be non-negative."
        raise ValueError(msg)
    if margin_db < 0:
        msg = "margin_db must be non-negative."
        raise ValueError(msg)

    return tx_oma_dbm - sum(losses_db) - margin_db
