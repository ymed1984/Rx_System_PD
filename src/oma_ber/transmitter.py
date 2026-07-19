"""External-laser and NRZ/OOK modulator models."""

from dataclasses import dataclass
from math import isfinite

from oma_ber.modulation import NRZLevels, nrz_levels_from_powers
from oma_ber.units import db_to_linear, dbm_to_watt, watt_to_dbm


@dataclass(frozen=True)
class ExternalLaser:
    """External laser output power.

    output_power_dbm is continuous-wave optical power in dBm at the modulator
    input reference plane.
    """

    output_power_dbm: float

    def __post_init__(self) -> None:
        if not isfinite(self.output_power_dbm):
            msg = "output_power_dbm must be finite."
            raise ValueError(msg)

    @property
    def output_power_w(self) -> float:
        """Return continuous-wave laser output power in watts (W)."""
        return dbm_to_watt(self.output_power_dbm)


@dataclass(frozen=True)
class OOKModulator:
    """NRZ/OOK modulator parameters using one-level insertion loss.

    insertion_loss_db is the optical power loss from laser input to the optical
    1 level in dB. er_db is the output extinction ratio in dB.
    specified_oma_dbm may hold an independent measured/specification OMA for
    consistency reporting; it does not override the laser/IL/ER-derived levels.
    """

    insertion_loss_db: float
    er_db: float
    specified_oma_dbm: float | None = None

    def __post_init__(self) -> None:
        if not isfinite(self.insertion_loss_db) or self.insertion_loss_db < 0:
            msg = "insertion_loss_db must be finite and non-negative."
            raise ValueError(msg)
        if not isfinite(self.er_db) or self.er_db <= 0:
            msg = "er_db must be finite and positive."
            raise ValueError(msg)
        if self.specified_oma_dbm is not None and not isfinite(self.specified_oma_dbm):
            msg = "specified_oma_dbm must be finite when provided."
            raise ValueError(msg)


@dataclass(frozen=True)
class ModulatorOutput:
    """Calculated NRZ/OOK modulator output and optional OMA consistency error."""

    levels: NRZLevels
    calculated_oma_dbm: float
    specified_oma_dbm: float | None
    oma_error_db: float | None


def modulator_output_from_laser(
    laser: ExternalLaser,
    modulator: OOKModulator,
) -> ModulatorOutput:
    """Calculate modulator output levels using optical-1 insertion loss.

    The optical 1 level is laser power attenuated by insertion_loss_db. The
    optical 0 level follows from er_db. All calculated level powers are in W.
    """
    one_level_transmission = 1 / db_to_linear(modulator.insertion_loss_db)
    p1_w = laser.output_power_w * one_level_transmission
    er_linear = db_to_linear(modulator.er_db)
    p0_w = p1_w / er_linear
    levels = nrz_levels_from_powers(p0_w=p0_w, p1_w=p1_w)
    calculated_oma_dbm = watt_to_dbm(levels.oma_w)
    oma_error_db = (
        calculated_oma_dbm - modulator.specified_oma_dbm
        if modulator.specified_oma_dbm is not None
        else None
    )
    return ModulatorOutput(
        levels=levels,
        calculated_oma_dbm=calculated_oma_dbm,
        specified_oma_dbm=modulator.specified_oma_dbm,
        oma_error_db=oma_error_db,
    )
