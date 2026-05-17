"""Optical modulation level calculations for NRZ/OOK signals."""

from dataclasses import dataclass

from oma_ber.units import db_to_linear


@dataclass(frozen=True)
class NRZLevels:
    """Optical NRZ/OOK levels in watts (W) for a given OMA and ER."""

    p0_w: float
    p1_w: float
    pavg_w: float
    oma_w: float
    er_linear: float
    er_db: float


def nrz_levels_from_oma_er(oma_w: float, er_db: float) -> NRZLevels:
    """Calculate NRZ/OOK optical levels from OMA in watts (W) and ER in dB."""
    if oma_w <= 0:
        msg = "oma_w must be positive."
        raise ValueError(msg)

    er_linear = db_to_linear(er_db)
    if er_linear <= 1:
        msg = "er_db must correspond to an extinction ratio larger than 1."
        raise ValueError(msg)

    p0_w = oma_w / (er_linear - 1)
    p1_w = er_linear * p0_w
    pavg_w = (p0_w + p1_w) / 2

    return NRZLevels(
        p0_w=p0_w,
        p1_w=p1_w,
        pavg_w=pavg_w,
        oma_w=oma_w,
        er_linear=er_linear,
        er_db=er_db,
    )
