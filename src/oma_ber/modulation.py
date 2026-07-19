"""Optical modulation level calculations for NRZ/OOK signals."""

from dataclasses import dataclass
from math import isfinite

from oma_ber.units import db_to_linear, linear_to_db


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


def nrz_levels_from_powers(p0_w: float, p1_w: float) -> NRZLevels:
    """Create NRZ/OOK optical levels from 0/1 powers in watts (W).

    p0_w may be zero for ideal extinction. p1_w must be strictly larger than
    p0_w. In the ideal-extinction case, er_linear and er_db are infinite.
    """
    if not isfinite(p0_w) or p0_w < 0:
        msg = "p0_w must be finite and non-negative."
        raise ValueError(msg)
    if not isfinite(p1_w) or p1_w <= p0_w:
        msg = "p1_w must be finite and larger than p0_w."
        raise ValueError(msg)

    oma_w = p1_w - p0_w
    pavg_w = (p0_w + p1_w) / 2
    if p0_w == 0:
        er_linear = float("inf")
        er_db = float("inf")
    else:
        er_linear = p1_w / p0_w
        er_db = linear_to_db(er_linear)

    return NRZLevels(
        p0_w=p0_w,
        p1_w=p1_w,
        pavg_w=pavg_w,
        oma_w=oma_w,
        er_linear=er_linear,
        er_db=er_db,
    )
