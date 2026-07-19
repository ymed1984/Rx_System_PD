"""Receiver-boundary helpers for comparing photodiode specifications."""

from collections.abc import Mapping
from typing import Any

from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.sweep import calculate_ber_from_optical_levels, required_oma_dbm


def compare_photodiodes_at_optical_levels(
    p0_w: float,
    p1_w: float,
    photodiodes: Mapping[str, Photodiode],
    rx: Receiver,
    optimize_threshold: bool = True,
) -> dict[str, dict[str, Any]]:
    """Compare PD variants at the same receiver-boundary optical levels in W."""
    _validate_photodiodes(photodiodes)
    return {
        name: calculate_ber_from_optical_levels(
            p0_w=p0_w,
            p1_w=p1_w,
            pd=pd,
            rx=rx,
            optimize_threshold=optimize_threshold,
        )
        for name, pd in photodiodes.items()
    }


def required_oma_by_photodiode(
    target_ber: float,
    er_db: float,
    photodiodes: Mapping[str, Photodiode],
    rx: Receiver,
    oma_dbm_min: float = -40.0,
    oma_dbm_max: float = 10.0,
) -> dict[str, float]:
    """Compare required PD-input OMA in dBm for named PD variants."""
    _validate_photodiodes(photodiodes)
    return {
        name: required_oma_dbm(
            target_ber=target_ber,
            er_db=er_db,
            pd=pd,
            rx=rx,
            oma_dbm_min=oma_dbm_min,
            oma_dbm_max=oma_dbm_max,
        )
        for name, pd in photodiodes.items()
    }


def _validate_photodiodes(photodiodes: Mapping[str, Photodiode]) -> None:
    if not photodiodes:
        msg = "photodiodes must contain at least one named Photodiode."
        raise ValueError(msg)
    if any(not isinstance(name, str) or not name.strip() for name in photodiodes):
        msg = "photodiode names must not be empty."
        raise ValueError(msg)
    if any(not isinstance(pd, Photodiode) for pd in photodiodes.values()):
        msg = "photodiodes values must be Photodiode instances."
        raise ValueError(msg)
