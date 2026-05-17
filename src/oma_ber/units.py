"""Unit conversions for optical power and logarithmic ratios."""

from math import log10


def dbm_to_watt(dbm: float) -> float:
    """Convert optical power from dBm to watts (W)."""
    return 1e-3 * 10 ** (dbm / 10)


def watt_to_dbm(watt: float) -> float:
    """Convert optical power from watts (W) to dBm."""
    if watt <= 0:
        msg = "watt must be positive to convert to dBm."
        raise ValueError(msg)
    return 10 * log10(watt / 1e-3)


def db_to_linear(db: float) -> float:
    """Convert a power ratio from decibels (dB) to linear scale."""
    return 10 ** (db / 10)


def linear_to_db(value: float) -> float:
    """Convert a positive linear power ratio to decibels (dB)."""
    if value <= 0:
        msg = "value must be positive to convert to dB."
        raise ValueError(msg)
    return 10 * log10(value)


def rin_db_per_hz_to_linear(rin_db_per_hz: float) -> float:
    """Convert RIN from dB/Hz to linear 1/Hz."""
    return 10 ** (rin_db_per_hz / 10)
