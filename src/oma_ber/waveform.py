"""Deterministic NRZ/OOK waveform generation helpers."""

import numpy as np


def samples_per_symbol(sample_rate_hz: float, symbol_rate_baud: float) -> int:
    """Calculate integer samples per symbol from rates in Hz and baud.

    sample_rate_hz is the waveform sample rate in samples/s. symbol_rate_baud
    is the NRZ/OOK symbol rate in symbols/s. The returned value is rounded to
    the nearest integer and must be at least 2 samples/symbol.
    """
    if sample_rate_hz <= 0:
        msg = "sample_rate_hz must be positive."
        raise ValueError(msg)
    if symbol_rate_baud <= 0:
        msg = "symbol_rate_baud must be positive."
        raise ValueError(msg)

    ratio = sample_rate_hz / symbol_rate_baud
    if ratio < 2:
        msg = "sample_rate_hz / symbol_rate_baud must be at least 2 samples/symbol."
        raise ValueError(msg)
    return int(round(ratio))


def nrz_bits_to_levels(
    bits: np.ndarray,
    low_level: float,
    high_level: float,
    samples_per_symbol: int,
) -> np.ndarray:
    """Convert 0/1 NRZ/OOK bits to repeated low/high level samples.

    low_level and high_level may be optical power in W or current in A, but
    both levels must use the same physical unit. samples_per_symbol is an
    integer number of waveform samples per NRZ/OOK symbol.
    """
    if high_level <= low_level:
        msg = "high_level must be greater than low_level."
        raise ValueError(msg)
    if samples_per_symbol < 2:
        msg = "samples_per_symbol must be at least 2."
        raise ValueError(msg)

    bit_values = np.asarray(bits)
    if bit_values.ndim != 1:
        msg = "bits must be a one-dimensional array."
        raise ValueError(msg)
    if bit_values.size == 0:
        msg = "bits must contain at least one bit."
        raise ValueError(msg)
    if not np.all((bit_values == 0) | (bit_values == 1)):
        msg = "bits must contain only 0 and 1."
        raise ValueError(msg)

    levels = np.where(bit_values == 0, low_level, high_level)
    return np.repeat(levels.astype(float), samples_per_symbol)


def prbs_bits(num_bits: int, seed: int = 1) -> np.ndarray:
    """Generate a deterministic pseudo-random binary sequence.

    num_bits is the number of 0/1 bits to return. seed controls the deterministic
    random generator. This helper is for repeatable waveform examples and does
    not claim compliance with a standards PRBS polynomial.
    """
    if num_bits <= 0:
        msg = "num_bits must be positive."
        raise ValueError(msg)
    return np.random.default_rng(seed).integers(0, 2, size=num_bits, dtype=np.int_)
