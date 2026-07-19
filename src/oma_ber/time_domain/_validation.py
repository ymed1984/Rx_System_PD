"""Shared validation helpers for time-domain analysis."""

import numpy as np


def validated_bits(bits: np.ndarray) -> np.ndarray:
    """Return a copied, read-only one-dimensional 0/1 bit array."""
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

    copied = np.array(bit_values, dtype=np.int_, copy=True)
    copied.setflags(write=False)
    return copied


def readonly_float_array(values: np.ndarray, name: str) -> np.ndarray:
    """Return a finite, copied, read-only one-dimensional float array."""
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        msg = f"{name} must be a one-dimensional array."
        raise ValueError(msg)
    if array.size == 0:
        msg = f"{name} must contain at least one sample."
        raise ValueError(msg)
    if not np.all(np.isfinite(array)):
        msg = f"{name} must contain only finite values."
        raise ValueError(msg)

    copied = np.array(array, dtype=float, copy=True)
    copied.setflags(write=False)
    return copied
