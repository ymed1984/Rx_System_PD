"""Standards-polynomial PRBS helpers for deterministic waveform analysis."""

from collections.abc import Sequence

import numpy as np


_PRBS_SECOND_TAP = {
    7: 6,
    9: 5,
    15: 14,
    31: 28,
}


def prbs_bits(
    order: int,
    num_bits: int,
    initial_state: Sequence[int] | None = None,
) -> np.ndarray:
    """Generate PRBS7/9/15/31 bits from a standard primitive polynomial.

    The supported polynomial form is x**order + x**second_tap + 1.  The
    default initial state is all ones.  initial_state must contain exactly
    order binary values and must not be the all-zero lock-up state.
    """
    if order not in _PRBS_SECOND_TAP:
        supported = ", ".join(str(value) for value in sorted(_PRBS_SECOND_TAP))
        msg = f"order must be one of: {supported}."
        raise ValueError(msg)
    if isinstance(num_bits, bool) or not isinstance(num_bits, int):
        msg = "num_bits must be an integer."
        raise ValueError(msg)
    if num_bits <= 0:
        msg = "num_bits must be positive."
        raise ValueError(msg)

    if initial_state is None:
        state = np.ones(order, dtype=np.int_)
    else:
        state = np.asarray(initial_state)
        if state.ndim != 1 or state.size != order:
            msg = "initial_state length must match order."
            raise ValueError(msg)
        if not np.all((state == 0) | (state == 1)):
            msg = "initial_state must contain only 0 and 1."
            raise ValueError(msg)
        if not np.any(state):
            msg = "initial_state must not be the all-zero state."
            raise ValueError(msg)
        state = np.array(state, dtype=np.int_, copy=True)

    second_tap = _PRBS_SECOND_TAP[order]
    sequence = np.empty(num_bits + order, dtype=np.int_)
    sequence[:order] = state
    for index in range(num_bits):
        sequence[index + order] = sequence[index + second_tap] ^ sequence[index]

    result = np.array(sequence[:num_bits], copy=True)
    result.setflags(write=False)
    return result
