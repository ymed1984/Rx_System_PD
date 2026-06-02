"""Eye-diagram trace preparation helpers."""

import numpy as np


def eye_traces(
    waveform: np.ndarray,
    samples_per_symbol: int,
    symbols_per_trace: int = 2,
    start_symbol: int = 0,
) -> np.ndarray:
    """Slice a one-dimensional waveform into overlapping eye-diagram traces.

    waveform uses arbitrary but consistent amplitude units, such as optical
    power in W or current in A. samples_per_symbol is in samples/symbol.
    symbols_per_trace is the number of unit intervals per eye trace.
    """
    waveform_values = np.asarray(waveform, dtype=float)
    if waveform_values.ndim != 1:
        msg = "waveform must be a one-dimensional array."
        raise ValueError(msg)
    if samples_per_symbol < 2:
        msg = "samples_per_symbol must be at least 2."
        raise ValueError(msg)
    if symbols_per_trace < 1:
        msg = "symbols_per_trace must be at least 1."
        raise ValueError(msg)
    if start_symbol < 0:
        msg = "start_symbol must be non-negative."
        raise ValueError(msg)

    trace_length = samples_per_symbol * symbols_per_trace
    start_sample = start_symbol * samples_per_symbol
    available_samples = waveform_values.size - start_sample
    if available_samples < trace_length:
        msg = "waveform must contain at least one complete eye trace."
        raise ValueError(msg)

    num_traces = available_samples // samples_per_symbol - symbols_per_trace + 1
    return np.vstack(
        [
            waveform_values[
                start_sample + trace_index * samples_per_symbol : start_sample
                + trace_index * samples_per_symbol
                + trace_length
            ]
            for trace_index in range(num_traces)
        ],
    )


def eye_unit_interval_axis(samples_per_symbol: int, symbols_per_trace: int = 2) -> np.ndarray:
    """Create an x-axis in unit intervals for eye-diagram traces."""
    if samples_per_symbol < 2:
        msg = "samples_per_symbol must be at least 2."
        raise ValueError(msg)
    if symbols_per_trace < 1:
        msg = "symbols_per_trace must be at least 1."
        raise ValueError(msg)

    return np.arange(samples_per_symbol * symbols_per_trace, dtype=float) / samples_per_symbol
