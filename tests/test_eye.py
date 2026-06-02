from collections.abc import Callable

import numpy as np
import pytest

from oma_ber.eye import eye_traces, eye_unit_interval_axis


def test_eye_traces_slices_overlapping_unit_intervals() -> None:
    waveform = np.arange(16, dtype=float)

    traces = eye_traces(waveform, samples_per_symbol=4, symbols_per_trace=2)

    assert traces.shape == (3, 8)
    np.testing.assert_allclose(traces[0], np.arange(0, 8))
    np.testing.assert_allclose(traces[1], np.arange(4, 12))
    np.testing.assert_allclose(traces[2], np.arange(8, 16))


def test_eye_traces_accepts_start_symbol() -> None:
    waveform = np.arange(20, dtype=float)

    traces = eye_traces(waveform, samples_per_symbol=4, symbols_per_trace=2, start_symbol=1)

    assert traces.shape == (3, 8)
    np.testing.assert_allclose(traces[0], np.arange(4, 12))


def test_eye_unit_interval_axis_uses_ui_units() -> None:
    axis = eye_unit_interval_axis(samples_per_symbol=4, symbols_per_trace=2)

    np.testing.assert_allclose(axis, np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75]))


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: eye_traces(np.array([[1.0]]), 2), "waveform must be a one-dimensional array"),
        (lambda: eye_traces(np.array([1.0, 2.0]), 1), "samples_per_symbol must be at least 2"),
        (lambda: eye_traces(np.array([1.0, 2.0]), 2, symbols_per_trace=0), "symbols_per_trace must be at least 1"),
        (lambda: eye_traces(np.array([1.0, 2.0]), 2, start_symbol=-1), "start_symbol must be non-negative"),
        (
            lambda: eye_traces(np.array([1.0, 2.0, 3.0]), 2, symbols_per_trace=2),
            "waveform must contain at least one complete eye trace",
        ),
        (lambda: eye_unit_interval_axis(1), "samples_per_symbol must be at least 2"),
        (lambda: eye_unit_interval_axis(2, symbols_per_trace=0), "symbols_per_trace must be at least 1"),
    ],
)
def test_eye_helpers_reject_invalid_inputs(call: Callable[[], object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
