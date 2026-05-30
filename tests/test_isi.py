from collections.abc import Callable

import numpy as np
import pytest

from oma_ber.isi import (
    apply_lti_filter,
    first_order_lowpass_impulse_response,
    sample_at_symbol_centers,
    sampled_eye_levels,
)


def test_first_order_lowpass_impulse_response_has_unity_dc_gain() -> None:
    impulse_response = first_order_lowpass_impulse_response(
        sample_rate_hz=100e9,
        bandwidth_3db_hz=25e9,
        num_taps=32,
    )

    assert len(impulse_response) == 32
    assert np.sum(impulse_response) == pytest.approx(1.0)
    assert np.all(impulse_response >= 0)


def test_apply_lti_filter_identity_impulse_returns_waveform() -> None:
    waveform = np.array([0.0, 1.0, 1.0, 0.0])

    filtered = apply_lti_filter(waveform, np.array([1.0]))

    np.testing.assert_allclose(filtered, waveform)


def test_apply_lti_filter_preserves_constant_level_with_normalized_impulse() -> None:
    waveform = np.ones(16)
    impulse_response = np.array([0.5, 0.25, 0.25])

    filtered = apply_lti_filter(waveform, impulse_response)

    assert filtered[-1] == pytest.approx(1.0)


def test_sample_at_symbol_centers_uses_default_center_offset() -> None:
    waveform = np.arange(12, dtype=float)

    samples = sample_at_symbol_centers(waveform, samples_per_symbol=4)

    np.testing.assert_allclose(samples, np.array([2.0, 6.0, 10.0]))


def test_sample_at_symbol_centers_accepts_explicit_offset() -> None:
    waveform = np.arange(12, dtype=float)

    samples = sample_at_symbol_centers(waveform, samples_per_symbol=4, timing_offset_samples=1)

    np.testing.assert_allclose(samples, np.array([1.0, 5.0, 9.0]))


def test_sampled_eye_levels_reports_zero_penalty_for_ideal_levels() -> None:
    samples = np.array([0.0, 1.0, 0.0, 1.0])
    bits = np.array([0, 1, 0, 1])

    metrics = sampled_eye_levels(samples, bits)

    assert metrics["mean_zero_level"] == pytest.approx(0.0)
    assert metrics["mean_one_level"] == pytest.approx(1.0)
    assert metrics["vertical_eye_opening"] == pytest.approx(1.0)
    assert metrics["isi_penalty_db"] == pytest.approx(0.0)


def test_sampled_eye_levels_reports_positive_penalty_for_eye_closure() -> None:
    samples = np.array([0.0, 0.8, 0.2, 1.0])
    bits = np.array([0, 1, 0, 1])

    metrics = sampled_eye_levels(samples, bits)

    assert metrics["vertical_eye_opening"] == pytest.approx(0.6)
    assert metrics["mean_level_separation"] == pytest.approx(0.8)
    assert metrics["isi_penalty_db"] > 0.0


def test_sampled_eye_levels_reports_infinite_penalty_for_closed_eye() -> None:
    samples = np.array([0.4, 0.5, 0.6, 0.5])
    bits = np.array([0, 1, 0, 1])

    metrics = sampled_eye_levels(samples, bits)

    assert metrics["vertical_eye_opening"] < 0.0
    assert np.isinf(metrics["isi_penalty_db"])


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (
            lambda: first_order_lowpass_impulse_response(0.0, 25e9, 8),
            "sample_rate_hz must be positive",
        ),
        (
            lambda: first_order_lowpass_impulse_response(100e9, 0.0, 8),
            "bandwidth_3db_hz must be positive",
        ),
        (
            lambda: first_order_lowpass_impulse_response(100e9, 25e9, 0),
            "num_taps must be positive",
        ),
        (lambda: apply_lti_filter(np.array([[1.0]]), np.array([1.0])), "waveform must be a one-dimensional array"),
        (
            lambda: apply_lti_filter(np.array([1.0]), np.array([[1.0]])),
            "impulse_response must be a one-dimensional array",
        ),
        (lambda: apply_lti_filter(np.array([]), np.array([1.0])), "waveform must contain at least one sample"),
        (
            lambda: apply_lti_filter(np.array([1.0]), np.array([])),
            "impulse_response must contain at least one sample",
        ),
        (
            lambda: sample_at_symbol_centers(np.array([[1.0]]), 2),
            "waveform must be a one-dimensional array",
        ),
        (
            lambda: sample_at_symbol_centers(np.array([1.0, 2.0]), 1),
            "samples_per_symbol must be at least 2",
        ),
        (
            lambda: sample_at_symbol_centers(np.array([1.0]), 2),
            "waveform must contain at least one complete symbol",
        ),
        (
            lambda: sample_at_symbol_centers(np.array([1.0, 2.0]), 2, timing_offset_samples=2),
            "timing_offset_samples must be in",
        ),
        (
            lambda: sampled_eye_levels(np.array([[0.0]]), np.array([0])),
            "samples must be a one-dimensional array",
        ),
        (
            lambda: sampled_eye_levels(np.array([0.0]), np.array([[0]])),
            "bits must be a one-dimensional array",
        ),
        (
            lambda: sampled_eye_levels(np.array([]), np.array([])),
            "samples must contain at least one sample",
        ),
        (
            lambda: sampled_eye_levels(np.array([0.0]), np.array([0, 1])),
            "samples length must match bits length",
        ),
        (
            lambda: sampled_eye_levels(np.array([0.0, 1.0]), np.array([0, 2])),
            "bits must contain only 0 and 1",
        ),
        (
            lambda: sampled_eye_levels(np.array([0.0, 0.0]), np.array([0, 0])),
            "bits must contain at least one 0 and one 1",
        ),
    ],
)
def test_isi_helpers_reject_invalid_inputs(call: Callable[[], object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
