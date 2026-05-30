from collections.abc import Callable

import numpy as np
import pytest

from oma_ber.saturation import (
    compressed_responsivity_a_per_w,
    compressed_responsivity_rational_a_per_w,
    saturated_photocurrent_a,
    saturated_photocurrent_exponential_a,
    saturated_photocurrent_rational_a,
    saturated_photocurrent_soft_clip_a,
    saturated_photocurrent_tanh_a,
)


def test_saturated_photocurrent_matches_linear_current_at_low_power() -> None:
    current_a = saturated_photocurrent_a(
        optical_power_w=1e-6,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    assert current_a == pytest.approx(0.8e-6, rel=1e-3)


def test_compressed_responsivity_decreases_with_power() -> None:
    optical_power_w = np.array([1e-6, 1e-4, 1e-3])

    responsivity = compressed_responsivity_a_per_w(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    assert responsivity[0] > responsivity[1] > responsivity[2]


def test_default_saturation_model_is_tanh() -> None:
    optical_power_w = np.array([0.0, 1e-6, 1e-3, 3e-3])

    default_current = saturated_photocurrent_a(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )
    tanh_current = saturated_photocurrent_tanh_a(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    np.testing.assert_allclose(default_current, tanh_current)


def test_compressed_responsivity_uses_small_signal_value_at_zero_power() -> None:
    responsivity = compressed_responsivity_a_per_w(
        optical_power_w=np.array([0.0, 1e-6]),
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    assert responsivity[0] == pytest.approx(0.8)
    assert responsivity[1] < 0.8


def test_original_rational_model_is_still_available() -> None:
    optical_power_w = 1e-3

    responsivity = compressed_responsivity_rational_a_per_w(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )
    current = saturated_photocurrent_rational_a(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    assert responsivity == pytest.approx(0.4)
    assert current == pytest.approx(0.4e-3)


def test_high_power_incremental_responsivity_decreases() -> None:
    optical_power_w = np.array([1e-6, 2e-6, 1e-3, 2e-3])
    current_a = saturated_photocurrent_a(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    low_power_slope = (current_a[1] - current_a[0]) / (optical_power_w[1] - optical_power_w[0])
    high_power_slope = (current_a[3] - current_a[2]) / (optical_power_w[3] - optical_power_w[2])

    assert high_power_slope < low_power_slope


def test_saturated_photocurrent_remains_non_negative_for_arrays() -> None:
    current_a = saturated_photocurrent_a(
        optical_power_w=np.array([0.0, 1e-6, 1e-3]),
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    assert np.all(current_a >= 0)


def test_saturation_helpers_return_float_for_scalar_input() -> None:
    responsivity = compressed_responsivity_a_per_w(
        optical_power_w=1e-6,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )
    current = saturated_photocurrent_a(
        optical_power_w=1e-6,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    assert isinstance(responsivity, float)
    assert isinstance(current, float)


@pytest.mark.parametrize(
    "model",
    [
        saturated_photocurrent_a,
        saturated_photocurrent_soft_clip_a,
        saturated_photocurrent_tanh_a,
        saturated_photocurrent_exponential_a,
    ],
)
def test_asymptotic_saturation_models_match_linear_current_at_low_power(model: Callable[..., object]) -> None:
    current_a = model(
        optical_power_w=1e-7,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    assert current_a == pytest.approx(0.8e-7, rel=1e-3)


@pytest.mark.parametrize(
    "model",
    [
        saturated_photocurrent_a,
        saturated_photocurrent_soft_clip_a,
        saturated_photocurrent_tanh_a,
        saturated_photocurrent_exponential_a,
    ],
)
def test_asymptotic_saturation_models_approach_saturation_current(model: Callable[..., object]) -> None:
    saturation_current_a = 0.8 * 1e-3
    current_a = model(
        optical_power_w=100e-3,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
    )

    assert current_a == pytest.approx(saturation_current_a, rel=1e-2)


def test_soft_clip_order_controls_transition_sharpness() -> None:
    optical_power_w = 1e-3

    gradual_current = saturated_photocurrent_soft_clip_a(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
        compression_order=1.0,
    )
    sharp_current = saturated_photocurrent_soft_clip_a(
        optical_power_w=optical_power_w,
        responsivity_a_per_w=0.8,
        saturation_power_w=1e-3,
        compression_order=4.0,
    )

    assert sharp_current > gradual_current


def test_asymptotic_saturation_models_return_float_for_scalar_input() -> None:
    assert isinstance(saturated_photocurrent_soft_clip_a(1e-6, 0.8, 1e-3), float)
    assert isinstance(saturated_photocurrent_tanh_a(1e-6, 0.8, 1e-3), float)
    assert isinstance(saturated_photocurrent_exponential_a(1e-6, 0.8, 1e-3), float)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (
            lambda: saturated_photocurrent_a(-1e-6, 0.8, 1e-3),
            "optical_power_w must be non-negative",
        ),
        (
            lambda: saturated_photocurrent_a(1e-6, 0.0, 1e-3),
            "responsivity_a_per_w must be positive",
        ),
        (
            lambda: saturated_photocurrent_a(1e-6, 0.8, 0.0),
            "saturation_power_w must be positive",
        ),
        (
            lambda: saturated_photocurrent_a(1e-6, 0.8, 1e-3, compression_order=0.0),
            "compression_order must be positive",
        ),
        (
            lambda: compressed_responsivity_a_per_w(np.array([0.0, -1e-6]), 0.8, 1e-3),
            "optical_power_w must be non-negative",
        ),
        (
            lambda: saturated_photocurrent_rational_a(1e-6, 0.8, 1e-3, compression_order=0.0),
            "compression_order must be positive",
        ),
        (
            lambda: saturated_photocurrent_soft_clip_a(1e-6, 0.8, 1e-3, compression_order=0.0),
            "compression_order must be positive",
        ),
        (
            lambda: saturated_photocurrent_tanh_a(1e-6, 0.8, 0.0),
            "saturation_power_w must be positive",
        ),
        (
            lambda: saturated_photocurrent_exponential_a(np.array([-1e-6]), 0.8, 1e-3),
            "optical_power_w must be non-negative",
        ),
    ],
)
def test_saturation_helpers_reject_invalid_inputs(call: Callable[[], object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
