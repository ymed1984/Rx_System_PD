from collections.abc import Callable

import numpy as np
import pytest

from oma_ber.ber import (
    ber_for_threshold,
    ber_from_gaussian_levels,
    ber_from_q,
    optimum_threshold,
    q_from_levels,
    qfunc,
)


def test_qfunc_zero_is_one_half() -> None:
    assert qfunc(0.0) == pytest.approx(0.5)


def test_qfunc_accepts_numpy_array() -> None:
    values = qfunc(np.array([0.0, 1.0]))

    assert values[0] == pytest.approx(0.5)
    assert values[1] < values[0]


def test_ber_for_threshold_uses_correct_error_directions() -> None:
    mu0_a = 0.0
    mu1_a = 2.0
    sigma0_a = 0.5
    sigma1_a = 0.5
    threshold_a = 1.0

    expected = 0.5 * (
        qfunc((threshold_a - mu0_a) / sigma0_a)
        + qfunc((mu1_a - threshold_a) / sigma1_a)
    )

    assert ber_for_threshold(mu0_a, mu1_a, sigma0_a, sigma1_a, threshold_a) == pytest.approx(expected)


def test_ber_decreases_when_current_separation_increases() -> None:
    close_ber = ber_for_threshold(
        mu0_a=0.0,
        mu1_a=1.0,
        sigma0_a=0.2,
        sigma1_a=0.2,
        threshold_a=0.5,
    )
    far_ber = ber_for_threshold(
        mu0_a=0.0,
        mu1_a=2.0,
        sigma0_a=0.2,
        sigma1_a=0.2,
        threshold_a=1.0,
    )

    assert far_ber < close_ber


def test_ber_increases_when_noise_increases() -> None:
    low_noise_ber = ber_for_threshold(
        mu0_a=0.0,
        mu1_a=1.0,
        sigma0_a=0.1,
        sigma1_a=0.1,
        threshold_a=0.5,
    )
    high_noise_ber = ber_for_threshold(
        mu0_a=0.0,
        mu1_a=1.0,
        sigma0_a=0.3,
        sigma1_a=0.3,
        threshold_a=0.5,
    )

    assert high_noise_ber > low_noise_ber


def test_optimum_threshold_lies_between_levels() -> None:
    threshold_a = optimum_threshold(mu0_a=0.0, mu1_a=1.0, sigma0_a=0.1, sigma1_a=0.3)

    assert 0.0 < threshold_a < 1.0


def test_optimum_threshold_handles_large_unequal_noise() -> None:
    threshold_a = optimum_threshold(mu0_a=0.0, mu1_a=1.0, sigma0_a=5.0, sigma1_a=1.0)

    assert 0.0 <= threshold_a <= 1.0


def test_optimum_threshold_is_accurate_at_microamp_scale() -> None:
    threshold_a = optimum_threshold(
        mu0_a=0.0,
        mu1_a=1e-6,
        sigma0_a=1e-7,
        sigma1_a=2e-7,
    )

    assert threshold_a == pytest.approx(3.470550625549095e-7, rel=1e-12)


@pytest.mark.parametrize("current_scale", [1e-6, 1e-3, 1.0, 1e3])
def test_optimum_threshold_is_invariant_to_current_scale(current_scale: float) -> None:
    mu0 = 1.0 * current_scale
    mu1 = 1.1 * current_scale
    sigma0 = 0.01 * current_scale
    sigma1 = 0.03 * current_scale

    threshold = optimum_threshold(mu0, mu1, sigma0, sigma1)
    normalized_threshold = (threshold - mu0) / (mu1 - mu0)

    assert normalized_threshold == pytest.approx(0.281624859661866, rel=1e-12)


def test_optimum_threshold_matches_dimensionless_reference_search() -> None:
    mu0_a = 2.5e-9
    mu1_a = 9.5e-9
    sigma0_a = 0.4e-9
    sigma1_a = 1.7e-9
    normalized_grid = np.linspace(0.0, 1.0, 1_000_001)
    threshold_grid_a = mu0_a + normalized_grid * (mu1_a - mu0_a)
    reference_bers = 0.5 * (
        qfunc((threshold_grid_a - mu0_a) / sigma0_a)
        + qfunc((mu1_a - threshold_grid_a) / sigma1_a)
    )
    reference_threshold_a = threshold_grid_a[np.argmin(reference_bers)]

    threshold_a = optimum_threshold(mu0_a, mu1_a, sigma0_a, sigma1_a)

    grid_step_a = threshold_grid_a[1] - threshold_grid_a[0]
    assert threshold_a == pytest.approx(reference_threshold_a, abs=grid_step_a)


def test_optimum_threshold_ber_is_no_worse_than_midpoint_for_extreme_noise_ratio() -> None:
    mu0_a = 0.0
    mu1_a = 1e-6
    sigma0_a = 5e-6
    sigma1_a = 0.1e-6
    threshold_a = optimum_threshold(mu0_a, mu1_a, sigma0_a, sigma1_a)
    midpoint_a = (mu0_a + mu1_a) / 2

    assert ber_for_threshold(mu0_a, mu1_a, sigma0_a, sigma1_a, threshold_a) <= ber_for_threshold(
        mu0_a,
        mu1_a,
        sigma0_a,
        sigma1_a,
        midpoint_a,
    )


def test_equal_noise_optimum_threshold_is_midpoint() -> None:
    assert optimum_threshold(
        mu0_a=1e-6,
        mu1_a=5e-6,
        sigma0_a=0.2e-6,
        sigma1_a=0.2e-6,
    ) == pytest.approx(3e-6)


def test_q_from_levels_matches_execplan_equation() -> None:
    assert q_from_levels(
        mu0_a=1e-6,
        mu1_a=5e-6,
        sigma0_a=0.5e-6,
        sigma1_a=1.5e-6,
    ) == pytest.approx(2.0)


def test_ber_from_q_decreases_with_q() -> None:
    assert ber_from_q(5.0) < ber_from_q(3.0)


def test_ber_from_gaussian_levels_returns_expected_keys() -> None:
    result = ber_from_gaussian_levels(
        mu0_a=0.0,
        mu1_a=1.0,
        sigma0_a=0.1,
        sigma1_a=0.2,
    )

    assert result["threshold_a"] == pytest.approx(optimum_threshold(0.0, 1.0, 0.1, 0.2))
    assert result["q_rx"] == pytest.approx(q_from_levels(0.0, 1.0, 0.1, 0.2))
    assert result["ber"] == pytest.approx(
        ber_for_threshold(0.0, 1.0, 0.1, 0.2, result["threshold_a"]),
    )
    assert result["ber_from_q"] == pytest.approx(ber_from_q(result["q_rx"]))


def test_ber_from_gaussian_levels_can_use_midpoint_threshold() -> None:
    result = ber_from_gaussian_levels(
        mu0_a=0.0,
        mu1_a=1.0,
        sigma0_a=0.1,
        sigma1_a=0.3,
        optimize_threshold=False,
    )

    assert result["threshold_a"] == pytest.approx(0.5)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: ber_for_threshold(1.0, 1.0, 0.1, 0.1, 1.0), "mu1_a must be larger"),
        (lambda: ber_for_threshold(0.0, 1.0, 0.0, 0.1, 0.5), "sigma0_a must be positive"),
        (lambda: ber_for_threshold(0.0, 1.0, 0.1, 0.0, 0.5), "sigma1_a must be positive"),
        (lambda: optimum_threshold(1.0, 0.0, 0.1, 0.1), "mu1_a must be larger"),
        (lambda: q_from_levels(0.0, 1.0, -0.1, 0.1), "sigma0_a must be positive"),
        (lambda: ber_from_q(-1.0), "q_rx must be non-negative"),
    ],
)
def test_ber_functions_reject_invalid_physical_inputs(call: Callable[[], float], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
