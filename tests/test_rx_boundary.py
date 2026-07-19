from collections.abc import Callable

import numpy as np
import pytest

from oma_ber.modulation import nrz_levels_from_powers
from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.sweep import calculate_ber_from_oma, calculate_ber_from_optical_levels


def example_pd() -> Photodiode:
    return Photodiode(responsivity_a_per_w=0.8, dark_current_a=1e-9)


def example_rx() -> Receiver:
    return Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )


def test_levels_from_powers_recovers_oma_er_and_average() -> None:
    levels = nrz_levels_from_powers(p0_w=1e-4, p1_w=4e-4)

    assert levels.oma_w == pytest.approx(3e-4)
    assert levels.pavg_w == pytest.approx(2.5e-4)
    assert levels.er_linear == pytest.approx(4.0)
    assert levels.er_db == pytest.approx(10 * np.log10(4.0))


def test_levels_from_powers_supports_ideal_extinction() -> None:
    levels = nrz_levels_from_powers(p0_w=0.0, p1_w=1e-4)

    assert levels.oma_w == pytest.approx(1e-4)
    assert np.isinf(levels.er_linear)
    assert np.isinf(levels.er_db)


def test_optical_level_receiver_api_matches_oma_wrapper() -> None:
    from_oma = calculate_ber_from_oma(-10.0, er_db=6.0, pd=example_pd(), rx=example_rx())
    from_levels = calculate_ber_from_optical_levels(
        p0_w=from_oma["p0_w"],
        p1_w=from_oma["p1_w"],
        pd=example_pd(),
        rx=example_rx(),
    )

    for key in ("oma_w", "er_db", "i0_a", "i1_a", "sigma0_a", "sigma1_a", "ber"):
        assert from_levels[key] == pytest.approx(from_oma[key])


def test_optical_level_receiver_api_supports_zero_level_power() -> None:
    result = calculate_ber_from_optical_levels(
        p0_w=0.0,
        p1_w=1e-4,
        pd=example_pd(),
        rx=example_rx(),
    )

    assert result["i0_a"] == 0.0
    assert result["responsivity0_effective_a_per_w"] == pytest.approx(
        example_pd().responsivity_a_per_w,
    )
    assert np.isinf(result["er_db"])


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: nrz_levels_from_powers(-1e-6, 1e-6), "p0_w must be finite and non-negative"),
        (
            lambda: nrz_levels_from_powers(float("nan"), 1e-6),
            "p0_w must be finite and non-negative",
        ),
        (
            lambda: nrz_levels_from_powers(1e-6, 1e-6),
            "p1_w must be finite and larger than p0_w",
        ),
        (
            lambda: nrz_levels_from_powers(2e-6, 1e-6),
            "p1_w must be finite and larger than p0_w",
        ),
    ],
)
def test_levels_from_powers_rejects_invalid_levels(
    call: Callable[[], object],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        call()
