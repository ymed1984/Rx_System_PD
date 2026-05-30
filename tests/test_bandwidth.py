from collections.abc import Callable

import numpy as np
import pytest

from oma_ber.bandwidth import (
    apply_oma_penalty_db,
    bandwidth_penalty_db,
    first_order_lowpass_mag,
)


def test_first_order_lowpass_mag_is_one_at_dc() -> None:
    assert first_order_lowpass_mag(0.0, bandwidth_3db_hz=25e9) == pytest.approx(1.0)


def test_first_order_lowpass_mag_is_minus_3db_at_bandwidth() -> None:
    assert first_order_lowpass_mag(25e9, bandwidth_3db_hz=25e9) == pytest.approx(1 / np.sqrt(2))


def test_first_order_lowpass_mag_accepts_array() -> None:
    magnitudes = first_order_lowpass_mag(np.array([0.0, 25e9]), bandwidth_3db_hz=25e9)

    assert magnitudes[0] == pytest.approx(1.0)
    assert magnitudes[1] == pytest.approx(1 / np.sqrt(2))


def test_bandwidth_penalty_db_is_zero_at_dc() -> None:
    assert bandwidth_penalty_db(0.0, bandwidth_3db_hz=25e9) == pytest.approx(0.0)


def test_bandwidth_penalty_db_is_about_3db_at_bandwidth() -> None:
    assert bandwidth_penalty_db(25e9, bandwidth_3db_hz=25e9) == pytest.approx(3.0103, rel=1e-4)


def test_apply_oma_penalty_db_subtracts_penalty() -> None:
    assert apply_oma_penalty_db(oma_dbm=-10.0, penalty_db=3.0) == pytest.approx(-13.0)


def test_bandwidth_penalty_worsens_oma() -> None:
    penalty_db = bandwidth_penalty_db(signal_frequency_hz=20e9, bandwidth_3db_hz=10e9)

    assert apply_oma_penalty_db(-10.0, penalty_db) < -10.0


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: first_order_lowpass_mag(1e9, bandwidth_3db_hz=0.0), "bandwidth_3db_hz must be positive"),
        (lambda: first_order_lowpass_mag(-1.0, bandwidth_3db_hz=25e9), "frequency_hz must be non-negative"),
        (lambda: bandwidth_penalty_db(-1.0, bandwidth_3db_hz=25e9), "frequency_hz must be non-negative"),
        (lambda: apply_oma_penalty_db(-10.0, penalty_db=-0.1), "penalty_db must be non-negative"),
    ],
)
def test_bandwidth_helpers_reject_invalid_inputs(call: Callable[[], float], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
