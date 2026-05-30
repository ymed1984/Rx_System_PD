from collections.abc import Callable

import numpy as np
import pytest

from oma_ber.waveform import nrz_bits_to_levels, prbs_bits, samples_per_symbol


def test_samples_per_symbol_from_rates() -> None:
    assert samples_per_symbol(sample_rate_hz=100e9, symbol_rate_baud=25e9) == 4


def test_nrz_bits_to_levels_repeats_known_sequence() -> None:
    bits = np.array([0, 1, 1, 0])

    waveform = nrz_bits_to_levels(bits, low_level=0.0, high_level=2.0, samples_per_symbol=3)

    np.testing.assert_allclose(
        waveform,
        np.array([0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0]),
    )


def test_nrz_bits_to_levels_output_length() -> None:
    bits = np.array([0, 1, 0, 1, 1])

    waveform = nrz_bits_to_levels(bits, low_level=1e-6, high_level=3e-6, samples_per_symbol=4)

    assert len(waveform) == len(bits) * 4


def test_prbs_bits_is_deterministic_for_fixed_seed() -> None:
    bits_a = prbs_bits(num_bits=32, seed=123)
    bits_b = prbs_bits(num_bits=32, seed=123)

    np.testing.assert_array_equal(bits_a, bits_b)
    assert len(bits_a) == 32
    assert set(bits_a.tolist()) <= {0, 1}


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: samples_per_symbol(sample_rate_hz=0.0, symbol_rate_baud=25e9), "sample_rate_hz must be positive"),
        (lambda: samples_per_symbol(sample_rate_hz=100e9, symbol_rate_baud=0.0), "symbol_rate_baud must be positive"),
        (
            lambda: samples_per_symbol(sample_rate_hz=25e9, symbol_rate_baud=25e9),
            "must be at least 2 samples/symbol",
        ),
        (
            lambda: nrz_bits_to_levels(np.array([0, 1]), 1.0, 1.0, 2),
            "high_level must be greater than low_level",
        ),
        (
            lambda: nrz_bits_to_levels(np.array([0, 1]), 0.0, 1.0, 1),
            "samples_per_symbol must be at least 2",
        ),
        (
            lambda: nrz_bits_to_levels(np.array([[0, 1]]), 0.0, 1.0, 2),
            "bits must be a one-dimensional array",
        ),
        (
            lambda: nrz_bits_to_levels(np.array([]), 0.0, 1.0, 2),
            "bits must contain at least one bit",
        ),
        (
            lambda: nrz_bits_to_levels(np.array([0, 2]), 0.0, 1.0, 2),
            "bits must contain only 0 and 1",
        ),
        (lambda: prbs_bits(num_bits=0), "num_bits must be positive"),
    ],
)
def test_waveform_helpers_reject_invalid_inputs(call: Callable[[], object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()
