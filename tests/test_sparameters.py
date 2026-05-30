from collections.abc import Callable

import numpy as np
import pytest

from oma_ber.sparameters import frequency_response_from_csv, impulse_response_from_frequency_response


def test_frequency_response_from_csv_reads_complex_response(tmp_path) -> None:
    csv_path = tmp_path / "response.csv"
    csv_path.write_text(
        "frequency_hz,magnitude_db,phase_deg\n"
        "1000000000,0,0\n"
        "2000000000,-6,90\n",
        encoding="utf-8",
    )

    frequency_hz, response_complex = frequency_response_from_csv(csv_path)

    np.testing.assert_allclose(frequency_hz, np.array([1e9, 2e9]))
    assert response_complex[0] == pytest.approx(1.0 + 0.0j)
    assert response_complex[1].real == pytest.approx(0.0, abs=1e-12)
    assert response_complex[1].imag == pytest.approx(10 ** (-6 / 20))


def test_flat_response_produces_near_delta_impulse() -> None:
    frequency_hz = np.array([1e9, 2e9, 3e9])
    response_complex = np.ones_like(frequency_hz, dtype=complex)

    impulse_response = impulse_response_from_frequency_response(
        frequency_hz,
        response_complex,
        sample_rate_hz=10e9,
        num_taps=16,
    )

    assert impulse_response[0] == pytest.approx(1.0)
    np.testing.assert_allclose(impulse_response[1:], 0.0, atol=1e-12)


def test_first_order_synthetic_response_is_low_pass() -> None:
    frequency_hz = np.linspace(0.5e9, 20e9, 64)
    bandwidth_3db_hz = 5e9
    response_complex = 1 / (1 + 1j * frequency_hz / bandwidth_3db_hz)

    impulse_response = impulse_response_from_frequency_response(
        frequency_hz,
        response_complex,
        sample_rate_hz=50e9,
        num_taps=128,
    )
    recovered_response = np.fft.rfft(impulse_response)
    recovered_frequency_hz = np.fft.rfftfreq(len(impulse_response), d=1 / 50e9)
    low_bin = np.argmin(np.abs(recovered_frequency_hz - 1e9))
    high_bin = np.argmin(np.abs(recovered_frequency_hz - 15e9))

    assert abs(recovered_response[low_bin]) > abs(recovered_response[high_bin])


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (
            lambda: impulse_response_from_frequency_response(np.array([[1e9]]), np.array([1.0 + 0j]), 10e9, 16),
            "frequency_hz must be a one-dimensional array",
        ),
        (
            lambda: impulse_response_from_frequency_response(np.array([1e9]), np.array([[1.0 + 0j]]), 10e9, 16),
            "response_complex must be a one-dimensional array",
        ),
        (
            lambda: impulse_response_from_frequency_response(np.array([]), np.array([], dtype=complex), 10e9, 16),
            "frequency_hz must contain at least one frequency",
        ),
        (
            lambda: impulse_response_from_frequency_response(np.array([1e9, 2e9]), np.array([1.0 + 0j]), 10e9, 16),
            "response length must match frequency_hz length",
        ),
        (
            lambda: impulse_response_from_frequency_response(np.array([0.0]), np.array([1.0 + 0j]), 10e9, 16),
            "frequency_hz values must be positive",
        ),
        (
            lambda: impulse_response_from_frequency_response(np.array([2e9, 1e9]), np.ones(2), 10e9, 16),
            "frequency_hz values must be strictly increasing",
        ),
        (
            lambda: impulse_response_from_frequency_response(np.array([1e9]), np.array([1.0 + 0j]), 0.0, 16),
            "sample_rate_hz must be positive",
        ),
        (
            lambda: impulse_response_from_frequency_response(np.array([1e9]), np.array([1.0 + 0j]), 10e9, 0),
            "num_taps must be positive",
        ),
        (
            lambda: impulse_response_from_frequency_response(np.array([5e9]), np.array([1.0 + 0j]), 10e9, 16),
            "sample_rate_hz must exceed twice the maximum frequency_hz",
        ),
    ],
)
def test_impulse_response_rejects_invalid_inputs(call: Callable[[], object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        call()


def test_frequency_response_from_csv_rejects_missing_columns(tmp_path) -> None:
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text("frequency_hz,magnitude_db\n1000000000,0\n", encoding="utf-8")

    with pytest.raises(ValueError, match="CSV is missing required columns: phase_deg"):
        frequency_response_from_csv(csv_path)
