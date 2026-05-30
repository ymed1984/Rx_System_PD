"""Frequency-response import and impulse-response conversion helpers."""

from pathlib import Path

import numpy as np


def frequency_response_from_csv(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a complex frequency response from CSV.

    The CSV file must contain columns named frequency_hz, magnitude_db, and
    phase_deg. frequency_hz is in Hz, magnitude_db is voltage/current amplitude
    gain in dB, and phase_deg is phase in degrees. The returned response is a
    complex linear amplitude response.
    """
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.dtype.names is None:
        msg = "CSV must contain a header row."
        raise ValueError(msg)

    required_columns = {"frequency_hz", "magnitude_db", "phase_deg"}
    missing_columns = required_columns - set(data.dtype.names)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        msg = f"CSV is missing required columns: {missing}."
        raise ValueError(msg)

    frequency_hz = np.atleast_1d(np.asarray(data["frequency_hz"], dtype=float))
    magnitude_db = np.atleast_1d(np.asarray(data["magnitude_db"], dtype=float))
    phase_deg = np.atleast_1d(np.asarray(data["phase_deg"], dtype=float))

    _validate_frequency_response(frequency_hz, magnitude_db)
    response_complex = 10 ** (magnitude_db / 20) * np.exp(1j * np.deg2rad(phase_deg))
    return frequency_hz, response_complex


def impulse_response_from_frequency_response(
    frequency_hz: np.ndarray,
    response_complex: np.ndarray,
    sample_rate_hz: float,
    num_taps: int,
) -> np.ndarray:
    """Convert a complex frequency response into a real impulse response.

    frequency_hz is in Hz and response_complex is linear amplitude gain.
    sample_rate_hz is in samples/s. The response is linearly interpolated onto
    the real-FFT frequency grid. DC uses the first supplied response value, and
    frequencies above the supplied maximum hold the last supplied response
    value up to Nyquist. This helper is intended for deterministic waveform
    filtering, not standards-compliance measurement.
    """
    frequency_values = np.asarray(frequency_hz, dtype=float)
    response_values = np.asarray(response_complex, dtype=complex)
    _validate_frequency_response(frequency_values, response_values)

    if sample_rate_hz <= 0:
        msg = "sample_rate_hz must be positive."
        raise ValueError(msg)
    if num_taps <= 0:
        msg = "num_taps must be positive."
        raise ValueError(msg)
    if sample_rate_hz <= 2 * float(np.max(frequency_values)):
        msg = "sample_rate_hz must exceed twice the maximum frequency_hz."
        raise ValueError(msg)

    fft_frequency_hz = np.fft.rfftfreq(num_taps, d=1 / sample_rate_hz)
    real_interp = np.interp(
        fft_frequency_hz,
        frequency_values,
        response_values.real,
        left=response_values.real[0],
        right=response_values.real[-1],
    )
    imag_interp = np.interp(
        fft_frequency_hz,
        frequency_values,
        response_values.imag,
        left=response_values.imag[0],
        right=response_values.imag[-1],
    )
    interpolated_response = real_interp + 1j * imag_interp
    return np.fft.irfft(interpolated_response, n=num_taps)


def _validate_frequency_response(frequency_hz: np.ndarray, response_values: np.ndarray) -> None:
    if frequency_hz.ndim != 1:
        msg = "frequency_hz must be a one-dimensional array."
        raise ValueError(msg)
    if response_values.ndim != 1:
        msg = "response_complex must be a one-dimensional array."
        raise ValueError(msg)
    if frequency_hz.size == 0:
        msg = "frequency_hz must contain at least one frequency."
        raise ValueError(msg)
    if frequency_hz.size != response_values.size:
        msg = "response length must match frequency_hz length."
        raise ValueError(msg)
    if np.any(frequency_hz <= 0):
        msg = "frequency_hz values must be positive."
        raise ValueError(msg)
    if not np.all(np.diff(frequency_hz) > 0):
        msg = "frequency_hz values must be strictly increasing."
        raise ValueError(msg)
