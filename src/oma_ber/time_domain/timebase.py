"""Exact time-grid definitions for symbol-rate waveform analysis."""

from dataclasses import dataclass
from math import isfinite

import numpy as np


@dataclass(frozen=True)
class TimeGrid:
    """Uniform waveform time grid with an exact integer samples/UI ratio.

    symbol_rate_baud is in symbols/s and samples_per_symbol is the exact
    integer number of waveform samples per unit interval (UI).  sample_rate_hz
    and sample_interval_s are derived, so no silent rate rounding occurs.
    """

    symbol_rate_baud: float
    samples_per_symbol: int

    def __post_init__(self) -> None:
        if not isfinite(self.symbol_rate_baud) or self.symbol_rate_baud <= 0:
            msg = "symbol_rate_baud must be finite and positive."
            raise ValueError(msg)
        if isinstance(self.samples_per_symbol, bool) or not isinstance(
            self.samples_per_symbol,
            int,
        ):
            msg = "samples_per_symbol must be an integer."
            raise ValueError(msg)
        if self.samples_per_symbol < 2:
            msg = "samples_per_symbol must be at least 2."
            raise ValueError(msg)

    @property
    def unit_interval_s(self) -> float:
        """Return one unit interval in seconds (s)."""
        return 1.0 / self.symbol_rate_baud

    @property
    def sample_rate_hz(self) -> float:
        """Return the exact waveform sample rate in samples/s."""
        return self.symbol_rate_baud * self.samples_per_symbol

    @property
    def sample_interval_s(self) -> float:
        """Return the waveform sample interval in seconds (s)."""
        return 1.0 / self.sample_rate_hz

    def sample_times_s(self, num_samples: int) -> np.ndarray:
        """Return sample times in seconds for num_samples waveform samples."""
        if isinstance(num_samples, bool) or not isinstance(num_samples, int):
            msg = "num_samples must be an integer."
            raise ValueError(msg)
        if num_samples <= 0:
            msg = "num_samples must be positive."
            raise ValueError(msg)
        return np.arange(num_samples, dtype=float) * self.sample_interval_s
