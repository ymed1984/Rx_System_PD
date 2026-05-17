"""Receiver parameter containers."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Receiver:
    """Receiver noise parameters with SI units.

    noise_bandwidth_hz is in Hz, input_current_noise_density_a_per_sqrt_hz
    is in A/sqrt(Hz), and rin_db_per_hz is in dB/Hz when provided.
    """

    noise_bandwidth_hz: float
    input_current_noise_density_a_per_sqrt_hz: float
    rin_db_per_hz: float | None = None

    def __post_init__(self) -> None:
        if self.noise_bandwidth_hz <= 0:
            msg = "noise_bandwidth_hz must be positive."
            raise ValueError(msg)
        if self.input_current_noise_density_a_per_sqrt_hz < 0:
            msg = "input_current_noise_density_a_per_sqrt_hz must be non-negative."
            raise ValueError(msg)
