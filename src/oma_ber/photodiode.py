"""Photodiode parameter containers."""

from dataclasses import dataclass


def _validate_optional_positive(value: float | None, name: str) -> None:
    if value is not None and value <= 0:
        msg = f"{name} must be positive when provided."
        raise ValueError(msg)


@dataclass(frozen=True)
class Photodiode:
    """Photodiode parameters with SI units.

    responsivity_a_per_w is in A/W, dark_current_a is in A,
    bandwidth_3db_hz is in Hz, capacitance_f is in F,
    saturation_power_w is in W, return_loss_db is in dB, bias_v is in V,
    shunt_resistance_ohm is in ohms, and temperature_k is in K.
    """

    responsivity_a_per_w: float
    dark_current_a: float
    bandwidth_3db_hz: float | None = None
    capacitance_f: float | None = None
    saturation_power_w: float | None = None
    return_loss_db: float | None = None
    bias_v: float | None = None
    shunt_resistance_ohm: float | None = None
    temperature_k: float = 300.0

    def __post_init__(self) -> None:
        if self.responsivity_a_per_w <= 0:
            msg = "responsivity_a_per_w must be positive."
            raise ValueError(msg)
        if self.dark_current_a < 0:
            msg = "dark_current_a must be non-negative."
            raise ValueError(msg)
        if self.temperature_k <= 0:
            msg = "temperature_k must be positive."
            raise ValueError(msg)

        _validate_optional_positive(self.bandwidth_3db_hz, "bandwidth_3db_hz")
        _validate_optional_positive(self.capacitance_f, "capacitance_f")
        _validate_optional_positive(self.saturation_power_w, "saturation_power_w")
        _validate_optional_positive(self.return_loss_db, "return_loss_db")
        _validate_optional_positive(self.bias_v, "bias_v")
        _validate_optional_positive(self.shunt_resistance_ohm, "shunt_resistance_ohm")
