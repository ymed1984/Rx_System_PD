"""Physical-unit waveform containers and photodetection helpers."""

from dataclasses import dataclass, field

import numpy as np

from oma_ber.modulation import nrz_levels_from_powers
from oma_ber.photodiode import Photodiode
from oma_ber.saturation import saturated_photocurrent_tanh_a
from oma_ber.time_domain._validation import readonly_float_array, validated_bits
from oma_ber.time_domain.timebase import TimeGrid


@dataclass(frozen=True)
class OpticalPowerWaveform:
    """Sampled optical power waveform in watts (W)."""

    power_w: np.ndarray = field(repr=False)
    time_grid: TimeGrid

    def __post_init__(self) -> None:
        values = readonly_float_array(self.power_w, "power_w")
        if np.any(values < 0):
            msg = "power_w must be non-negative."
            raise ValueError(msg)
        object.__setattr__(self, "power_w", values)


@dataclass(frozen=True)
class PhotocurrentWaveform:
    """Sampled photodiode current waveforms in amperes (A).

    signal_current_a excludes the common dark-current mean. total_current_a
    includes it. Both arrays are before the photodiode electrical response.
    """

    signal_current_a: np.ndarray = field(repr=False)
    total_current_a: np.ndarray = field(repr=False)
    time_grid: TimeGrid

    def __post_init__(self) -> None:
        signal = readonly_float_array(self.signal_current_a, "signal_current_a")
        total = readonly_float_array(self.total_current_a, "total_current_a")
        if signal.size != total.size:
            msg = "signal_current_a and total_current_a must have the same length."
            raise ValueError(msg)
        if np.any(signal < 0) or np.any(total < 0):
            msg = "photocurrent waveforms must be non-negative."
            raise ValueError(msg)
        object.__setattr__(self, "signal_current_a", signal)
        object.__setattr__(self, "total_current_a", total)


@dataclass(frozen=True)
class VoltageWaveform:
    """Sampled receiver voltage waveform in volts (V)."""

    voltage_v: np.ndarray = field(repr=False)
    time_grid: TimeGrid

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "voltage_v",
            readonly_float_array(self.voltage_v, "voltage_v"),
        )


def nrz_optical_power_waveform(
    bits: np.ndarray,
    p0_w: float,
    p1_w: float,
    time_grid: TimeGrid,
) -> OpticalPowerWaveform:
    """Create a rectangular NRZ/OOK optical-power waveform in watts (W)."""
    bit_values = validated_bits(bits)
    levels = nrz_levels_from_powers(p0_w=p0_w, p1_w=p1_w)
    symbol_levels_w = np.where(bit_values == 0, levels.p0_w, levels.p1_w)
    power_w = np.repeat(symbol_levels_w, time_grid.samples_per_symbol)
    return OpticalPowerWaveform(power_w=power_w, time_grid=time_grid)


def photodetect_power_waveform(
    optical_waveform: OpticalPowerWaveform,
    pd: Photodiode,
) -> PhotocurrentWaveform:
    """Convert optical power in W to pre-bandwidth PD currents in A.

    Linear operation uses i_ph(t) = responsivity_a_per_w * power_w(t). If the
    Photodiode has saturation_power_w, its documented static tanh compression
    model is applied sample by sample. The common dark-current mean is added
    only to total_current_a.
    """
    power_w = optical_waveform.power_w
    if pd.saturation_power_w is None:
        signal_current_a = pd.responsivity_a_per_w * power_w
    else:
        signal_current_a = saturated_photocurrent_tanh_a(
            optical_power_w=power_w,
            responsivity_a_per_w=pd.responsivity_a_per_w,
            saturation_power_w=pd.saturation_power_w,
        )
    total_current_a = np.asarray(signal_current_a) + pd.dark_current_a
    return PhotocurrentWaveform(
        signal_current_a=np.asarray(signal_current_a),
        total_current_a=total_current_a,
        time_grid=optical_waveform.time_grid,
    )
