"""Deterministic PD/TIA waveform simulation at the receiver boundary."""

from dataclasses import dataclass, field
from math import ceil, isclose

import numpy as np

from oma_ber.modulation import nrz_levels_from_powers
from oma_ber.photodiode import Photodiode
from oma_ber.time_domain._validation import readonly_float_array, validated_bits
from oma_ber.time_domain.timebase import TimeGrid
from oma_ber.time_domain.transfer import (
    DiscreteTransferFunction,
    apply_transfer,
    identity_transfer,
)
from oma_ber.time_domain.waveforms import (
    nrz_optical_power_waveform,
    photodetect_power_waveform,
)


@dataclass(frozen=True)
class ReceiverWaveformResult:
    """Deterministic receiver waveforms with explicit physical units."""

    time_grid: TimeGrid
    bits: np.ndarray = field(repr=False)
    p0_w: float
    p1_w: float
    warmup_symbols: int
    photodiode: Photodiode
    pd_current_response: DiscreteTransferFunction
    tia_transimpedance_response: DiscreteTransferFunction | None
    optical_power_w: np.ndarray = field(repr=False)
    raw_signal_current_a: np.ndarray = field(repr=False)
    raw_total_current_a: np.ndarray = field(repr=False)
    pd_output_current_a: np.ndarray = field(repr=False)
    tia_output_voltage_v: np.ndarray | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        bits = validated_bits(self.bits)
        nrz_levels_from_powers(p0_w=self.p0_w, p1_w=self.p1_w)
        if not isinstance(self.photodiode, Photodiode):
            msg = "photodiode must be a Photodiode."
            raise ValueError(msg)
        if self.pd_current_response.response_kind != "dimensionless":
            msg = "pd_current_response must have response_kind='dimensionless'."
            raise ValueError(msg)
        _validate_matching_sample_rate(
            self.pd_current_response,
            self.time_grid.sample_rate_hz,
            "pd_current_response",
        )
        if self.tia_transimpedance_response is not None:
            if self.tia_transimpedance_response.response_kind != "transimpedance_ohm":
                msg = (
                    "tia_transimpedance_response must have "
                    "response_kind='transimpedance_ohm'."
                )
                raise ValueError(msg)
            _validate_matching_sample_rate(
                self.tia_transimpedance_response,
                self.time_grid.sample_rate_hz,
                "tia_transimpedance_response",
            )
        if isinstance(self.warmup_symbols, bool) or not isinstance(
            self.warmup_symbols, int
        ):
            msg = "warmup_symbols must be an integer."
            raise ValueError(msg)
        if self.warmup_symbols < 0:
            msg = "warmup_symbols must be non-negative."
            raise ValueError(msg)
        arrays = {
            "optical_power_w": self.optical_power_w,
            "raw_signal_current_a": self.raw_signal_current_a,
            "raw_total_current_a": self.raw_total_current_a,
            "pd_output_current_a": self.pd_output_current_a,
        }
        expected_samples = bits.size * self.time_grid.samples_per_symbol
        for name, values in arrays.items():
            array = readonly_float_array(values, name)
            if array.size != expected_samples:
                msg = f"{name} length must equal bits * samples_per_symbol."
                raise ValueError(msg)
            object.__setattr__(self, name, array)

        if self.tia_output_voltage_v is not None:
            voltage = readonly_float_array(
                self.tia_output_voltage_v, "tia_output_voltage_v"
            )
            if voltage.size != expected_samples:
                msg = (
                    "tia_output_voltage_v length must equal bits * samples_per_symbol."
                )
                raise ValueError(msg)
            object.__setattr__(self, "tia_output_voltage_v", voltage)
        object.__setattr__(self, "bits", bits)


def simulate_pd_tia_waveform(
    bits: np.ndarray,
    p0_w: float,
    p1_w: float,
    time_grid: TimeGrid,
    pd: Photodiode,
    pd_current_response: DiscreteTransferFunction | None = None,
    tia_transimpedance_response: DiscreteTransferFunction | None = None,
    warmup_symbols: int | None = None,
    settling_tolerance: float = 1e-9,
) -> ReceiverWaveformResult:
    """Simulate deterministic PD-current and optional TIA-voltage waveforms.

    p0_w and p1_w are optical powers at the PD input reference plane. The PD
    static conversion and saturation are applied before the dimensionless PD
    electrical response. The optional TIA response must have transimpedance
    gain in ohms and produces a voltage waveform in V.

    A periodic continuation of bits is prepended as warm-up. If warmup_symbols
    is omitted, it is estimated from the supplied causal filter poles and the
    requested settling_tolerance.
    """
    bit_values = validated_bits(bits)
    nrz_levels_from_powers(p0_w=p0_w, p1_w=p1_w)
    sample_rate_hz = time_grid.sample_rate_hz

    pd_response = pd_current_response or identity_transfer(sample_rate_hz)
    if pd_response.response_kind != "dimensionless":
        msg = "pd_current_response must have response_kind='dimensionless'."
        raise ValueError(msg)
    _validate_matching_sample_rate(pd_response, sample_rate_hz, "pd_current_response")

    if tia_transimpedance_response is not None:
        if tia_transimpedance_response.response_kind != "transimpedance_ohm":
            msg = (
                "tia_transimpedance_response must have "
                "response_kind='transimpedance_ohm'."
            )
            raise ValueError(msg)
        _validate_matching_sample_rate(
            tia_transimpedance_response,
            sample_rate_hz,
            "tia_transimpedance_response",
        )

    if warmup_symbols is None:
        settling_samples = pd_response.settling_samples(settling_tolerance)
        if tia_transimpedance_response is not None:
            # Cascaded filters can settle more slowly than either stage alone.
            # Summing their individual conservative estimates avoids retaining
            # a zero-state transient in the evaluated symbols.
            settling_samples += tia_transimpedance_response.settling_samples(
                settling_tolerance,
            )
        resolved_warmup_symbols = ceil(settling_samples / time_grid.samples_per_symbol)
    else:
        if isinstance(warmup_symbols, bool) or not isinstance(warmup_symbols, int):
            msg = "warmup_symbols must be an integer when provided."
            raise ValueError(msg)
        if warmup_symbols < 0:
            msg = "warmup_symbols must be non-negative."
            raise ValueError(msg)
        resolved_warmup_symbols = warmup_symbols

    warmup_bits = _periodic_prefix(bit_values, resolved_warmup_symbols)
    full_bits = np.concatenate((warmup_bits, bit_values))
    optical = nrz_optical_power_waveform(full_bits, p0_w, p1_w, time_grid)
    photocurrent = photodetect_power_waveform(optical, pd)
    pd_output_current_a = apply_transfer(photocurrent.total_current_a, pd_response)
    tia_output_voltage_v = (
        apply_transfer(pd_output_current_a, tia_transimpedance_response)
        if tia_transimpedance_response is not None
        else None
    )

    start = resolved_warmup_symbols * time_grid.samples_per_symbol
    return ReceiverWaveformResult(
        time_grid=time_grid,
        bits=bit_values,
        p0_w=p0_w,
        p1_w=p1_w,
        warmup_symbols=resolved_warmup_symbols,
        photodiode=pd,
        pd_current_response=pd_response,
        tia_transimpedance_response=tia_transimpedance_response,
        optical_power_w=optical.power_w[start:],
        raw_signal_current_a=photocurrent.signal_current_a[start:],
        raw_total_current_a=photocurrent.total_current_a[start:],
        pd_output_current_a=pd_output_current_a[start:],
        tia_output_voltage_v=(
            tia_output_voltage_v[start:] if tia_output_voltage_v is not None else None
        ),
    )


def _periodic_prefix(bits: np.ndarray, num_symbols: int) -> np.ndarray:
    if num_symbols == 0:
        return np.empty(0, dtype=np.int_)
    indices = np.arange(-num_symbols, 0) % bits.size
    return np.asarray(bits[indices], dtype=np.int_)


def _validate_matching_sample_rate(
    response: DiscreteTransferFunction,
    sample_rate_hz: float,
    name: str,
) -> None:
    if not isclose(response.sample_rate_hz, sample_rate_hz, rel_tol=1e-12, abs_tol=0.0):
        msg = f"{name}.sample_rate_hz must match time_grid.sample_rate_hz."
        raise ValueError(msg)
