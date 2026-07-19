"""PSD-based Gaussian noise propagation for the time-domain PD/TIA path."""

from dataclasses import dataclass, field
from math import isfinite

import numpy as np
from scipy.signal import lfilter

from oma_ber.noise import K_B, Q_E
from oma_ber.time_domain._validation import readonly_float_array
from oma_ber.time_domain.receiver import ReceiverWaveformResult
from oma_ber.time_domain.transfer import transfer_impulse_response
from oma_ber.units import rin_db_per_hz_to_linear


@dataclass(frozen=True)
class TimeDomainNoiseModel:
    """Flat input-noise spectra used by the Phase 2 Gaussian model.

    tia_input_current_noise_density_a_per_sqrt_hz is the one-sided TIA
    input-referred current-noise amplitude density in A/sqrt(Hz).
    rin_db_per_hz is the one-sided laser RIN power spectral density in dB/Hz.
    Photocurrent shot noise, dark-current shot noise, and optional PD shunt
    thermal noise use parameters already stored by Photodiode.
    """

    tia_input_current_noise_density_a_per_sqrt_hz: float
    rin_db_per_hz: float | None = None
    include_shunt_thermal_noise: bool = True

    def __post_init__(self) -> None:
        density = self.tia_input_current_noise_density_a_per_sqrt_hz
        if not isfinite(density) or density < 0:
            msg = (
                "tia_input_current_noise_density_a_per_sqrt_hz must be "
                "finite and non-negative."
            )
            raise ValueError(msg)
        if self.rin_db_per_hz is not None and not isfinite(self.rin_db_per_hz):
            msg = "rin_db_per_hz must be finite when provided."
            raise ValueError(msg)
        if not isinstance(self.include_shunt_thermal_noise, bool):
            msg = "include_shunt_thermal_noise must be a bool."
            raise ValueError(msg)


@dataclass(frozen=True)
class TiaOutputNoiseResult:
    """Pattern-dependent TIA-output Gaussian noise moments in volts."""

    photocurrent_shot_variance_v2: np.ndarray = field(repr=False)
    dark_current_shot_variance_v2: np.ndarray = field(repr=False)
    tia_variance_v2: np.ndarray = field(repr=False)
    rin_variance_v2: np.ndarray = field(repr=False)
    shunt_thermal_variance_v2: np.ndarray = field(repr=False)
    total_variance_v2: np.ndarray = field(repr=False)
    total_lag1_covariance_v2: np.ndarray = field(repr=False)
    total_sigma_v: np.ndarray = field(repr=False)
    pd_tia_noise_bandwidth_hz: float
    tia_noise_bandwidth_hz: float
    impulse_tolerance: float

    def __post_init__(self) -> None:
        component_names = (
            "photocurrent_shot_variance_v2",
            "dark_current_shot_variance_v2",
            "tia_variance_v2",
            "rin_variance_v2",
            "shunt_thermal_variance_v2",
            "total_variance_v2",
            "total_lag1_covariance_v2",
            "total_sigma_v",
        )
        expected_size: int | None = None
        for name in component_names:
            values = readonly_float_array(getattr(self, name), name)
            if expected_size is None:
                expected_size = values.size
            elif values.size != expected_size:
                msg = "all TIA-output noise arrays must have the same length."
                raise ValueError(msg)
            object.__setattr__(self, name, values)

        variance_names = component_names[:6]
        if any(np.any(getattr(self, name) < 0) for name in variance_names):
            msg = "noise variance arrays must be non-negative."
            raise ValueError(msg)
        if (
            not isfinite(self.pd_tia_noise_bandwidth_hz)
            or self.pd_tia_noise_bandwidth_hz <= 0
        ):
            msg = "pd_tia_noise_bandwidth_hz must be finite and positive."
            raise ValueError(msg)
        if (
            not isfinite(self.tia_noise_bandwidth_hz)
            or self.tia_noise_bandwidth_hz <= 0
        ):
            msg = "tia_noise_bandwidth_hz must be finite and positive."
            raise ValueError(msg)
        if not isfinite(self.impulse_tolerance) or not 0 < self.impulse_tolerance < 1:
            msg = "impulse_tolerance must be finite and between 0 and 1."
            raise ValueError(msg)


def calculate_tia_output_noise(
    waveforms: ReceiverWaveformResult,
    noise_model: TimeDomainNoiseModel,
    impulse_tolerance: float = 1e-12,
) -> TiaOutputNoiseResult:
    """Propagate one-sided PD/TIA noise PSDs to TIA-output moments.

    Photocurrent shot, dark-current shot, shunt thermal, and RIN noise originate
    before the PD electrical response and therefore pass through H_PD*Z_TIA.
    TIA input-current noise originates after H_PD and passes through Z_TIA
    only. For a one-sided white PSD S1 sampled at fs, the unfiltered sample
    variance is S1*fs/2. Pattern-dependent PSDs are periodically extended so
    the evaluated result has no artificial noise-startup transient.
    """
    if not isinstance(noise_model, TimeDomainNoiseModel):
        msg = "noise_model must be a TimeDomainNoiseModel."
        raise ValueError(msg)
    if not isfinite(impulse_tolerance) or not 0 < impulse_tolerance < 1:
        msg = "impulse_tolerance must be finite and between 0 and 1."
        raise ValueError(msg)
    tia_response = waveforms.tia_transimpedance_response
    if tia_response is None or waveforms.tia_output_voltage_v is None:
        msg = "TIA-output noise requires a tia_transimpedance_response."
        raise ValueError(msg)

    pd = waveforms.photodiode
    sample_rate_hz = waveforms.time_grid.sample_rate_hz
    pd_impulse = transfer_impulse_response(
        waveforms.pd_current_response,
        tolerance=impulse_tolerance,
    )
    tia_impulse = transfer_impulse_response(tia_response, tolerance=impulse_tolerance)
    pd_tia_impulse = np.convolve(pd_impulse, tia_impulse)

    photocurrent_shot_psd_a2_per_hz = 2 * Q_E * waveforms.raw_signal_current_a
    dark_current_shot_psd_a2_per_hz = np.full_like(
        waveforms.raw_signal_current_a,
        2 * Q_E * pd.dark_current_fano_factor * pd.dark_current_a,
    )
    tia_psd_a2_per_hz = np.full_like(
        waveforms.raw_signal_current_a,
        noise_model.tia_input_current_noise_density_a_per_sqrt_hz**2,
    )
    rin_psd_a2_per_hz = _rin_current_psd_a2_per_hz(waveforms, noise_model)
    shunt_thermal_psd_a2_per_hz = np.full_like(
        waveforms.raw_signal_current_a,
        _shunt_thermal_psd_a2_per_hz(waveforms, noise_model),
    )

    photo_variance, photo_covariance = _propagate_white_psd(
        photocurrent_shot_psd_a2_per_hz,
        pd_tia_impulse,
        sample_rate_hz,
    )
    dark_variance, dark_covariance = _propagate_white_psd(
        dark_current_shot_psd_a2_per_hz,
        pd_tia_impulse,
        sample_rate_hz,
    )
    tia_variance, tia_covariance = _propagate_white_psd(
        tia_psd_a2_per_hz,
        tia_impulse,
        sample_rate_hz,
    )
    rin_variance, rin_covariance = _propagate_white_psd(
        rin_psd_a2_per_hz,
        pd_tia_impulse,
        sample_rate_hz,
    )
    thermal_variance, thermal_covariance = _propagate_white_psd(
        shunt_thermal_psd_a2_per_hz,
        pd_tia_impulse,
        sample_rate_hz,
    )

    total_variance_v2 = (
        photo_variance + dark_variance + tia_variance + rin_variance + thermal_variance
    )
    total_covariance_v2 = (
        photo_covariance
        + dark_covariance
        + tia_covariance
        + rin_covariance
        + thermal_covariance
    )
    pd_tia_dc_gain = waveforms.pd_current_response.dc_gain * tia_response.dc_gain
    return TiaOutputNoiseResult(
        photocurrent_shot_variance_v2=photo_variance,
        dark_current_shot_variance_v2=dark_variance,
        tia_variance_v2=tia_variance,
        rin_variance_v2=rin_variance,
        shunt_thermal_variance_v2=thermal_variance,
        total_variance_v2=total_variance_v2,
        total_lag1_covariance_v2=total_covariance_v2,
        total_sigma_v=np.sqrt(total_variance_v2),
        pd_tia_noise_bandwidth_hz=_noise_bandwidth_from_impulse(
            pd_tia_impulse,
            sample_rate_hz,
            pd_tia_dc_gain,
        ),
        tia_noise_bandwidth_hz=_noise_bandwidth_from_impulse(
            tia_impulse,
            sample_rate_hz,
            tia_response.dc_gain,
        ),
        impulse_tolerance=impulse_tolerance,
    )


def _rin_current_psd_a2_per_hz(
    waveforms: ReceiverWaveformResult,
    noise_model: TimeDomainNoiseModel,
) -> np.ndarray:
    if noise_model.rin_db_per_hz is None:
        return np.zeros_like(waveforms.raw_signal_current_a)

    pd = waveforms.photodiode
    power_w = waveforms.optical_power_w
    if pd.saturation_power_w is None:
        rin_amplitude_a = pd.responsivity_a_per_w * power_w
    else:
        normalized_power = power_w / pd.saturation_power_w
        local_slope_a_per_w = pd.responsivity_a_per_w * (
            1 - np.tanh(normalized_power) ** 2
        )
        rin_amplitude_a = local_slope_a_per_w * power_w
    return rin_amplitude_a**2 * rin_db_per_hz_to_linear(noise_model.rin_db_per_hz)


def _shunt_thermal_psd_a2_per_hz(
    waveforms: ReceiverWaveformResult,
    noise_model: TimeDomainNoiseModel,
) -> float:
    pd = waveforms.photodiode
    if not noise_model.include_shunt_thermal_noise or pd.shunt_resistance_ohm is None:
        return 0.0
    return 4 * K_B * pd.temperature_k / pd.shunt_resistance_ohm


def _propagate_white_psd(
    one_sided_psd_per_hz: np.ndarray,
    impulse_response: np.ndarray,
    sample_rate_hz: float,
) -> tuple[np.ndarray, np.ndarray]:
    psd = readonly_float_array(one_sided_psd_per_hz, "one_sided_psd_per_hz")
    impulse = readonly_float_array(impulse_response, "impulse_response")
    if np.any(psd < 0):
        msg = "one_sided_psd_per_hz must be non-negative."
        raise ValueError(msg)
    if not isfinite(sample_rate_hz) or sample_rate_hz <= 0:
        msg = "sample_rate_hz must be finite and positive."
        raise ValueError(msg)

    memory = impulse.size - 1
    if memory:
        prefix_indices = np.arange(-memory, 0) % psd.size
        extended_psd = np.concatenate((psd[prefix_indices], psd))
    else:
        extended_psd = psd
    input_variance = extended_psd * sample_rate_hz / 2

    variance_extended = lfilter(impulse**2, np.array([1.0]), input_variance)
    variance = np.asarray(variance_extended[memory : memory + psd.size], dtype=float)

    if impulse.size == 1:
        covariance = np.zeros_like(variance)
    else:
        lag1_kernel = impulse[:-1] * impulse[1:]
        covariance_extended = lfilter(
            lag1_kernel,
            np.array([1.0]),
            input_variance,
        )
        covariance = np.asarray(
            covariance_extended[memory : memory + psd.size],
            dtype=float,
        )
    return variance, covariance


def _noise_bandwidth_from_impulse(
    impulse_response: np.ndarray,
    sample_rate_hz: float,
    dc_gain: float,
) -> float:
    if dc_gain == 0:
        msg = "noise bandwidth requires a non-zero DC gain."
        raise ValueError(msg)
    return float(
        sample_rate_hz
        / 2
        * np.sum(np.asarray(impulse_response) ** 2)
        / abs(dc_gain) ** 2
    )
