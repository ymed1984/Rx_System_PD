"""OMA-to-BER and photodiode receiver analysis tools."""

from oma_ber.ber import (
    ber_for_threshold,
    ber_from_gaussian_levels,
    ber_from_q,
    optimum_threshold,
    q_from_levels,
    qfunc,
)
from oma_ber.bandwidth import apply_oma_penalty_db, bandwidth_penalty_db, first_order_lowpass_mag
from oma_ber.eye import eye_traces, eye_unit_interval_axis
from oma_ber.isi import (
    apply_lti_filter,
    first_order_lowpass_impulse_response,
    sample_at_symbol_centers,
    sampled_eye_levels,
)
from oma_ber.link_budget import pd_input_oma_dbm
from oma_ber.modulation import NRZLevels, nrz_levels_from_oma_er
from oma_ber.noise import (
    K_B,
    Q_E,
    rin_noise_rms_a,
    shot_noise_rms_a,
    thermal_noise_rms_a,
    tia_noise_rms_a,
    total_noise_rms_a,
)
from oma_ber.photodiode import Photodiode
from oma_ber.photodiode import photocurrent_a
from oma_ber.receiver import Receiver
from oma_ber.saturation import (
    compressed_responsivity_a_per_w,
    compressed_responsivity_rational_a_per_w,
    saturated_photocurrent_a,
    saturated_photocurrent_exponential_a,
    saturated_photocurrent_rational_a,
    saturated_photocurrent_soft_clip_a,
    saturated_photocurrent_tanh_a,
)
from oma_ber.sparameters import frequency_response_from_csv, impulse_response_from_frequency_response
from oma_ber.sweep import calculate_ber_from_oma, required_oma_dbm, sweep_oma
from oma_ber.units import (
    db_to_linear,
    dbm_to_watt,
    linear_to_db,
    rin_db_per_hz_to_linear,
    watt_to_dbm,
)
from oma_ber.waveform import nrz_bits_to_levels, prbs_bits, samples_per_symbol

__all__ = [
    "NRZLevels",
    "Photodiode",
    "K_B",
    "Q_E",
    "Receiver",
    "apply_oma_penalty_db",
    "apply_lti_filter",
    "ber_for_threshold",
    "ber_from_gaussian_levels",
    "ber_from_q",
    "calculate_ber_from_oma",
    "compressed_responsivity_a_per_w",
    "compressed_responsivity_rational_a_per_w",
    "db_to_linear",
    "dbm_to_watt",
    "eye_traces",
    "eye_unit_interval_axis",
    "bandwidth_penalty_db",
    "first_order_lowpass_impulse_response",
    "first_order_lowpass_mag",
    "frequency_response_from_csv",
    "impulse_response_from_frequency_response",
    "linear_to_db",
    "nrz_levels_from_oma_er",
    "nrz_bits_to_levels",
    "optimum_threshold",
    "pd_input_oma_dbm",
    "photocurrent_a",
    "prbs_bits",
    "q_from_levels",
    "qfunc",
    "rin_noise_rms_a",
    "rin_db_per_hz_to_linear",
    "required_oma_dbm",
    "shot_noise_rms_a",
    "samples_per_symbol",
    "sample_at_symbol_centers",
    "sampled_eye_levels",
    "saturated_photocurrent_a",
    "saturated_photocurrent_exponential_a",
    "saturated_photocurrent_rational_a",
    "saturated_photocurrent_soft_clip_a",
    "saturated_photocurrent_tanh_a",
    "sweep_oma",
    "thermal_noise_rms_a",
    "tia_noise_rms_a",
    "total_noise_rms_a",
    "watt_to_dbm",
]
