"""OMA-to-BER and photodiode receiver analysis tools."""

from oma_ber.ber import (
    ber_for_threshold,
    ber_from_gaussian_levels,
    ber_from_q,
    optimum_threshold,
    q_from_levels,
    qfunc,
)
from oma_ber.link_budget import pd_input_oma_dbm
from oma_ber.modulation import NRZLevels, nrz_levels_from_oma_er
from oma_ber.noise import (
    Q_E,
    rin_noise_rms_a,
    shot_noise_rms_a,
    tia_noise_rms_a,
    total_noise_rms_a,
)
from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.sweep import calculate_ber_from_oma, required_oma_dbm, sweep_oma
from oma_ber.units import (
    db_to_linear,
    dbm_to_watt,
    linear_to_db,
    rin_db_per_hz_to_linear,
    watt_to_dbm,
)

__all__ = [
    "NRZLevels",
    "Photodiode",
    "Q_E",
    "Receiver",
    "ber_for_threshold",
    "ber_from_gaussian_levels",
    "ber_from_q",
    "calculate_ber_from_oma",
    "db_to_linear",
    "dbm_to_watt",
    "linear_to_db",
    "nrz_levels_from_oma_er",
    "optimum_threshold",
    "pd_input_oma_dbm",
    "q_from_levels",
    "qfunc",
    "rin_noise_rms_a",
    "rin_db_per_hz_to_linear",
    "required_oma_dbm",
    "shot_noise_rms_a",
    "sweep_oma",
    "tia_noise_rms_a",
    "total_noise_rms_a",
    "watt_to_dbm",
]
