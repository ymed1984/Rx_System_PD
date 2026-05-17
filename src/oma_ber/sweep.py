"""Top-level OMA-to-BER calculations and OMA sweeps."""

from typing import Any

import numpy as np
from scipy.optimize import brentq

from oma_ber.ber import ber_from_gaussian_levels
from oma_ber.modulation import nrz_levels_from_oma_er
from oma_ber.noise import total_noise_rms_a
from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.units import dbm_to_watt


def calculate_ber_from_oma(
    oma_dbm: float,
    er_db: float,
    pd: Photodiode,
    rx: Receiver,
    optimize_threshold: bool = True,
) -> dict[str, Any]:
    """Calculate NRZ/OOK BER from OMA in dBm and ER in dB.

    The returned optical powers are in W, currents and RMS noise values are in A,
    and threshold_a is the receiver decision threshold in A.
    """
    oma_w = dbm_to_watt(oma_dbm)
    levels = nrz_levels_from_oma_er(oma_w=oma_w, er_db=er_db)

    i0_a = pd.responsivity_a_per_w * levels.p0_w
    i1_a = pd.responsivity_a_per_w * levels.p1_w
    delta_i_a = i1_a - i0_a

    sigma0_a = total_noise_rms_a(
        optical_power_w=levels.p0_w,
        photocurrent_a=i0_a,
        pd=pd,
        rx=rx,
    )
    sigma1_a = total_noise_rms_a(
        optical_power_w=levels.p1_w,
        photocurrent_a=i1_a,
        pd=pd,
        rx=rx,
    )
    ber_result = ber_from_gaussian_levels(
        mu0_a=i0_a,
        mu1_a=i1_a,
        sigma0_a=sigma0_a,
        sigma1_a=sigma1_a,
        optimize_threshold=optimize_threshold,
    )

    return {
        "oma_dbm": oma_dbm,
        "oma_w": oma_w,
        "er_db": er_db,
        "p0_w": levels.p0_w,
        "p1_w": levels.p1_w,
        "pavg_w": levels.pavg_w,
        "i0_a": i0_a,
        "i1_a": i1_a,
        "delta_i_a": delta_i_a,
        "sigma0_a": sigma0_a,
        "sigma1_a": sigma1_a,
        "threshold_a": ber_result["threshold_a"],
        "q_rx": ber_result["q_rx"],
        "ber": ber_result["ber"],
        "ber_from_q": ber_result["ber_from_q"],
    }


def sweep_oma(
    oma_dbm_values: np.ndarray,
    er_db: float,
    pd: Photodiode,
    rx: Receiver,
) -> list[dict[str, Any]]:
    """Calculate BER results for an array of OMA values in dBm."""
    return [
        calculate_ber_from_oma(
            oma_dbm=float(oma_dbm),
            er_db=er_db,
            pd=pd,
            rx=rx,
        )
        for oma_dbm in oma_dbm_values
    ]


def required_oma_dbm(
    target_ber: float,
    er_db: float,
    pd: Photodiode,
    rx: Receiver,
    oma_dbm_min: float = -40.0,
    oma_dbm_max: float = 10.0,
) -> float:
    """Find the OMA in dBm required to reach target_ber."""
    if not 0 < target_ber < 0.5:
        msg = "target_ber must be between 0 and 0.5."
        raise ValueError(msg)
    if oma_dbm_min >= oma_dbm_max:
        msg = "oma_dbm_min must be smaller than oma_dbm_max."
        raise ValueError(msg)

    ber_at_min = calculate_ber_from_oma(oma_dbm_min, er_db, pd, rx)["ber"]
    ber_at_max = calculate_ber_from_oma(oma_dbm_max, er_db, pd, rx)["ber"]

    if ber_at_min <= target_ber:
        return oma_dbm_min
    if ber_at_max > target_ber:
        msg = (
            f"target_ber={target_ber:g} cannot be reached between "
            f"{oma_dbm_min:g} dBm and {oma_dbm_max:g} dBm; "
            f"BER at max OMA is {ber_at_max:g}."
        )
        raise ValueError(msg)

    def residual(oma_dbm: float) -> float:
        return calculate_ber_from_oma(oma_dbm, er_db, pd, rx)["ber"] - target_ber

    return float(brentq(residual, oma_dbm_min, oma_dbm_max))
