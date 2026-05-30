"""Example: filter an NRZ/OOK waveform with a CSV frequency response."""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import dbm_to_watt, nrz_levels_from_oma_er
from oma_ber.isi import apply_lti_filter, sample_at_symbol_centers, sampled_eye_levels
from oma_ber.sparameters import impulse_response_from_frequency_response
from oma_ber.waveform import nrz_bits_to_levels, prbs_bits, samples_per_symbol


def main() -> None:
    sample_rate_hz = 100e9
    symbol_rate_baud = 25e9
    bandwidth_3db_hz = 18e9
    oma_dbm = -12.0
    er_db = 6.0

    sps = samples_per_symbol(sample_rate_hz, symbol_rate_baud)
    bits = prbs_bits(num_bits=64, seed=9)
    levels = nrz_levels_from_oma_er(oma_w=dbm_to_watt(oma_dbm), er_db=er_db)
    waveform_w = nrz_bits_to_levels(bits, levels.p0_w, levels.p1_w, sps)

    frequency_hz = np.linspace(0.5e9, 40e9, 128)
    response_complex = 1 / (1 + 1j * frequency_hz / bandwidth_3db_hz)
    impulse_response = impulse_response_from_frequency_response(
        frequency_hz=frequency_hz,
        response_complex=response_complex,
        sample_rate_hz=sample_rate_hz,
        num_taps=128,
    )

    filtered_w = apply_lti_filter(waveform_w, impulse_response)
    sampled_w = sample_at_symbol_centers(filtered_w, sps)
    metrics = sampled_eye_levels(sampled_w, bits)

    print("Example CSV-style frequency-response filtering")
    print(f"Sample rate: {sample_rate_hz:.3e} Hz")
    print(f"Symbol rate: {symbol_rate_baud:.3e} baud")
    print(f"Synthetic 3 dB bandwidth: {bandwidth_3db_hz:.3e} Hz")
    print(f"Frequency-response points: {len(frequency_hz)}")
    print(f"Impulse-response taps: {len(impulse_response)}")
    print(f"Vertical eye opening: {metrics['vertical_eye_opening']:.3e} W")
    print(f"ISI penalty indicator: {metrics['isi_penalty_db']:.3f} dB")


if __name__ == "__main__":
    main()
