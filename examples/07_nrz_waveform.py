"""Example: create a deterministic NRZ/OOK optical waveform."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import dbm_to_watt, nrz_levels_from_oma_er
from oma_ber.waveform import nrz_bits_to_levels, prbs_bits, samples_per_symbol


def main() -> None:
    sample_rate_hz = 100e9
    symbol_rate_baud = 25e9
    oma_dbm = -12.0
    er_db = 6.0

    sps = samples_per_symbol(
        sample_rate_hz=sample_rate_hz,
        symbol_rate_baud=symbol_rate_baud,
    )
    levels = nrz_levels_from_oma_er(oma_w=dbm_to_watt(oma_dbm), er_db=er_db)
    bits = prbs_bits(num_bits=16, seed=7)
    waveform_w = nrz_bits_to_levels(
        bits=bits,
        low_level=levels.p0_w,
        high_level=levels.p1_w,
        samples_per_symbol=sps,
    )

    print("Example deterministic NRZ/OOK waveform")
    print(f"Sample rate: {sample_rate_hz:.3e} Hz")
    print(f"Symbol rate: {symbol_rate_baud:.3e} baud")
    print(f"Samples per symbol: {sps}")
    print(f"OMA: {oma_dbm:.2f} dBm")
    print(f"ER: {er_db:.2f} dB")
    print(f"Bits: {bits.tolist()}")
    print(f"Waveform length: {len(waveform_w)} samples")
    print(f"First 12 optical power samples: {waveform_w[:12].tolist()} W")


if __name__ == "__main__":
    main()
