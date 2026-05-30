"""Example: apply a simple first-order bandwidth penalty to OMA."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, Receiver, calculate_ber_from_oma
from oma_ber.bandwidth import apply_oma_penalty_db, bandwidth_penalty_db


def main() -> None:
    pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=1e-9,
        bandwidth_3db_hz=40e9,
    )
    rx = Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )

    oma_dbm = -14.0
    er_db = 6.0
    signal_frequency_hz = 20e9
    bandwidth_3db_hz = 25e9

    penalty_db = bandwidth_penalty_db(
        signal_frequency_hz=signal_frequency_hz,
        bandwidth_3db_hz=bandwidth_3db_hz,
    )
    effective_oma_dbm = apply_oma_penalty_db(oma_dbm=oma_dbm, penalty_db=penalty_db)

    baseline = calculate_ber_from_oma(oma_dbm=oma_dbm, er_db=er_db, pd=pd, rx=rx)
    penalized = calculate_ber_from_oma(oma_dbm=effective_oma_dbm, er_db=er_db, pd=pd, rx=rx)

    print("Example first-order bandwidth penalty")
    print(f"Input OMA: {oma_dbm:.2f} dBm")
    print(f"Signal frequency: {signal_frequency_hz:.3e} Hz")
    print(f"3 dB bandwidth: {bandwidth_3db_hz:.3e} Hz")
    print(f"Bandwidth penalty: {penalty_db:.3f} dB")
    print(f"Effective OMA: {effective_oma_dbm:.2f} dBm")
    print(f"Baseline BER: {baseline['ber']:.3e}")
    print(f"Penalized BER: {penalized['ber']:.3e}")


if __name__ == "__main__":
    main()
