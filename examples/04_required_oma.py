"""Example: solve required OMA for a target BER."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, Receiver, calculate_ber_from_oma, required_oma_dbm


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
    target_ber = 1e-6

    oma_dbm = required_oma_dbm(
        target_ber=target_ber,
        er_db=6.0,
        pd=pd,
        rx=rx,
        oma_dbm_min=-40.0,
        oma_dbm_max=0.0,
    )
    result = calculate_ber_from_oma(oma_dbm=oma_dbm, er_db=6.0, pd=pd, rx=rx)

    print("Required OMA for example NRZ/OOK scalar Gaussian model")
    print(f"Target BER: {target_ber:.1e}")
    print(f"Required OMA: {oma_dbm:.2f} dBm")
    print(f"Computed BER: {result['ber']:.3e}")


if __name__ == "__main__":
    main()
