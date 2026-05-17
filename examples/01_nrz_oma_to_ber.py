"""Example: calculate NRZ/OOK BER from one OMA value."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, Receiver, calculate_ber_from_oma


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

    result = calculate_ber_from_oma(
        oma_dbm=-10.0,
        er_db=6.0,
        pd=pd,
        rx=rx,
    )

    print("Example NRZ/OOK scalar Gaussian model")
    print(f"OMA: {result['oma_dbm']:.2f} dBm ({result['oma_w']:.3e} W)")
    print(f"P0/P1: {result['p0_w']:.3e} W / {result['p1_w']:.3e} W")
    print(f"I0/I1: {result['i0_a']:.3e} A / {result['i1_a']:.3e} A")
    print(f"Q_rx: {result['q_rx']:.3f}")
    print(f"BER: {result['ber']:.3e}")


if __name__ == "__main__":
    main()
