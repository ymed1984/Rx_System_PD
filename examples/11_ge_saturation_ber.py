"""Example: calculate NRZ/OOK BER with simplified tanh Ge PD saturation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, Receiver, calculate_ber_from_oma


def main() -> None:
    linear_pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=1e-9,
        bandwidth_3db_hz=40e9,
    )
    ge_saturated_pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=1e-9,
        bandwidth_3db_hz=40e9,
        saturation_power_w=1e-4,
    )
    rx = Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )

    linear = calculate_ber_from_oma(
        oma_dbm=-5.0,
        er_db=6.0,
        pd=linear_pd,
        rx=rx,
    )
    saturated = calculate_ber_from_oma(
        oma_dbm=-5.0,
        er_db=6.0,
        pd=ge_saturated_pd,
        rx=rx,
    )

    print("Example simplified tanh Ge PD saturation in BER calculation")
    print(f"Linear Delta_I: {linear['delta_i_a']:.3e} A")
    print(f"Saturated Delta_I: {saturated['delta_i_a']:.3e} A")
    print(f"Current OMA compression: {saturated['oma_current_compression_db']:.3f} dB")
    print(f"Linear BER: {linear['ber']:.3e}")
    print(f"Saturated BER: {saturated['ber']:.3e}")


if __name__ == "__main__":
    main()
