"""Example: sweep photodiode responsivity at a fixed OMA."""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, Receiver, calculate_ber_from_oma


def main() -> None:
    rx = Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )

    print("Responsivity sweep for example NRZ/OOK scalar Gaussian model")
    for responsivity_a_per_w in np.linspace(0.5, 1.0, 6):
        pd = Photodiode(
            responsivity_a_per_w=float(responsivity_a_per_w),
            dark_current_a=1e-9,
            bandwidth_3db_hz=40e9,
        )
        result = calculate_ber_from_oma(
            oma_dbm=-12.0,
            er_db=6.0,
            pd=pd,
            rx=rx,
        )
        print(
            f"Rpd {responsivity_a_per_w:.2f} A/W  "
            f"Q {result['q_rx']:6.3f}  BER {result['ber']:.3e}"
        )


if __name__ == "__main__":
    main()
