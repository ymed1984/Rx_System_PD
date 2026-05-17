"""Example: sweep OMA and print BER values."""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, Receiver, sweep_oma


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

    print("OMA sweep for example NRZ/OOK scalar Gaussian model")
    for result in sweep_oma(np.linspace(-22.0, -8.0, 8), er_db=6.0, pd=pd, rx=rx):
        print(
            f"OMA {result['oma_dbm']:6.2f} dBm  "
            f"Q {result['q_rx']:6.3f}  BER {result['ber']:.3e}"
        )


if __name__ == "__main__":
    main()
