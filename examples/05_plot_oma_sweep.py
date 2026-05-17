"""Example: plot BER and Q versus OMA for the scalar NRZ/OOK MVP."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, Receiver, sweep_oma
from oma_ber.plotting import plot_oma_sweep, plot_q_vs_oma


def build_example_results() -> list[dict]:
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

    return sweep_oma(
        oma_dbm_values=np.linspace(-24.0, -8.0, 33),
        er_db=6.0,
        pd=pd,
        rx=rx,
    )


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_save_path = repo_root / "tmp" / "oma_sweep.png"

    parser = argparse.ArgumentParser(description="Plot example OMA sweep results.")
    parser.add_argument(
        "--save",
        nargs="?",
        const=default_save_path,
        type=Path,
        help="Optional path to save the PNG figure. Defaults to ./tmp/oma_sweep.png.",
    )
    parser.add_argument("--show", action="store_true", help="Show the figure in a matplotlib window.")
    args = parser.parse_args()

    results = build_example_results()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    plot_oma_sweep(results, ax=axes[0])
    plot_q_vs_oma(results, ax=axes[1])
    fig.suptitle("Example NRZ/OOK Scalar Gaussian OMA Sweep")

    if args.save is not None:
        save_path = args.save
        if not save_path.is_absolute():
            save_path = repo_root / save_path
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        print(f"Saved figure to {save_path}")

    if args.show:
        plt.show()

    if args.save is None and not args.show:
        first = results[0]
        last = results[-1]
        print(
            "Generated OMA sweep figure in memory: "
            f"{first['oma_dbm']:.1f} to {last['oma_dbm']:.1f} dBm, "
            f"BER {first['ber']:.3e} to {last['ber']:.3e}."
        )

    if not args.show:
        plt.close(fig)


if __name__ == "__main__":
    main()
