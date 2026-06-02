"""Example: plot deterministic NRZ/OOK eye diagrams with bandwidth-limited ISI."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import dbm_to_watt, nrz_levels_from_oma_er
from oma_ber.isi import (
    apply_lti_filter,
    first_order_lowpass_impulse_response,
    sample_at_symbol_centers,
    sampled_eye_levels,
)
from oma_ber.plotting import plot_eye_diagram
from oma_ber.waveform import nrz_bits_to_levels, prbs_bits, samples_per_symbol


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_save_path = repo_root / "tmp" / "eye_diagram.png"

    parser = argparse.ArgumentParser(description="Plot deterministic NRZ/OOK eye diagrams.")
    parser.add_argument(
        "--save",
        nargs="?",
        const=default_save_path,
        type=Path,
        help="Optional path to save the PNG figure. Defaults to ./tmp/eye_diagram.png.",
    )
    parser.add_argument("--show", action="store_true", help="Show the figure in a matplotlib window.")
    args = parser.parse_args()

    sample_rate_hz = 200e9
    symbol_rate_baud = 25e9
    bandwidth_3db_hz = 14e9
    oma_dbm = -12.0
    er_db = 6.0

    sps = samples_per_symbol(sample_rate_hz, symbol_rate_baud)
    bits = prbs_bits(num_bits=256, seed=12)
    levels = nrz_levels_from_oma_er(oma_w=dbm_to_watt(oma_dbm), er_db=er_db)
    waveform_w = nrz_bits_to_levels(bits, levels.p0_w, levels.p1_w, sps)
    impulse_response = first_order_lowpass_impulse_response(
        sample_rate_hz=sample_rate_hz,
        bandwidth_3db_hz=bandwidth_3db_hz,
        num_taps=96,
    )
    filtered_w = apply_lti_filter(waveform_w, impulse_response)

    ideal_samples_w = sample_at_symbol_centers(waveform_w, sps)
    filtered_samples_w = sample_at_symbol_centers(filtered_w, sps)
    ideal_metrics = sampled_eye_levels(ideal_samples_w, bits)
    filtered_metrics = sampled_eye_levels(filtered_samples_w, bits)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    plot_eye_diagram(
        waveform_w,
        samples_per_symbol=sps,
        sample_rate_hz=sample_rate_hz,
        max_traces=120,
        eye_metrics=ideal_metrics,
        ax=axes[0],
        ylabel="Optical power [W]",
    )
    axes[0].set_title("Ideal NRZ/OOK")
    plot_eye_diagram(
        filtered_w,
        samples_per_symbol=sps,
        sample_rate_hz=sample_rate_hz,
        max_traces=120,
        eye_metrics=filtered_metrics,
        ax=axes[1],
        ylabel="Optical power [W]",
    )
    axes[1].set_title(f"First-order LPF, f3dB={bandwidth_3db_hz / 1e9:.1f} GHz")
    fig.suptitle(f"Deterministic 2 UI Eye, OMA={oma_dbm:.1f} dBm, ER={er_db:.1f} dB")

    if args.save is not None:
        save_path = args.save
        if not save_path.is_absolute():
            save_path = repo_root / save_path
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        print(f"Saved figure to {save_path}")

    print("Sampled-eye metrics")
    print(f"Ideal vertical eye opening: {ideal_metrics['vertical_eye_opening']:.3e} W")
    print(f"Filtered vertical eye opening: {filtered_metrics['vertical_eye_opening']:.3e} W")
    print(f"Filtered ISI penalty indicator: {filtered_metrics['isi_penalty_db']:.3f} dB")

    if args.show:
        plt.show()

    if not args.show:
        plt.close(fig)


if __name__ == "__main__":
    main()
