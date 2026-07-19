"""Example: generate a physically separated PD-current and TIA-voltage eye."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, dbm_to_watt, nrz_levels_from_oma_er
from oma_ber.plotting import plot_deterministic_eye_analysis
from oma_ber.time_domain import (
    TimeGrid,
    analyze_deterministic_eye,
    first_order_lowpass_transfer,
    prbs_bits,
    simulate_pd_tia_waveform,
)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_save_path = repo_root / "tmp" / "physical_rx_eye.png"

    parser = argparse.ArgumentParser(
        description="Plot deterministic PD-current and TIA-voltage eyes.",
    )
    parser.add_argument(
        "--save",
        nargs="?",
        const=default_save_path,
        type=Path,
        help="Optional PNG path. Defaults to ./tmp/physical_rx_eye.png.",
    )
    parser.add_argument(
        "--show", action="store_true", help="Show the matplotlib window."
    )
    args = parser.parse_args()

    symbol_rate_baud = 25e9
    samples_per_symbol = 16
    pd_bandwidth_3db_hz = 14e9
    tia_bandwidth_3db_hz = 18e9
    tia_transimpedance_ohm = 1.5e3
    oma_dbm = -12.0
    er_db = 6.0

    time_grid = TimeGrid(
        symbol_rate_baud=symbol_rate_baud,
        samples_per_symbol=samples_per_symbol,
    )
    bits = prbs_bits(order=7, num_bits=2**7 - 1)
    levels = nrz_levels_from_oma_er(oma_w=dbm_to_watt(oma_dbm), er_db=er_db)
    pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=1e-9,
        bandwidth_3db_hz=pd_bandwidth_3db_hz,
    )
    pd_response = first_order_lowpass_transfer(
        sample_rate_hz=time_grid.sample_rate_hz,
        bandwidth_3db_hz=pd_bandwidth_3db_hz,
        response_kind="dimensionless",
    )
    tia_response = first_order_lowpass_transfer(
        sample_rate_hz=time_grid.sample_rate_hz,
        bandwidth_3db_hz=tia_bandwidth_3db_hz,
        dc_gain=tia_transimpedance_ohm,
        response_kind="transimpedance_ohm",
    )
    receiver_waveforms = simulate_pd_tia_waveform(
        bits=bits,
        p0_w=levels.p0_w,
        p1_w=levels.p1_w,
        time_grid=time_grid,
        pd=pd,
        pd_current_response=pd_response,
        tia_transimpedance_response=tia_response,
    )
    current_eye = analyze_deterministic_eye(
        waveform=receiver_waveforms.pd_output_current_a,
        bits=bits,
        time_grid=time_grid,
        amplitude_unit="A",
    )
    if receiver_waveforms.tia_output_voltage_v is None:
        msg = "TIA response did not produce a voltage waveform."
        raise RuntimeError(msg)
    voltage_eye = analyze_deterministic_eye(
        waveform=receiver_waveforms.tia_output_voltage_v,
        bits=bits,
        time_grid=time_grid,
        amplitude_unit="V",
    )

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    plot_deterministic_eye_analysis(current_eye, max_traces=127, ax=axes[0])
    axes[0].set_title("PD output current")
    plot_deterministic_eye_analysis(voltage_eye, max_traces=127, ax=axes[1])
    axes[1].set_title("TIA output voltage")
    fig.suptitle("Clock-centered deterministic NRZ/OOK receiver eye")

    if args.save is not None:
        save_path = args.save if args.save.is_absolute() else repo_root / args.save
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        print(f"Saved figure to {save_path}")

    print("Physical RX eye example (deterministic ISI only)")
    print(f"PD-input P0: {levels.p0_w:.6e} W")
    print(f"PD-input P1: {levels.p1_w:.6e} W")
    print(f"Exact sample rate: {time_grid.sample_rate_hz:.6e} Hz")
    print(f"Automatic warm-up: {receiver_waveforms.warmup_symbols} symbols")
    print(
        "Current-eye optimum phase: "
        f"{current_eye.optimum_metrics.sampling_phase_ui:.4f} UI",
    )
    print(
        "Current vertical opening: "
        f"{current_eye.optimum_metrics.vertical_eye_opening:.6e} A",
    )
    print(
        "Voltage-eye optimum phase: "
        f"{voltage_eye.optimum_metrics.sampling_phase_ui:.4f} UI",
    )
    print(
        "Voltage vertical opening: "
        f"{voltage_eye.optimum_metrics.vertical_eye_opening:.6e} V",
    )
    print("Random noise, jitter, CDR, and BER contours are not included in Phase 1.")

    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
