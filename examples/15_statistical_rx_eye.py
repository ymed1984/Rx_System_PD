"""Example: PSD-based Gaussian-mixture statistical TIA eye and BER."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, dbm_to_watt, nrz_levels_from_oma_er
from oma_ber.plotting import (
    plot_deterministic_eye_analysis,
    plot_statistical_ber_vs_phase,
    plot_statistical_eye_density,
)
from oma_ber.time_domain import (
    TimeDomainNoiseModel,
    TimeGrid,
    analyze_deterministic_eye,
    analyze_statistical_tia_eye,
    calculate_tia_output_noise,
    first_order_lowpass_transfer,
    prbs_bits,
    simulate_pd_tia_waveform,
)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_save_path = repo_root / "tmp" / "statistical_rx_eye.png"

    parser = argparse.ArgumentParser(
        description="Plot PSD-based Gaussian-mixture statistical RX eye results.",
    )
    parser.add_argument(
        "--save",
        nargs="?",
        const=default_save_path,
        type=Path,
        help="Optional PNG path. Defaults to ./tmp/statistical_rx_eye.png.",
    )
    parser.add_argument(
        "--show", action="store_true", help="Show the matplotlib window."
    )
    args = parser.parse_args()

    symbol_rate_baud = 25e9
    samples_per_symbol = 32
    pd_bandwidth_3db_hz = 14e9
    tia_bandwidth_3db_hz = 18e9
    tia_transimpedance_ohm = 1.5e3
    oma_dbm = -18.0
    er_db = 6.0

    time_grid = TimeGrid(symbol_rate_baud, samples_per_symbol)
    bits = prbs_bits(order=7, num_bits=2**7 - 1)
    levels = nrz_levels_from_oma_er(dbm_to_watt(oma_dbm), er_db)
    pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=2e-9,
        bandwidth_3db_hz=pd_bandwidth_3db_hz,
        shunt_resistance_ohm=1e8,
        temperature_k=300.0,
    )
    pd_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        pd_bandwidth_3db_hz,
        response_kind="dimensionless",
    )
    tia_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        tia_bandwidth_3db_hz,
        dc_gain=tia_transimpedance_ohm,
        response_kind="transimpedance_ohm",
    )
    waveforms = simulate_pd_tia_waveform(
        bits=bits,
        p0_w=levels.p0_w,
        p1_w=levels.p1_w,
        time_grid=time_grid,
        pd=pd,
        pd_current_response=pd_response,
        tia_transimpedance_response=tia_response,
    )
    if waveforms.tia_output_voltage_v is None:
        msg = "TIA response did not produce a voltage waveform."
        raise RuntimeError(msg)
    deterministic_eye = analyze_deterministic_eye(
        waveforms.tia_output_voltage_v,
        bits,
        time_grid,
        amplitude_unit="V",
    )
    noise = calculate_tia_output_noise(
        waveforms,
        TimeDomainNoiseModel(
            tia_input_current_noise_density_a_per_sqrt_hz=10e-12,
            rin_db_per_hz=-150.0,
        ),
    )
    statistical_eye = analyze_statistical_tia_eye(waveforms, noise)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    plot_deterministic_eye_analysis(deterministic_eye, max_traces=127, ax=axes[0])
    axes[0].set_title("Deterministic TIA eye")
    plot_statistical_eye_density(statistical_eye, ax=axes[1])
    axes[1].set_title("Gaussian-mixture density")
    plot_statistical_ber_vs_phase(statistical_eye, ax=axes[2])
    axes[2].set_title("BER-optimized threshold")
    fig.suptitle("PSD-based statistical NRZ/OOK receiver eye")

    if args.save is not None:
        save_path = args.save if args.save.is_absolute() else repo_root / args.save
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        print(f"Saved figure to {save_path}")

    optimum = statistical_eye.optimum_result
    print("Physical RX statistical eye example (Gaussian mixture)")
    print(f"PD-input P0: {levels.p0_w:.6e} W")
    print(f"PD-input P1: {levels.p1_w:.6e} W")
    print(f"PD+TIA noise bandwidth: {noise.pd_tia_noise_bandwidth_hz:.6e} Hz")
    print(f"TIA noise bandwidth: {noise.tia_noise_bandwidth_hz:.6e} Hz")
    print(f"BER-optimum phase: {optimum.sampling_phase_ui:.5f} UI")
    print(f"Optimum threshold: {optimum.threshold_v:.6e} V")
    print(f"Effective Q indicator: {optimum.effective_q:.6f}")
    print(f"Gaussian-mixture BER: {optimum.ber:.6e}")
    print(f"Bit-0 RMS noise: {optimum.rms_zero_noise_v:.6e} V")
    print(f"Bit-1 RMS noise: {optimum.rms_one_noise_v:.6e} V")
    print(
        "Jitter, CDR, optical-field dispersion, and standards masks are not included."
    )

    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
