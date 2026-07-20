"""Example: residual RJ/SJ averaging and fixed-threshold BER bathtub."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, dbm_to_watt, nrz_levels_from_oma_er
from oma_ber.plotting import (
    plot_deterministic_eye_analysis,
    plot_jittered_ber_bathtub,
)
from oma_ber.time_domain import (
    ResidualTimingJitter,
    TimeDomainNoiseModel,
    TimeGrid,
    analyze_deterministic_eye,
    analyze_jittered_tia_eye,
    analyze_statistical_tia_eye,
    calculate_tia_output_noise,
    first_order_lowpass_transfer,
    prbs_bits,
    simulate_pd_tia_waveform,
)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_save_path = repo_root / "tmp" / "residual_jitter_bathtub.png"
    parser = argparse.ArgumentParser(
        description="Plot fixed-threshold BER under residual RJ and sinusoidal jitter.",
    )
    parser.add_argument(
        "--save",
        nargs="?",
        const=default_save_path,
        type=Path,
        help="Optional PNG path. Defaults to ./tmp/residual_jitter_bathtub.png.",
    )
    parser.add_argument(
        "--show", action="store_true", help="Show the matplotlib window."
    )
    args = parser.parse_args()

    time_grid = TimeGrid(symbol_rate_baud=25e9, samples_per_symbol=16)
    bits = prbs_bits(order=7, num_bits=2**7 - 1)
    levels = nrz_levels_from_oma_er(dbm_to_watt(-18.0), er_db=6.0)
    pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=2e-9,
        shunt_resistance_ohm=1e8,
        temperature_k=300.0,
    )
    pd_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        bandwidth_3db_hz=14e9,
        response_kind="dimensionless",
    )
    tia_response = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        bandwidth_3db_hz=18e9,
        dc_gain=1.5e3,
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
    # The causal PD/TIA latency is separated from the phase within one UI so
    # the bathtub is centered near 0.5 UI. This changes the phase coordinate,
    # not the physical sampled waveform.
    decision_delay_s = 0.4375 * time_grid.unit_interval_s
    phase2 = analyze_statistical_tia_eye(
        waveforms,
        noise,
        decision_delay_s=decision_delay_s,
    )
    residual_jitter = ResidualTimingJitter(
        random_jitter_rms_s=2.5e-12,
        sinusoidal_jitter_peak_s=4.0e-12,
    )
    phase4 = analyze_jittered_tia_eye(
        waveforms,
        noise,
        residual_jitter,
        decision_delay_s=decision_delay_s,
        random_quadrature_points=9,
        sinusoidal_phase_points=16,
        threshold_grid_points=129,
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    plot_deterministic_eye_analysis(deterministic_eye, max_traces=127, ax=axes[0])
    axes[0].set_title("Deterministic TIA eye")
    plot_jittered_ber_bathtub(phase4, ax=axes[1])
    phase2_phases = np.array(
        [result.sampling_phase_ui for result in phase2.phase_results]
    )
    phase2_ber = np.array(
        [max(result.ber, np.finfo(float).tiny) for result in phase2.phase_results]
    )
    axes[1].semilogy(
        phase2_phases,
        phase2_ber,
        color="C2",
        marker=".",
        label="Phase 2: no timing jitter",
    )
    axes[1].legend()
    axes[1].set_title("Fixed-threshold residual-jitter bathtub")
    fig.suptitle("Phase 4 residual timing-jitter analysis")

    if args.save is not None:
        save_path = args.save if args.save.is_absolute() else repo_root / args.save
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        print(f"Saved figure to {save_path}")

    optimum = phase4.optimum_result
    print("Phase 4 residual timing-jitter example")
    print("Inputs are residual sampler jitter after any clock tracking.")
    print(f"Unit interval: {time_grid.unit_interval_s * 1e12:.6f} ps")
    print(f"Residual RJ RMS: {residual_jitter.random_jitter_rms_s * 1e12:.6f} ps")
    print(
        "Residual sinusoidal jitter peak: "
        f"{residual_jitter.sinusoidal_jitter_peak_s * 1e12:.6f} ps"
    )
    print(f"Residual RJ RMS: {phase4.random_jitter_rms_ui:.6f} UI")
    print(f"Residual SJ peak: {phase4.sinusoidal_jitter_peak_ui:.6f} UI")
    print(f"Quadrature nodes: {optimum.quadrature_nodes}")
    print(f"Optimum nominal phase: {optimum.nominal_sampling_phase_ui:.6f} UI")
    print(f"Fixed optimum threshold: {optimum.threshold_v:.6e} V")
    print(f"Jitter-averaged BER: {optimum.ber:.6e}")
    print(
        "The voltage threshold is fixed across all instantaneous jitter nodes at "
        "each nominal phase."
    )
    print(
        "CDR transfer, dual-Dirac extrapolation, and standards masks are not included."
    )

    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
