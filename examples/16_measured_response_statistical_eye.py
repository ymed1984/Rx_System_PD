"""Example: compare analytic and synthetic measured-like PD/TIA responses."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import Photodiode, dbm_to_watt, nrz_levels_from_oma_er
from oma_ber.plotting import (
    plot_deterministic_eye_analysis,
    plot_statistical_ber_comparison,
    plot_transfer_magnitude_comparison,
)
from oma_ber.time_domain import (
    DeterministicEyeAnalysis,
    DiscreteTransferFunction,
    MeasuredFrequencyResponse,
    StatisticalEyeAnalysis,
    TiaOutputNoiseResult,
    TimeDomainNoiseModel,
    TimeGrid,
    analyze_deterministic_eye,
    analyze_statistical_tia_eye,
    calculate_tia_output_noise,
    fir_transfer_from_measured_response,
    first_order_lowpass_transfer,
    prbs_bits,
    simulate_pd_tia_waveform,
)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_save_path = repo_root / "tmp" / "measured_response_statistical_eye.png"
    parser = argparse.ArgumentParser(
        description="Compare one-pole and synthetic measured-like RX responses.",
    )
    parser.add_argument(
        "--save",
        nargs="?",
        const=default_save_path,
        type=Path,
        help="Optional PNG path. Defaults to ./tmp/measured_response_statistical_eye.png.",
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
    noise_model = TimeDomainNoiseModel(
        tia_input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150.0,
    )

    analytic_pd = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        bandwidth_3db_hz=14e9,
        response_kind="dimensionless",
    )
    analytic_tia = first_order_lowpass_transfer(
        time_grid.sample_rate_hz,
        bandwidth_3db_hz=18e9,
        dc_gain=1.5e3,
        response_kind="transimpedance_ohm",
    )

    # These causal second-order responses stand in for imported measured data.
    # Replace MeasuredFrequencyResponse below with
    # measured_frequency_response_from_csv() when calibrated CSV data exists.
    measured_source_pd = _butterworth_transfer(
        sample_rate_hz=time_grid.sample_rate_hz,
        bandwidth_3db_hz=14e9,
        dc_gain=1.0,
        response_kind="dimensionless",
    )
    measured_source_tia = _butterworth_transfer(
        sample_rate_hz=time_grid.sample_rate_hz,
        bandwidth_3db_hz=18e9,
        dc_gain=1.5e3,
        response_kind="transimpedance_ohm",
    )
    fft_size = 4096
    fir_length_samples = 512
    measured_pd, pd_diagnostics = fir_transfer_from_measured_response(
        _sample_as_measured_response(
            measured_source_pd,
            fft_size,
            "PD optical conversion output -> TIA input",
        ),
        time_grid.sample_rate_hz,
        fft_size,
        fir_length_samples,
        causality_tolerance=1e-10,
        tail_energy_tolerance=1e-10,
    )
    measured_tia, tia_diagnostics = fir_transfer_from_measured_response(
        _sample_as_measured_response(
            measured_source_tia,
            fft_size,
            "TIA input current -> TIA output voltage",
        ),
        time_grid.sample_rate_hz,
        fft_size,
        fir_length_samples,
        causality_tolerance=1e-10,
        tail_energy_tolerance=1e-10,
    )

    # Use one explicit time reference for both responses. This separates their
    # causal latency from the phase coordinate without independently centering
    # either result.
    decision_delay_s = 0.4375 * time_grid.unit_interval_s
    analytic_eye, analytic_statistical, analytic_noise = _analyze_case(
        bits,
        levels.p0_w,
        levels.p1_w,
        time_grid,
        pd,
        analytic_pd,
        analytic_tia,
        noise_model,
        decision_delay_s,
    )
    measured_eye, measured_statistical, measured_noise = _analyze_case(
        bits,
        levels.p0_w,
        levels.p1_w,
        time_grid,
        pd,
        measured_pd,
        measured_tia,
        noise_model,
        decision_delay_s,
    )

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    plot_transfer_magnitude_comparison(
        {"one-pole": analytic_pd, "measured-like FIR": measured_pd},
        maximum_frequency_hz=60e9,
        ax=axes[0, 0],
    )
    axes[0, 0].set_title("PD current response")
    plot_transfer_magnitude_comparison(
        {"one-pole": analytic_tia, "measured-like FIR": measured_tia},
        maximum_frequency_hz=60e9,
        ax=axes[0, 1],
    )
    axes[0, 1].set_title("TIA transimpedance response")
    plot_deterministic_eye_analysis(measured_eye, max_traces=127, ax=axes[1, 0])
    axes[1, 0].set_title("Measured-like deterministic TIA eye")
    plot_statistical_ber_comparison(
        {
            "one-pole": analytic_statistical,
            "measured-like FIR": measured_statistical,
        },
        ax=axes[1, 1],
    )
    axes[1, 1].set_title("Response-dependent BER")
    fig.suptitle("Phase 3 complex-response integration")

    if args.save is not None:
        save_path = args.save if args.save.is_absolute() else repo_root / args.save
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
        print(f"Saved figure to {save_path}")

    print("Phase 3 measured-response integration example")
    print(
        "Synthetic second-order responses are used; this is not measured hardware data."
    )
    print(f"PD FIR pre-echo energy ratio: {pd_diagnostics.pre_echo_energy_ratio:.3e}")
    print(
        "PD FIR discarded-tail energy ratio: "
        f"{pd_diagnostics.discarded_tail_energy_ratio:.3e}"
    )
    print(f"TIA FIR pre-echo energy ratio: {tia_diagnostics.pre_echo_energy_ratio:.3e}")
    print(
        "TIA FIR discarded-tail energy ratio: "
        f"{tia_diagnostics.discarded_tail_energy_ratio:.3e}"
    )
    _print_case("one-pole", analytic_statistical, analytic_noise)
    _print_case("measured-like FIR", measured_statistical, measured_noise)
    print("Jitter, CDR, de-embedding, and standards masks are not included.")

    if args.show:
        plt.show()
    else:
        plt.close(fig)


def _butterworth_transfer(
    sample_rate_hz: float,
    bandwidth_3db_hz: float,
    dc_gain: float,
    response_kind: str,
) -> DiscreteTransferFunction:
    numerator, denominator = butter(
        2,
        bandwidth_3db_hz,
        btype="lowpass",
        fs=sample_rate_hz,
    )
    numerator = numerator * (dc_gain / (np.sum(numerator) / np.sum(denominator)))
    return DiscreteTransferFunction(
        numerator=numerator,
        denominator=denominator,
        sample_rate_hz=sample_rate_hz,
        response_kind=response_kind,  # type: ignore[arg-type]
    )


def _sample_as_measured_response(
    transfer: DiscreteTransferFunction,
    fft_size: int,
    reference_plane: str,
) -> MeasuredFrequencyResponse:
    frequency_hz = np.fft.rfftfreq(fft_size, d=1 / transfer.sample_rate_hz)
    z_inverse = np.exp(-1j * 2 * np.pi * frequency_hz / transfer.sample_rate_hz)
    numerator = np.polynomial.polynomial.polyval(z_inverse, transfer.numerator)
    denominator = np.polynomial.polynomial.polyval(z_inverse, transfer.denominator)
    response_complex = numerator / denominator
    response_complex[-1] = response_complex[-1].real
    return MeasuredFrequencyResponse(
        frequency_hz=frequency_hz,
        response_complex=response_complex,
        response_kind=transfer.response_kind,
        reference_plane=reference_plane,
    )


def _analyze_case(
    bits: np.ndarray,
    p0_w: float,
    p1_w: float,
    time_grid: TimeGrid,
    pd: Photodiode,
    pd_response: DiscreteTransferFunction,
    tia_response: DiscreteTransferFunction,
    noise_model: TimeDomainNoiseModel,
    decision_delay_s: float,
) -> tuple[DeterministicEyeAnalysis, StatisticalEyeAnalysis, TiaOutputNoiseResult]:
    waveforms = simulate_pd_tia_waveform(
        bits=bits,
        p0_w=p0_w,
        p1_w=p1_w,
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
    noise = calculate_tia_output_noise(waveforms, noise_model)
    statistical_eye = analyze_statistical_tia_eye(
        waveforms,
        noise,
        decision_delay_s=decision_delay_s,
    )
    return deterministic_eye, statistical_eye, noise


def _print_case(
    label: str,
    statistical_eye: StatisticalEyeAnalysis,
    noise: TiaOutputNoiseResult,
) -> None:
    optimum = statistical_eye.optimum_result
    print(f"{label} PD+TIA ENBW: {noise.pd_tia_noise_bandwidth_hz:.6e} Hz")
    print(f"{label} optimum phase: {optimum.sampling_phase_ui:.5f} UI")
    print(f"{label} optimum threshold: {optimum.threshold_v:.6e} V")
    print(f"{label} Gaussian-mixture BER: {optimum.ber:.6e}")


if __name__ == "__main__":
    main()
