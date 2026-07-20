from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

from oma_ber import Photodiode
from oma_ber.time_domain import (
    DiscreteTransferFunction,
    MeasuredFrequencyResponse,
    TimeDomainNoiseModel,
    TimeGrid,
    analyze_statistical_tia_eye,
    calculate_tia_output_noise,
    equivalent_noise_bandwidth_hz,
    fir_transfer_from_measured_response,
    first_order_lowpass_transfer,
    measured_frequency_response_from_csv,
    prbs_bits,
    simulate_pd_tia_waveform,
    transfer_impulse_response,
)


def test_measured_frequency_response_loads_csv_with_physical_metadata(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "tia_response.csv"
    csv_path.write_text(
        "frequency_hz,magnitude_db,phase_deg\n"
        "0,63.5218251811,0\n"
        "50000000000,63.5218251811,0\n",
        encoding="utf-8",
    )

    response = measured_frequency_response_from_csv(
        csv_path,
        response_kind="transimpedance_ohm",
        reference_plane="TIA input current -> TIA output voltage",
    )

    assert response.response_kind == "transimpedance_ohm"
    assert response.reference_plane == "TIA input current -> TIA output voltage"
    assert response.source_path == csv_path
    assert response.frequency_hz[0] == 0.0
    np.testing.assert_allclose(response.response_complex, 1500.0, rtol=1e-11)
    assert not response.frequency_hz.flags.writeable
    assert not response.response_complex.flags.writeable


def test_flat_dimensionless_response_converts_to_identity_fir() -> None:
    sample_rate_hz = 100e9
    response = _flat_response(sample_rate_hz, 1.0, "dimensionless")

    transfer, diagnostics = fir_transfer_from_measured_response(
        response,
        sample_rate_hz=sample_rate_hz,
        fft_size=128,
        fir_length_samples=1,
    )

    np.testing.assert_allclose(transfer.numerator, np.array([1.0]), atol=1e-14)
    assert transfer.response_kind == "dimensionless"
    assert diagnostics.measured_dc_gain == pytest.approx(1.0)
    assert diagnostics.fir_dc_gain == pytest.approx(1.0)
    assert diagnostics.dc_gain_relative_error == pytest.approx(0.0)
    assert diagnostics.equivalent_noise_bandwidth_hz == pytest.approx(50e9)
    assert diagnostics.group_delay_s == pytest.approx(0.0)
    assert diagnostics.pre_echo_energy_ratio == pytest.approx(0.0)
    assert diagnostics.discarded_tail_energy_ratio == pytest.approx(0.0)


def test_flat_transimpedance_response_preserves_ohm_gain() -> None:
    sample_rate_hz = 100e9
    response = _flat_response(sample_rate_hz, -1500.0, "transimpedance_ohm")

    transfer, diagnostics = fir_transfer_from_measured_response(
        response,
        sample_rate_hz,
        fft_size=128,
        fir_length_samples=1,
    )

    assert transfer.response_kind == "transimpedance_ohm"
    assert transfer.dc_gain == pytest.approx(-1500.0)
    assert diagnostics.fir_dc_gain == pytest.approx(-1500.0)


def test_synthetic_first_order_response_matches_analytic_transfer() -> None:
    sample_rate_hz = 400e9
    analytic = first_order_lowpass_transfer(sample_rate_hz, 15e9)
    fft_size = 4096
    frequency_hz, response_complex = _sample_transfer_response(analytic, fft_size)
    response = MeasuredFrequencyResponse(
        frequency_hz,
        response_complex,
        "dimensionless",
        "synthetic PD current response",
    )

    measured, diagnostics = fir_transfer_from_measured_response(
        response,
        sample_rate_hz,
        fft_size=fft_size,
        fir_length_samples=512,
        causality_tolerance=1e-12,
        tail_energy_tolerance=1e-12,
    )

    expected_impulse = transfer_impulse_response(analytic)
    np.testing.assert_allclose(
        measured.numerator[: expected_impulse.size],
        expected_impulse,
        rtol=1e-11,
        atol=1e-13,
    )
    assert diagnostics.equivalent_noise_bandwidth_hz == pytest.approx(
        equivalent_noise_bandwidth_hz(analytic),
        rel=1e-11,
    )


def test_known_pure_delay_can_be_removed_explicitly() -> None:
    sample_rate_hz = 100e9
    fft_size = 128
    delay_samples = 4
    delay_s = delay_samples / sample_rate_hz
    frequency_hz = np.fft.rfftfreq(fft_size, d=1 / sample_rate_hz)
    delayed_response = np.exp(-1j * 2 * np.pi * frequency_hz * delay_s)
    response = MeasuredFrequencyResponse(
        frequency_hz,
        delayed_response,
        "dimensionless",
        "known pure delay",
    )

    delayed_transfer, delayed_diagnostics = fir_transfer_from_measured_response(
        response,
        sample_rate_hz,
        fft_size,
        fir_length_samples=delay_samples + 1,
    )
    deembedded_transfer, deembedded_diagnostics = fir_transfer_from_measured_response(
        response,
        sample_rate_hz,
        fft_size,
        fir_length_samples=1,
        removed_reference_delay_s=delay_s,
    )

    expected_delayed = np.zeros(delay_samples + 1)
    expected_delayed[-1] = 1.0
    np.testing.assert_allclose(delayed_transfer.numerator, expected_delayed, atol=1e-14)
    np.testing.assert_allclose(
        deembedded_transfer.numerator, np.array([1.0]), atol=1e-14
    )
    assert delayed_diagnostics.group_delay_s == pytest.approx(delay_s)
    assert deembedded_diagnostics.group_delay_s == pytest.approx(0.0, abs=1e-24)
    assert deembedded_diagnostics.removed_reference_delay_s == pytest.approx(delay_s)


def test_noncausal_pre_echo_is_rejected() -> None:
    sample_rate_hz = 100e9
    fft_size = 128
    frequency_hz = np.fft.rfftfreq(fft_size, d=1 / sample_rate_hz)
    advance_s = 2 / sample_rate_hz
    response = MeasuredFrequencyResponse(
        frequency_hz,
        np.exp(1j * 2 * np.pi * frequency_hz * advance_s),
        "dimensionless",
        "noncausal advance",
    )

    with pytest.raises(ValueError, match="pre_echo_energy_ratio"):
        fir_transfer_from_measured_response(
            response,
            sample_rate_hz,
            fft_size,
            fir_length_samples=1,
        )


def test_discarded_positive_time_tail_is_rejected() -> None:
    sample_rate_hz = 100e9
    fft_size = 128
    frequency_hz = np.fft.rfftfreq(fft_size, d=1 / sample_rate_hz)
    delay_s = 10 / sample_rate_hz
    response = MeasuredFrequencyResponse(
        frequency_hz,
        np.exp(-1j * 2 * np.pi * frequency_hz * delay_s),
        "dimensionless",
        "causal delayed response",
    )

    with pytest.raises(ValueError, match="discarded_tail_energy_ratio"):
        fir_transfer_from_measured_response(
            response,
            sample_rate_hz,
            fft_size,
            fir_length_samples=5,
        )


def test_measured_firs_feed_the_phase2_signal_and_noise_paths() -> None:
    time_grid = TimeGrid(10e9, 8)
    fft_size = 128
    pd_response, _ = fir_transfer_from_measured_response(
        _flat_response(time_grid.sample_rate_hz, 1.0, "dimensionless"),
        time_grid.sample_rate_hz,
        fft_size,
        fir_length_samples=1,
    )
    tia_response, _ = fir_transfer_from_measured_response(
        _flat_response(time_grid.sample_rate_hz, 1e3, "transimpedance_ohm"),
        time_grid.sample_rate_hz,
        fft_size,
        fir_length_samples=1,
    )
    waveforms = simulate_pd_tia_waveform(
        bits=prbs_bits(7, 127),
        p0_w=2e-6,
        p1_w=8e-6,
        time_grid=time_grid,
        pd=Photodiode(0.8, 1e-9),
        pd_current_response=pd_response,
        tia_transimpedance_response=tia_response,
    )

    noise = calculate_tia_output_noise(
        waveforms,
        TimeDomainNoiseModel(8e-12),
    )
    analysis = analyze_statistical_tia_eye(waveforms, noise)

    assert np.isfinite(analysis.optimum_result.ber)
    assert analysis.optimum_result.threshold_v > 0
    assert noise.tia_noise_bandwidth_hz == pytest.approx(
        diagnostics_bandwidth_hz := time_grid.sample_rate_hz / 2
    )
    assert diagnostics_bandwidth_hz > 0


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (
            lambda: MeasuredFrequencyResponse(
                np.array([1e9, 2e9]),
                np.ones(2, dtype=complex),
                "dimensionless",
                "PD response",
            ),
            "include DC",
        ),
        (
            lambda: MeasuredFrequencyResponse(
                np.array([0.0, 1e9]),
                np.array([1.0, 1.0j]),
                "dimensionless",
                "",
            ),
            "reference_plane",
        ),
        (
            lambda: fir_transfer_from_measured_response(
                MeasuredFrequencyResponse(
                    np.array([0.0, 49e9]),
                    np.ones(2, dtype=complex),
                    "dimensionless",
                    "PD response",
                ),
                sample_rate_hz=100e9,
                fft_size=128,
                fir_length_samples=1,
            ),
            "cover the Nyquist frequency",
        ),
        (
            lambda: fir_transfer_from_measured_response(
                MeasuredFrequencyResponse(
                    np.array([0.0, 50e9]),
                    np.array([1.0, np.exp(-1j * np.pi / 4)]),
                    "dimensionless",
                    "PD response",
                ),
                sample_rate_hz=100e9,
                fft_size=128,
                fir_length_samples=1,
            ),
            "Nyquist response must be real",
        ),
    ],
)
def test_measured_response_rejects_invalid_physical_inputs(
    call: Callable[[], object],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        call()


def _flat_response(
    sample_rate_hz: float,
    gain: float,
    response_kind: str,
) -> MeasuredFrequencyResponse:
    return MeasuredFrequencyResponse(
        frequency_hz=np.array([0.0, sample_rate_hz / 2]),
        response_complex=np.array([gain, gain], dtype=complex),
        response_kind=response_kind,  # type: ignore[arg-type]
        reference_plane="synthetic flat response",
    )


def _sample_transfer_response(
    transfer: DiscreteTransferFunction,
    fft_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    frequency_hz = np.fft.rfftfreq(fft_size, d=1 / transfer.sample_rate_hz)
    z_inverse = np.exp(-1j * 2 * np.pi * frequency_hz / transfer.sample_rate_hz)
    numerator = np.polynomial.polynomial.polyval(z_inverse, transfer.numerator)
    denominator = np.polynomial.polynomial.polyval(z_inverse, transfer.denominator)
    response_complex = numerator / denominator
    response_complex[-1] = response_complex[-1].real
    return frequency_hz, response_complex
