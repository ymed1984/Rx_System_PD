import pytest

from oma_ber.pd_analysis import compare_photodiodes_at_optical_levels, required_oma_by_photodiode
from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver


def example_rx() -> Receiver:
    return Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
    )


def pd_variants() -> dict[str, Photodiode]:
    return {
        "low_R": Photodiode(responsivity_a_per_w=0.5, dark_current_a=1e-9),
        "high_R": Photodiode(responsivity_a_per_w=1.0, dark_current_a=1e-9),
    }


def test_compare_photodiodes_uses_identical_rx_input_levels() -> None:
    results = compare_photodiodes_at_optical_levels(
        p0_w=1e-6,
        p1_w=4e-6,
        photodiodes=pd_variants(),
        rx=example_rx(),
    )

    assert results["low_R"]["p0_w"] == pytest.approx(results["high_R"]["p0_w"])
    assert results["low_R"]["p1_w"] == pytest.approx(results["high_R"]["p1_w"])
    assert results["high_R"]["q_rx"] > results["low_R"]["q_rx"]


def test_required_oma_comparison_reports_pd_input_requirement() -> None:
    required = required_oma_by_photodiode(
        target_ber=1e-6,
        er_db=6.0,
        photodiodes=pd_variants(),
        rx=example_rx(),
        oma_dbm_min=-40.0,
        oma_dbm_max=0.0,
    )

    assert required["high_R"] < required["low_R"]


def test_pd_comparison_rejects_empty_variants() -> None:
    with pytest.raises(ValueError, match="photodiodes must contain at least one"):
        compare_photodiodes_at_optical_levels(1e-6, 2e-6, {}, example_rx())
