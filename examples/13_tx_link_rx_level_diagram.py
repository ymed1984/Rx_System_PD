"""Example: per-wavelength WDM TX -> fiber link -> RX level diagrams."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber import (
    ExternalLaser,
    OOKModulator,
    OpticalFiberLink,
    Photodiode,
    Receiver,
    WdmChannel,
    analyze_wdm_link,
    compare_photodiodes_at_optical_levels,
)


def main() -> None:
    channels = [
        WdmChannel(
            name="lambda_1310",
            wavelength_nm=1310.0,
            laser=ExternalLaser(output_power_dbm=-4.0),
            modulator=OOKModulator(
                insertion_loss_db=3.0,
                er_db=6.0,
                specified_oma_dbm=-9.5,
            ),
            laser_to_modulator_coupling_loss_db=1.0,
            awg_split_loss_db=2.0,
            tx_fiber_coupling_loss_db=1.0,
            rx_fiber_coupling_loss_db=1.0,
        ),
        WdmChannel(
            name="lambda_1330",
            wavelength_nm=1330.0,
            laser=ExternalLaser(output_power_dbm=-3.5),
            modulator=OOKModulator(
                insertion_loss_db=3.2,
                er_db=6.0,
                specified_oma_dbm=-9.3,
            ),
            laser_to_modulator_coupling_loss_db=1.1,
            awg_split_loss_db=2.3,
            tx_fiber_coupling_loss_db=1.0,
            rx_fiber_coupling_loss_db=1.1,
            fiber_attenuation_db_per_km=0.34,
        ),
    ]
    fiber_link = OpticalFiberLink(
        length_km=10.0,
        attenuation_db_per_km=0.35,
        additional_loss_db=0.2,
    )
    pd = Photodiode(
        responsivity_a_per_w=0.8,
        dark_current_a=1e-9,
        dark_current_fano_factor=1.0,
    )
    rx = Receiver(
        noise_bandwidth_hz=25e9,
        input_current_noise_density_a_per_sqrt_hz=10e-12,
        rin_db_per_hz=-150,
    )

    results = analyze_wdm_link(
        channels=channels,
        fiber_link=fiber_link,
        margin_db=0.5,
        pd=pd,
        rx=rx,
    )

    print("WDM per-channel level diagrams")
    print(
        f"Common fiber: {fiber_link.length_km:.1f} km, "
        f"default loss={fiber_link.propagation_loss_db:.2f} dB",
    )
    for name, result in results.items():
        channel = result.channel
        print(
            f"\n{name}: wavelength={channel.wavelength_nm:.1f} nm, "
            f"fiber loss={result.fiber_propagation_loss_db:.2f} dB",
        )
        print(
            f"Calculated MOD OMA={result.modulator_output.calculated_oma_dbm:.2f} dBm, "
            f"calculated-specified={result.modulator_output.oma_error_db:+.2f} dB",
        )
        print("Section           Stage                    Inc.loss  Cum.loss  Optical levels")
        for point in result.level_diagram:
            if point.is_modulated:
                levels = (
                    f"P0={point.p0_w * 1e3:.5f} mW, "
                    f"P1={point.p1_w * 1e3:.5f} mW, "
                    f"OMA={point.oma_dbm:.2f} dBm"
                )
            else:
                levels = f"CW={point.cw_power_dbm:.2f} dBm"
            print(
                f"{point.section:17s} {point.name:24s} "
                f"{point.incremental_loss_db:8.2f} {point.cumulative_loss_db:9.2f}  {levels}",
            )
        receiver = result.receiver_result
        print(
            f"PD input: OMA={result.pd_input.oma_dbm:.2f} dBm, "
            f"Q={receiver['q_rx']:.3f}, BER={receiver['ber']:.3e}",
        )

    reference = results["lambda_1310"]
    pd_comparison = compare_photodiodes_at_optical_levels(
        p0_w=reference.pd_input.p0_w,
        p1_w=reference.pd_input.p1_w,
        photodiodes={
            "baseline": pd,
            "lower responsivity": Photodiode(
                responsivity_a_per_w=0.6,
                dark_current_a=1e-9,
            ),
        },
        rx=rx,
    )
    print("\nRX-only PD comparison for lambda_1310 at identical P0/P1")
    for name, pd_result in pd_comparison.items():
        print(f"{name:20s} Q={pd_result['q_rx']:.3f}, BER={pd_result['ber']:.3e}")


if __name__ == "__main__":
    main()
