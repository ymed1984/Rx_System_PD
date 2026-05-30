"""Example: sweep a simplified static photodiode saturation model."""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oma_ber.saturation import (
    compressed_responsivity_a_per_w,
    saturated_photocurrent_a,
    saturated_photocurrent_exponential_a,
    saturated_photocurrent_rational_a,
    saturated_photocurrent_soft_clip_a,
)


def main() -> None:
    responsivity_a_per_w = 0.8
    saturation_power_w = 1e-3
    optical_power_w_values = np.array([1e-6, 1e-5, 1e-4, 1e-3, 3e-3])

    responsivity_values = compressed_responsivity_a_per_w(
        optical_power_w=optical_power_w_values,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
    )
    tanh_current_values = saturated_photocurrent_a(
        optical_power_w=optical_power_w_values,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
    )
    rational_current_values = saturated_photocurrent_rational_a(
        optical_power_w=optical_power_w_values,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
    )
    soft_clip_current_values = saturated_photocurrent_soft_clip_a(
        optical_power_w=optical_power_w_values,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
        compression_order=2.0,
    )
    exponential_current_values = saturated_photocurrent_exponential_a(
        optical_power_w=optical_power_w_values,
        responsivity_a_per_w=responsivity_a_per_w,
        saturation_power_w=saturation_power_w,
    )

    print("Example simplified static photodiode saturation sweep")
    print(f"Small-signal responsivity: {responsivity_a_per_w:.3f} A/W")
    print(f"Saturation power scale: {saturation_power_w:.3e} W")
    print(
        "Optical power [W], tanh default current [A], rational current [A], "
        "soft-clip current [A], exponential current [A], tanh R_eff [A/W]"
    )
    for optical_power_w, tanh_current_a, rational_current_a, soft_clip_current_a, exponential_current_a, responsivity in zip(
        optical_power_w_values,
        tanh_current_values,
        rational_current_values,
        soft_clip_current_values,
        exponential_current_values,
        responsivity_values,
        strict=True,
    ):
        print(
            f"{optical_power_w:.3e}, {tanh_current_a:.3e}, {rational_current_a:.3e}, "
            f"{soft_clip_current_a:.3e}, {exponential_current_a:.3e}, {responsivity:.3e}"
        )


if __name__ == "__main__":
    main()
