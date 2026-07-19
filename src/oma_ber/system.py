"""Integrated external-laser, modulator, optical-link, and receiver analysis."""

from dataclasses import dataclass
from typing import Any, Sequence

from oma_ber.link_budget import OpticalLevelPoint, OpticalLoss, build_level_diagram
from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.sweep import calculate_ber_from_optical_levels
from oma_ber.transmitter import ExternalLaser, ModulatorOutput, OOKModulator, modulator_output_from_laser


@dataclass(frozen=True)
class TxLinkRxResult:
    """Integrated TX/link/RX result with explicit optical reference planes."""

    laser: ExternalLaser
    modulator: OOKModulator
    modulator_output: ModulatorOutput
    level_diagram: tuple[OpticalLevelPoint, ...]
    receiver_result: dict[str, Any]

    @property
    def pd_input(self) -> OpticalLevelPoint:
        """Return the final PD-input optical level point."""
        return self.level_diagram[-1]


def analyze_laser_to_receiver(
    laser: ExternalLaser,
    modulator: OOKModulator,
    losses: Sequence[OpticalLoss],
    pd: Photodiode,
    rx: Receiver,
    margin_db: float = 0.0,
    optimize_threshold: bool = True,
) -> TxLinkRxResult:
    """Analyze External Laser -> OOK modulator -> passive link -> PD/TIA.

    Modulator insertion loss is referenced to the optical 1 level. Link losses
    and margin are level-independent passive optical power losses in dB.
    """
    modulator_output = modulator_output_from_laser(laser, modulator)
    diagram = build_level_diagram(
        initial_levels=modulator_output.levels,
        losses=losses,
        margin_db=margin_db,
    )
    final_point = diagram[-1]
    if final_point.name != "PD input":
        diagram.append(
            OpticalLevelPoint(
                name="PD input",
                incremental_loss_db=0.0,
                cumulative_loss_db=final_point.cumulative_loss_db,
                p0_w=final_point.p0_w,
                p1_w=final_point.p1_w,
                pavg_w=final_point.pavg_w,
                oma_w=final_point.oma_w,
                oma_dbm=final_point.oma_dbm,
                er_linear=final_point.er_linear,
                er_db=final_point.er_db,
            ),
        )
    pd_input = diagram[-1]
    receiver_result = calculate_ber_from_optical_levels(
        p0_w=pd_input.p0_w,
        p1_w=pd_input.p1_w,
        pd=pd,
        rx=rx,
        optimize_threshold=optimize_threshold,
    )
    return TxLinkRxResult(
        laser=laser,
        modulator=modulator,
        modulator_output=modulator_output,
        level_diagram=tuple(diagram),
        receiver_result=receiver_result,
    )
