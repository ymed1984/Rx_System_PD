"""Per-wavelength WDM TX, fiber-link, and RX level-diagram analysis."""

from collections.abc import Sequence
from dataclasses import dataclass
from math import isfinite
from typing import Any

from oma_ber.link_budget import apply_optical_loss
from oma_ber.modulation import NRZLevels
from oma_ber.photodiode import Photodiode
from oma_ber.receiver import Receiver
from oma_ber.sweep import calculate_ber_from_optical_levels
from oma_ber.transmitter import ExternalLaser, ModulatorOutput, OOKModulator, modulator_output_from_laser
from oma_ber.units import db_to_linear, watt_to_dbm


def _validate_loss_db(value: float, name: str) -> None:
    if not isfinite(value) or value < 0:
        msg = f"{name} must be finite and non-negative."
        raise ValueError(msg)


@dataclass(frozen=True)
class OpticalFiberLink:
    """Optical-fiber propagation loss parameters.

    length_km is communication distance in km, attenuation_db_per_km is fiber
    power attenuation in dB/km, and additional_loss_db is an optional fixed
    passive loss in dB assigned to the fiber-link section.
    """

    length_km: float
    attenuation_db_per_km: float
    additional_loss_db: float = 0.0

    def __post_init__(self) -> None:
        if not isfinite(self.length_km) or self.length_km < 0:
            msg = "length_km must be finite and non-negative."
            raise ValueError(msg)
        _validate_loss_db(self.attenuation_db_per_km, "attenuation_db_per_km")
        _validate_loss_db(self.additional_loss_db, "additional_loss_db")

    @property
    def propagation_loss_db(self) -> float:
        """Return total optical-fiber link loss in dB."""
        return self.length_km * self.attenuation_db_per_km + self.additional_loss_db


@dataclass(frozen=True)
class WdmChannel:
    """One wavelength channel and its TX/RX coupling-loss parameters.

    The modeled chain is External Laser -> laser-to-MOD coupling -> MOD -> AWG
    split loss -> TX fiber coupling -> optical-fiber link -> RX fiber coupling.
    All losses are optical power losses in dB. Modulator insertion loss remains
    referenced to the optical 1 level.
    """

    name: str
    wavelength_nm: float
    laser: ExternalLaser
    modulator: OOKModulator
    laser_to_modulator_coupling_loss_db: float
    awg_split_loss_db: float
    tx_fiber_coupling_loss_db: float
    rx_fiber_coupling_loss_db: float
    fiber_attenuation_db_per_km: float | None = None

    def __post_init__(self) -> None:
        if not self.name.strip():
            msg = "name must not be empty."
            raise ValueError(msg)
        if not isfinite(self.wavelength_nm) or self.wavelength_nm <= 0:
            msg = "wavelength_nm must be finite and positive."
            raise ValueError(msg)
        _validate_loss_db(
            self.laser_to_modulator_coupling_loss_db,
            "laser_to_modulator_coupling_loss_db",
        )
        _validate_loss_db(self.awg_split_loss_db, "awg_split_loss_db")
        _validate_loss_db(self.tx_fiber_coupling_loss_db, "tx_fiber_coupling_loss_db")
        _validate_loss_db(self.rx_fiber_coupling_loss_db, "rx_fiber_coupling_loss_db")
        if self.fiber_attenuation_db_per_km is not None:
            _validate_loss_db(
                self.fiber_attenuation_db_per_km,
                "fiber_attenuation_db_per_km",
            )


@dataclass(frozen=True)
class WdmLevelPoint:
    """CW or NRZ/OOK optical levels at one WDM-channel reference plane.

    cumulative_loss_db is referenced to the External Laser and follows the
    optical-1 path after modulation. CW power fields are populated before the
    modulator; P0/P1/Pavg/OMA fields are populated from MOD output onward.
    """

    name: str
    section: str
    wavelength_nm: float
    incremental_loss_db: float
    cumulative_loss_db: float
    cw_power_w: float | None = None
    cw_power_dbm: float | None = None
    p0_w: float | None = None
    p1_w: float | None = None
    pavg_w: float | None = None
    oma_w: float | None = None
    oma_dbm: float | None = None
    er_db: float | None = None

    @property
    def is_modulated(self) -> bool:
        """Return True when this point contains NRZ/OOK P0/P1 levels."""
        return self.p0_w is not None


@dataclass(frozen=True)
class WdmChannelResult:
    """Level diagram and receiver result for one WDM wavelength channel."""

    channel: WdmChannel
    fiber_link: OpticalFiberLink
    fiber_propagation_loss_db: float
    modulator_output: ModulatorOutput
    level_diagram: tuple[WdmLevelPoint, ...]
    receiver_result: dict[str, Any]

    @property
    def pd_input(self) -> WdmLevelPoint:
        """Return the final PD-input level point."""
        return self.level_diagram[-1]


def _cw_point(
    name: str,
    section: str,
    wavelength_nm: float,
    incremental_loss_db: float,
    cumulative_loss_db: float,
    cw_power_w: float,
) -> WdmLevelPoint:
    return WdmLevelPoint(
        name=name,
        section=section,
        wavelength_nm=wavelength_nm,
        incremental_loss_db=incremental_loss_db,
        cumulative_loss_db=cumulative_loss_db,
        cw_power_w=cw_power_w,
        cw_power_dbm=watt_to_dbm(cw_power_w),
    )


def _modulated_point(
    name: str,
    section: str,
    wavelength_nm: float,
    incremental_loss_db: float,
    cumulative_loss_db: float,
    levels: NRZLevels,
) -> WdmLevelPoint:
    return WdmLevelPoint(
        name=name,
        section=section,
        wavelength_nm=wavelength_nm,
        incremental_loss_db=incremental_loss_db,
        cumulative_loss_db=cumulative_loss_db,
        p0_w=levels.p0_w,
        p1_w=levels.p1_w,
        pavg_w=levels.pavg_w,
        oma_w=levels.oma_w,
        oma_dbm=watt_to_dbm(levels.oma_w),
        er_db=levels.er_db,
    )


def analyze_wdm_channel(
    channel: WdmChannel,
    fiber_link: OpticalFiberLink,
    pd: Photodiode,
    rx: Receiver,
    margin_db: float = 0.0,
    optimize_threshold: bool = True,
) -> WdmChannelResult:
    """Analyze one WDM channel from External Laser through the RX PD input."""
    _validate_loss_db(margin_db, "margin_db")
    wavelength_nm = channel.wavelength_nm
    cumulative_loss_db = 0.0
    points = [
        _cw_point(
            name="External Laser",
            section="TX",
            wavelength_nm=wavelength_nm,
            incremental_loss_db=0.0,
            cumulative_loss_db=0.0,
            cw_power_w=channel.laser.output_power_w,
        ),
    ]

    coupling_loss_db = channel.laser_to_modulator_coupling_loss_db
    cumulative_loss_db += coupling_loss_db
    modulator_input_power_w = channel.laser.output_power_w / db_to_linear(coupling_loss_db)
    points.append(
        _cw_point(
            name="Laser-to-MOD coupling",
            section="TX",
            wavelength_nm=wavelength_nm,
            incremental_loss_db=coupling_loss_db,
            cumulative_loss_db=cumulative_loss_db,
            cw_power_w=modulator_input_power_w,
        ),
    )

    effective_laser = ExternalLaser(output_power_dbm=watt_to_dbm(modulator_input_power_w))
    modulator_output = modulator_output_from_laser(effective_laser, channel.modulator)
    cumulative_loss_db += channel.modulator.insertion_loss_db
    current_levels = modulator_output.levels
    points.append(
        _modulated_point(
            name="MOD output",
            section="TX",
            wavelength_nm=wavelength_nm,
            incremental_loss_db=channel.modulator.insertion_loss_db,
            cumulative_loss_db=cumulative_loss_db,
            levels=current_levels,
        ),
    )

    fiber_attenuation_db_per_km = (
        channel.fiber_attenuation_db_per_km
        if channel.fiber_attenuation_db_per_km is not None
        else fiber_link.attenuation_db_per_km
    )
    fiber_propagation_loss_db = (
        fiber_link.length_km * fiber_attenuation_db_per_km + fiber_link.additional_loss_db
    )
    link_stages = (
        ("AWG split loss", "TX", channel.awg_split_loss_db),
        ("TX fiber coupling", "TX", channel.tx_fiber_coupling_loss_db),
        (
            f"Optical fiber ({fiber_link.length_km:g} km)",
            "OpticalFiberLink",
            fiber_propagation_loss_db,
        ),
        ("RX fiber coupling", "RX", channel.rx_fiber_coupling_loss_db),
    )
    for name, section, loss_db in link_stages:
        current_levels = apply_optical_loss(current_levels, loss_db)
        cumulative_loss_db += loss_db
        points.append(
            _modulated_point(
                name=name,
                section=section,
                wavelength_nm=wavelength_nm,
                incremental_loss_db=loss_db,
                cumulative_loss_db=cumulative_loss_db,
                levels=current_levels,
            ),
        )

    if margin_db > 0:
        current_levels = apply_optical_loss(current_levels, margin_db)
        cumulative_loss_db += margin_db
        points.append(
            _modulated_point(
                name="Design margin",
                section="RX",
                wavelength_nm=wavelength_nm,
                incremental_loss_db=margin_db,
                cumulative_loss_db=cumulative_loss_db,
                levels=current_levels,
            ),
        )
    points.append(
        _modulated_point(
            name="PD input",
            section="RX",
            wavelength_nm=wavelength_nm,
            incremental_loss_db=0.0,
            cumulative_loss_db=cumulative_loss_db,
            levels=current_levels,
        ),
    )

    receiver_result = calculate_ber_from_optical_levels(
        p0_w=current_levels.p0_w,
        p1_w=current_levels.p1_w,
        pd=pd,
        rx=rx,
        optimize_threshold=optimize_threshold,
    )
    return WdmChannelResult(
        channel=channel,
        fiber_link=fiber_link,
        fiber_propagation_loss_db=fiber_propagation_loss_db,
        modulator_output=modulator_output,
        level_diagram=tuple(points),
        receiver_result=receiver_result,
    )


def analyze_wdm_link(
    channels: Sequence[WdmChannel],
    fiber_link: OpticalFiberLink,
    pd: Photodiode,
    rx: Receiver,
    margin_db: float = 0.0,
    optimize_threshold: bool = True,
) -> dict[str, WdmChannelResult]:
    """Analyze multiple independent WDM channels over a common fiber link.

    AWG crosstalk and inter-channel optical noise are not included. The same PD
    and receiver parameters are used for every channel.
    """
    if not channels:
        msg = "channels must contain at least one WdmChannel."
        raise ValueError(msg)
    if any(not isinstance(channel, WdmChannel) for channel in channels):
        msg = "channels must contain only WdmChannel values."
        raise ValueError(msg)
    names = [channel.name for channel in channels]
    if len(set(names)) != len(names):
        msg = "WDM channel names must be unique."
        raise ValueError(msg)
    return {
        channel.name: analyze_wdm_channel(
            channel=channel,
            fiber_link=fiber_link,
            pd=pd,
            rx=rx,
            margin_db=margin_db,
            optimize_threshold=optimize_threshold,
        )
        for channel in channels
    }
