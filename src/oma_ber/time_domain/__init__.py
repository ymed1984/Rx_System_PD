"""Physically separated time-domain NRZ/OOK receiver analysis.

This package keeps optical power in W, photodiode current in A, and TIA output
in V as distinct data types. It covers deterministic single-channel IM/DD
receiver waveforms plus PSD-based Gaussian statistical eye/BER analysis.
Jitter, CDR, and complex optical-field propagation remain separate extensions.
"""

from oma_ber.time_domain.eye import (
    DeterministicEyeAnalysis,
    EyeLevelMetrics,
    analyze_deterministic_eye,
    sample_eye_at_phase,
)
from oma_ber.time_domain.noise import (
    TiaOutputNoiseResult,
    TimeDomainNoiseModel,
    calculate_tia_output_noise,
)
from oma_ber.time_domain.patterns import prbs_bits
from oma_ber.time_domain.receiver import (
    ReceiverWaveformResult,
    simulate_pd_tia_waveform,
)
from oma_ber.time_domain.timebase import TimeGrid
from oma_ber.time_domain.statistical import (
    StatisticalEyeAnalysis,
    StatisticalEyePhaseResult,
    analyze_statistical_tia_eye,
    gaussian_mixture_ber_for_threshold,
    optimize_gaussian_mixture_threshold,
)
from oma_ber.time_domain.transfer import (
    DiscreteTransferFunction,
    equivalent_noise_bandwidth_hz,
    first_order_lowpass_transfer,
    identity_transfer,
    transfer_impulse_response,
)
from oma_ber.time_domain.waveforms import (
    OpticalPowerWaveform,
    PhotocurrentWaveform,
    VoltageWaveform,
    nrz_optical_power_waveform,
    photodetect_power_waveform,
)

__all__ = [
    "DeterministicEyeAnalysis",
    "DiscreteTransferFunction",
    "EyeLevelMetrics",
    "OpticalPowerWaveform",
    "PhotocurrentWaveform",
    "ReceiverWaveformResult",
    "StatisticalEyeAnalysis",
    "StatisticalEyePhaseResult",
    "TiaOutputNoiseResult",
    "TimeGrid",
    "TimeDomainNoiseModel",
    "VoltageWaveform",
    "analyze_deterministic_eye",
    "analyze_statistical_tia_eye",
    "calculate_tia_output_noise",
    "equivalent_noise_bandwidth_hz",
    "first_order_lowpass_transfer",
    "gaussian_mixture_ber_for_threshold",
    "identity_transfer",
    "nrz_optical_power_waveform",
    "photodetect_power_waveform",
    "prbs_bits",
    "sample_eye_at_phase",
    "simulate_pd_tia_waveform",
    "optimize_gaussian_mixture_threshold",
    "transfer_impulse_response",
]
