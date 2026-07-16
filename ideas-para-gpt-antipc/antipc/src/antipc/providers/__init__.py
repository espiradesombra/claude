"""Compute Providers — Wave, Mechanical (R001–R003 gptcomputing)."""

from .compute_provider import ComputeProvider, ProviderResult
from .mechanical_provider import MechanicalProvider
from .wave_provider import WaveProvider

__all__ = [
    "ComputeProvider",
    "ProviderResult",
    "WaveProvider",
    "MechanicalProvider",
]