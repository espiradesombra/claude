"""ComputeProvider — abstracción gptcomputing (CPU/GPU/Wave/Mechanical)."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

from runtime.plugin import Capability, Cost


@dataclass
class ProviderResult:
    payload: bytes
    metadata: dict[str, str | int | float]


class ComputeProvider(ABC):
    @abstractmethod
    def capability(self) -> Capability: ...

    @abstractmethod
    def estimate(self) -> Cost: ...

    @abstractmethod
    def execute(self) -> ProviderResult: ...