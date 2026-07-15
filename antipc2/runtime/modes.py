"""Runtime modes: Local, Battery, Network, Research, Benchmark."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto


class RuntimeModeName(Enum):
    LOCAL = auto()
    BATTERY = auto()
    NETWORK = auto()
    RESEARCH = auto()
    BENCHMARK = auto()


@dataclass
class CachePolicy:
    reuse_knowledge: bool = True
    prefer_hot: bool = True


@dataclass
class NetworkPolicy:
    enabled: bool = False
    require_permission: bool = True
    allow_discovery: bool = False


@dataclass
class EnergyPolicy:
    batch_queries: bool = False
    defer_compute: bool = False


class RuntimeMode(ABC):
    @property
    @abstractmethod
    def name(self) -> RuntimeModeName: ...

    @abstractmethod
    def cache_policy(self) -> CachePolicy: ...

    @abstractmethod
    def network_policy(self) -> NetworkPolicy: ...

    @abstractmethod
    def energy_policy(self) -> EnergyPolicy: ...

    def priority_boost(self, cacheable: bool, reusable_hit: bool) -> float:
        return 0.0


class LocalMode(RuntimeMode):
    @property
    def name(self) -> RuntimeModeName:
        return RuntimeModeName.LOCAL

    def cache_policy(self) -> CachePolicy:
        return CachePolicy(reuse_knowledge=True, prefer_hot=True)

    def network_policy(self) -> NetworkPolicy:
        return NetworkPolicy(enabled=False)

    def energy_policy(self) -> EnergyPolicy:
        return EnergyPolicy()


class BatteryMode(RuntimeMode):
    """Planner: agrupa, reutiliza, no despierta núcleos innecesarios."""

    @property
    def name(self) -> RuntimeModeName:
        return RuntimeModeName.BATTERY

    def cache_policy(self) -> CachePolicy:
        return CachePolicy(reuse_knowledge=True, prefer_hot=True)

    def network_policy(self) -> NetworkPolicy:
        return NetworkPolicy(enabled=False)

    def energy_policy(self) -> EnergyPolicy:
        return EnergyPolicy(batch_queries=True, defer_compute=True)

    def priority_boost(self, cacheable: bool, reusable_hit: bool) -> float:
        return 2.0 if reusable_hit else (0.5 if cacheable else 0.0)


class NetworkMode(RuntimeMode):
    @property
    def name(self) -> RuntimeModeName:
        return RuntimeModeName.NETWORK

    def cache_policy(self) -> CachePolicy:
        return CachePolicy(reuse_knowledge=True)

    def network_policy(self) -> NetworkPolicy:
        return NetworkPolicy(enabled=True, require_permission=True, allow_discovery=True)

    def energy_policy(self) -> EnergyPolicy:
        return EnergyPolicy()


class ResearchMode(RuntimeMode):
    """No reutilizar — ejecutar siempre para validar hipótesis."""

    @property
    def name(self) -> RuntimeModeName:
        return RuntimeModeName.RESEARCH

    def cache_policy(self) -> CachePolicy:
        return CachePolicy(reuse_knowledge=False)

    def network_policy(self) -> NetworkPolicy:
        return NetworkPolicy(enabled=False)

    def energy_policy(self) -> EnergyPolicy:
        return EnergyPolicy()


class BenchmarkMode(RuntimeMode):
    @property
    def name(self) -> RuntimeModeName:
        return RuntimeModeName.BENCHMARK

    def cache_policy(self) -> CachePolicy:
        return CachePolicy(reuse_knowledge=True)

    def network_policy(self) -> NetworkPolicy:
        return NetworkPolicy(enabled=True, require_permission=False)

    def energy_policy(self) -> EnergyPolicy:
        return EnergyPolicy()