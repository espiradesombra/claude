"""Policy Engine — decide si un KOP cuesta demasiado o ya existe."""

from __future__ import annotations

from dataclasses import dataclass
from .kop import KOPId
from .runtime_context import KOPCost, RuntimeContext


@dataclass(frozen=True)
class PolicyDecision:
    allowed: bool
    reason: str = ""


class PolicyEngine:
    """Política Balanced por defecto — extensible sin tocar el núcleo."""

    def __init__(self, max_cpu: float = 1000.0, allow_network: bool = False) -> None:
        self.max_cpu = max_cpu
        self.allow_network = allow_network
        self.name = "Balanced"

    def allow(self, kop_id: KOPId, cost: KOPCost, ctx: RuntimeContext) -> PolicyDecision:
        if cost.cpu_units > self.max_cpu:
            return PolicyDecision(False, f"CPU {cost.cpu_units} > max {self.max_cpu}")
        if cost.network_units > 0 and not self.allow_network:
            return PolicyDecision(False, "network disabled in local policy")
        return PolicyDecision(True, "ok")

    def should_promote(self, hits: int, trust: float) -> bool:
        return hits >= 3 and trust >= 0.8

    def should_expire(self, hits: int, tier_cold: bool) -> bool:
        return hits == 0 and tier_cold