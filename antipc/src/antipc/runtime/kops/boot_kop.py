"""BootKOP — primer KOP ejecutable (health check)."""

from __future__ import annotations

from ..kop import KOPId
from ..runtime_context import KOPCost, KOPResult, KOPStatus, RuntimeContext


class BootKOP:
    kop_id = KOPId.IDENTITY

    def cost(self) -> KOPCost:
        return KOPCost(cpu_units=0.1)

    def requires(self) -> tuple[KOPId, ...]:
        return ()

    def produces(self) -> tuple[KOPId, ...]:
        return (KOPId.IDENTITY,)

    def execute(self, ctx: RuntimeContext, payload: bytes | None = None) -> KOPResult:
        return KOPResult(
            KOPStatus.OK,
            int(self.kop_id),
            "BootKOP health check",
            data={"registry": ctx.registry.count()},
        )