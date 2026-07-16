"""KOP Executor — ISA del microkernel: Load → Execute → Return."""

from __future__ import annotations

from typing import Any, Protocol

from .kop import KOPId, KOP_DEFINITIONS
from .kop_registry import KOPRegistry
from .runtime_context import KOPCost, KOPResult, KOPStatus, RuntimeContext


class IKOPHandler(Protocol):
    kop_id: KOPId

    def cost(self) -> KOPCost: ...
    def requires(self) -> tuple[KOPId, ...]: ...
    def produces(self) -> tuple[KOPId, ...]: ...
    def execute(self, ctx: RuntimeContext, payload: bytes | None = None) -> KOPResult: ...


class KOPExecutor:
    def __init__(self) -> None:
        self._handlers: dict[KOPId, IKOPHandler] = {}

    def register_handler(self, handler: IKOPHandler) -> None:
        self._handlers[handler.kop_id] = handler
        KOPRegistry.register(handler.kop_id)

    def get_handler(self, kop_id: KOPId) -> IKOPHandler | None:
        return self._handlers.get(kop_id)

    def execute(
        self,
        ctx: RuntimeContext,
        kop_id: KOPId,
        payload: bytes | None = None,
    ) -> KOPResult:
        handler = self._handlers.get(kop_id)
        if handler is None:
            return KOPResult(KOPStatus.ERROR, int(kop_id), f"KOP {kop_id} no registrado")

        cost = handler.cost()
        decision = ctx.policy.allow(kop_id, cost, ctx)
        if not decision.allowed:
            ctx.metrics.record_kop_denied(kop_id)
            return KOPResult(KOPStatus.DENIED, int(kop_id), decision.reason, cost=cost)

        if ctx.metrics.already_computed(handler, payload, ctx):
            ctx.metrics.record_reuse(kop_id)
            return KOPResult(KOPStatus.SKIPPED, int(kop_id), "reuse", cost=KOPCost())

        result = handler.execute(ctx, payload)
        ctx.metrics.record_kop_executed(kop_id, result.cost)
        ctx.ledger.append(
            action="EXECUTE_KOP",
            kop_id=kop_id,
            detail=result.message,
            trust=1.0,
        )
        return result

    def describe(self, kop_id: KOPId) -> dict[str, Any]:
        defn = KOP_DEFINITIONS.get(kop_id, {})
        handler = self._handlers.get(kop_id)
        return {
            "id": int(kop_id),
            "name": defn.get("name", "?"),
            "registered": handler is not None,
            "requires": [int(r) for r in handler.requires()] if handler else [],
            "produces": [int(p) for p in handler.produces()] if handler else [],
        }