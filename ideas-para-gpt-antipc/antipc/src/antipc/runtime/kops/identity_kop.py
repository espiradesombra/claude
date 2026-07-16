"""IdentityKOP — crea KnowledgeObject y lo registra."""

from __future__ import annotations

from ..kop import KOPId
from ..knowledge_object import KnowledgeObject
from ..runtime_context import KOPCost, KOPResult, KOPStatus, RuntimeContext
from ..signature import make_signature


class IdentityKOP:
    kop_id = KOPId.IDENTITY

    def knowledge_signature(self, payload: bytes | None) -> str | None:
        if not payload:
            return None
        label = payload.decode("utf-8", errors="replace")[:32]
        return make_signature("KOP_IDENTITY", [payload.hex(), label])

    def cost(self) -> KOPCost:
        return KOPCost(cpu_units=2.0, ram_bytes=512, trust_gain=0.1)

    def requires(self) -> tuple[KOPId, ...]:
        return ()

    def produces(self) -> tuple[KOPId, ...]:
        return (KOPId.IDENTITY, KOPId.DNA, KOPId.SIGNATURE)

    def execute(self, ctx: RuntimeContext, payload: bytes | None = None) -> KOPResult:
        if not payload:
            return KOPResult(KOPStatus.ERROR, int(self.kop_id), "payload vacio")

        sig = self.knowledge_signature(payload)
        assert sig is not None
        label = payload.decode("utf-8", errors="replace")[:32]
        existing = ctx.knowledge.lookup(sig)
        if existing is not None:
            return KOPResult(
                KOPStatus.SKIPPED,
                int(self.kop_id),
                f"reuse {existing.ref_id}",
                data=existing,
            )

        ko = KnowledgeObject.create(payload, label=label)
        ko.ref = type(ko.ref)(
            ref_id=ko.ref.ref_id,
            signature=sig,
            metadata=ko.ref.metadata,
            payload=ko.ref.payload,
            parents=ko.ref.parents,
        )
        ctx.registry.register(ko)
        ctx.knowledge.publish(ko.ref.signature, ko.ref)
        blob = ctx.knowledge.publish_with_kop_blob(
            ko.ref.signature, ko.ref, producer=ko.producer, trust=ko.trust
        )
        ctx.ledger.append(
            action="CREATE",
            knowledge_id=ko.id(),
            producer=ko.producer,
            trust=ko.trust,
            kop_id=self.kop_id,
            detail=f"blob={len(blob)}B",
        )
        from ..event_bus import Event, EventType

        ctx.events.emit(
            Event(EventType.REFERENCE_PUBLISHED, ko.id(), ko.ref.signature)
        )
        return KOPResult(
            KOPStatus.OK,
            int(self.kop_id),
            f"KnowledgeObject {ko.id()}",
            data=ko,
            cost=self.cost(),
        )