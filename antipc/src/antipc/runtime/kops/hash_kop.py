"""HashKOP — K3 hash como operación de conocimiento del microkernel."""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from hash_engine import get_backend, k3_hash_buffer

from ..event_bus import Event, EventType
from ..kop import KOPId
from ..runtime_context import KOPCost, KOPResult, KOPStatus, RuntimeContext
from ..reference import Reference
from ..signature import make_signature


class HashKOP:
    kop_id = KOPId.HASH

    def knowledge_signature(self, payload: bytes | None) -> str | None:
        if not payload:
            return None
        return make_signature("K3_HASH", [payload.hex()])

    def cost(self) -> KOPCost:
        return KOPCost(cpu_units=1.0, ram_bytes=256)

    def cost_for(self, payload: bytes) -> KOPCost:
        units = max(1.0, len(payload) / 32.0)
        return KOPCost(cpu_units=units, ram_bytes=len(payload))

    def requires(self) -> tuple[KOPId, ...]:
        return ()

    def produces(self) -> tuple[KOPId, ...]:
        return (KOPId.SIGNATURE,)

    def execute(self, ctx: RuntimeContext, payload: bytes | None = None) -> KOPResult:
        if not payload:
            return KOPResult(KOPStatus.ERROR, int(self.kop_id), "payload vacio")

        sig = self.knowledge_signature(payload)
        assert sig is not None
        existing = ctx.knowledge.lookup(sig)
        if existing is not None:
            digest_hex = existing.metadata.get("digest_hex", "")
            return KOPResult(
                KOPStatus.SKIPPED,
                int(self.kop_id),
                f"reuse {digest_hex}",
                data=existing,
                cost=KOPCost(),
            )

        digest = k3_hash_buffer(payload)
        ref = Reference.create(
            signature=sig,
            payload=digest.to_bytes(4, "little"),
            metadata={
                "plugin": "HashKOP",
                "digest_hex": f"{digest:08X}",
                "backend": get_backend(),
                "bytes": len(payload),
            },
        )
        ctx.knowledge.publish(sig, ref)
        ctx.knowledge.publish_with_kop_blob(sig, ref, producer="HashKOP", trust=1.0)
        ctx.ledger.append(
            action="HASH",
            knowledge_id=ref.ref_id,
            producer="HashKOP",
            trust=1.0,
            kop_id=self.kop_id,
            detail=f"digest={digest:08X}",
        )
        ctx.events.emit(Event(EventType.REFERENCE_PUBLISHED, ref.ref_id, sig))
        return KOPResult(
            KOPStatus.OK,
            int(self.kop_id),
            f"K3 {digest:08X}",
            data=ref,
            cost=self.cost_for(payload),
        )