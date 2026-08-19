"""Knowledge Buffer K = (S, R) with HOT/WARM/COLD tiers."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto

from .reference import Reference
from .kop import build_knowledge_blob, lookup_signature_from_blob


class KnowledgeTier(Enum):
    HOT = auto()
    WARM = auto()
    COLD = auto()


@dataclass
class KnowledgeEntry:
    signature: str
    reference: Reference
    tier: KnowledgeTier = KnowledgeTier.HOT
    hits: int = 0


@dataclass
class KnowledgeBuffer:
    """Stores resolved transformations before ALU execution."""

    _by_signature: dict[str, KnowledgeEntry] = field(default_factory=dict)
    queries: int = 0
    hits: int = 0

    def lookup(self, signature: str) -> Reference | None:
        self.queries += 1
        entry = self._by_signature.get(signature)
        if entry is None:
            return None
        entry.hits += 1
        self.hits += 1
        if entry.hits > 100:
            entry.tier = KnowledgeTier.HOT
        elif entry.hits > 10:
            entry.tier = KnowledgeTier.WARM
        else:
            entry.tier = KnowledgeTier.COLD
        return entry.reference

    def publish(self, signature: str, reference: Reference) -> None:
        self._by_signature[signature] = KnowledgeEntry(signature, reference)

    def publish_with_kop_blob(
        self,
        signature: str,
        reference: Reference,
        producer: str = "antipc",
        trust: float = 1.0,
    ) -> bytes:
        """Publica y adjunta blob K3MK para lookup binario rapido."""
        payload = reference.payload or b""
        blob = build_knowledge_blob(
            reference.ref_id, signature, payload, producer=producer, trust=trust
        )
        meta = dict(reference.metadata)
        meta["kop_blob"] = blob
        meta["kop_signature"] = lookup_signature_from_blob(blob)
        updated = Reference(
            ref_id=reference.ref_id,
            signature=reference.signature,
            metadata=meta,
            payload=reference.payload,
            parents=reference.parents,
        )
        self._by_signature[signature] = KnowledgeEntry(signature, updated)
        return blob

    def lookup_by_kop_blob(self, blob: bytes) -> str | None:
        """Resuelve signature desde blob sin tocar el indice JSON."""
        return lookup_signature_from_blob(blob)