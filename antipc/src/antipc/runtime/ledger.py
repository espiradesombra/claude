"""Knowledge Ledger — append-only, entradas como extension KOP 005 / K3MK."""

from __future__ import annotations

import struct
import time
from dataclasses import dataclass, field

from .kop import KOPId, K3MicroKernel


@dataclass(frozen=True)
class LedgerEntry:
    tick: int
    action: str
    knowledge_id: str
    producer: str
    parents: str
    trust: float
    policy: str
    node: str
    kop_id: int | None = None
    detail: str = ""

    def to_kop_payload(self) -> bytes:
        line = (
            f"tick={self.tick}|action={self.action}|kid={self.knowledge_id}|"
            f"prod={self.producer}|par={self.parents}|trust={self.trust:.2f}|"
            f"pol={self.policy}|node={self.node}|kop={self.kop_id}|{self.detail}"
        )
        return line.encode("utf-8")[:256]


@dataclass
class KnowledgeLedger:
    """Libro mayor interno — NO blockchain."""

    policy_name: str = "Balanced"
    node: str = "local"
    _tick: int = 0
    _entries: list[LedgerEntry] = field(default_factory=list)
    _blobs: list[bytes] = field(default_factory=list)

    def append(
        self,
        action: str,
        knowledge_id: str = "",
        producer: str = "runtime",
        parents: str = "0",
        trust: float = 1.0,
        kop_id: KOPId | None = None,
        detail: str = "",
    ) -> LedgerEntry:
        self._tick += 1
        entry = LedgerEntry(
            tick=self._tick,
            action=action,
            knowledge_id=knowledge_id,
            producer=producer,
            parents=parents,
            trust=trust,
            policy=self.policy_name,
            node=self.node,
            kop_id=int(kop_id) if kop_id is not None else None,
            detail=detail,
        )
        self._entries.append(entry)
        blob = self._entry_to_k3mk(entry)
        self._blobs.append(blob)
        return entry

    def _entry_to_k3mk(self, entry: LedgerEntry) -> bytes:
        mk = K3MicroKernel()
        mk.set_kop(KOPId.HISTORY, entry.to_kop_payload())
        mk.set_kop(KOPId.IDENTITY, entry.knowledge_id.encode("utf-8")[:64] or b"ledger")
        mk.set_kop(KOPId.DNA, f"producer={entry.producer}|parents={entry.parents}".encode())
        mk.set_kop(KOPId.TRUST, struct.pack("<d", entry.trust))
        return mk.pack()

    def count(self) -> int:
        return len(self._entries)

    def last(self) -> LedgerEntry | None:
        return self._entries[-1] if self._entries else None

    def export_blobs(self) -> list[bytes]:
        return list(self._blobs)