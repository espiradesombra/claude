"""Metrics — P_util, knowledge hits, ALU ahorrada."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from .kop import KOPId, KOP_DEFINITIONS
from .runtime_context import KOPCost, RuntimeContext

if TYPE_CHECKING:
    from .kop_executor import IKOPHandler


@dataclass
class MetricsCollector:
    kop_executed: int = 0
    kop_denied: int = 0
    kop_reused: int = 0
    knowledge_queries: int = 0
    knowledge_hits: int = 0
    alu_saved_units: float = 0.0
    _kop_cpu: dict[int, float] = field(default_factory=dict)
    _start: float = field(default_factory=time.perf_counter)

    def record_kop_executed(self, kop_id: KOPId, cost: KOPCost) -> None:
        self.kop_executed += 1
        self._kop_cpu[int(kop_id)] = self._kop_cpu.get(int(kop_id), 0.0) + cost.cpu_units

    def record_kop_denied(self, kop_id: KOPId) -> None:
        self.kop_denied += 1

    def record_reuse(self, kop_id: KOPId) -> None:
        self.kop_reused += 1
        self.knowledge_hits += 1
        self.alu_saved_units += self._kop_cpu.get(int(kop_id), 1.0)

    def record_query(self, hit: bool) -> None:
        self.knowledge_queries += 1
        if hit:
            self.knowledge_hits += 1

    def already_computed(
        self,
        handler: IKOPHandler,
        payload: bytes | None,
        ctx: RuntimeContext,
    ) -> bool:
        if payload is None:
            return False
        sig_fn = getattr(handler, "knowledge_signature", None)
        if sig_fn is None:
            return False
        sig = sig_fn(payload)
        if sig is None:
            return False
        self.knowledge_queries += 1
        if ctx.knowledge.lookup(sig) is not None:
            self.knowledge_hits += 1
            return True
        return False

    @property
    def p_util(self) -> float:
        if self.knowledge_queries == 0:
            return 0.0
        return self.knowledge_hits / self.knowledge_queries

    @property
    def alu_saved_pct(self) -> float:
        total = self.kop_executed + self.kop_reused
        if total == 0:
            return 0.0
        return 100.0 * self.kop_reused / total

    def kop_breakdown(self) -> dict[str, dict[str, int | float | str]]:
        out: dict[str, dict[str, int | float | str]] = {}
        for kid, cpu in self._kop_cpu.items():
            kop_id = KOPId(kid)
            defn = KOP_DEFINITIONS.get(kop_id, {})
            out[str(kid)] = {
                "name": defn.get("name", "?"),
                "cpu_units": round(cpu, 2),
            }
        return out

    def to_dict(self, version: str = "0.1.2-alpha") -> dict:
        return {
            "version": version,
            "kop_executed": self.kop_executed,
            "kop_denied": self.kop_denied,
            "kop_reused": self.kop_reused,
            "knowledge_queries": self.knowledge_queries,
            "knowledge_hits": self.knowledge_hits,
            "p_util": round(self.p_util, 4),
            "alu_saved_pct": round(self.alu_saved_pct, 2),
            "alu_saved_units": round(self.alu_saved_units, 2),
            "elapsed_sec": round(time.perf_counter() - self._start, 4),
            "kop_breakdown": self.kop_breakdown(),
        }

    def to_json(self, version: str = "0.1.2-alpha") -> str:
        return json.dumps(self.to_dict(version=version), indent=2)

    def export(self, path: Path | str, version: str = "0.1.2-alpha") -> Path:
        dest = Path(path)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(self.to_json(version=version), encoding="utf-8")
        return dest