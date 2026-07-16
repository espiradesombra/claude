"""Grafcet state engine with existence matrix and mathematical redundancy."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto

try:
    from .hash_engine import k3_hash
except ImportError:
    from hash_engine import k3_hash


class GrafcetStep(Enum):
    IDLE = auto()
    ACCUMULATE = auto()
    REDUNDANCY = auto()
    EXISTENCE = auto()
    PROCESS = auto()
    EMIT = auto()


@dataclass
class Packet:
    seq: int
    payload: bytes


@dataclass
class ExistenceMatrix:
    """Rows = packet batches, columns = redundant hash channels."""

    rows: list[list[int]] = field(default_factory=list)

    def append_row(self, hashes: list[int]) -> None:
        self.rows.append(hashes)

    def row_exists(self, row_idx: int, threshold: int = 2) -> bool:
        if row_idx >= len(self.rows):
            return False
        row = self.rows[row_idx]
        if len(row) < 2:
            return False
        parities = [h & 1 for h in row]
        return max(parities.count(0), parities.count(1)) >= threshold


def redundant_hash(payload: bytes, replicas: int = 3) -> list[int]:
    seeds = (0xA5A5A5A5, 0x5A5A5A5A, 0x12345678)
    return [k3_hash(payload, seed=seeds[i % len(seeds)]) for i in range(replicas)]


def majority_vote(values: list[int]) -> int:
    return max(set(values), key=values.count)


def _heavy_hash(payload: bytes) -> int:
    digest = k3_hash(payload)
    for _ in range(4):
        digest = k3_hash(digest.to_bytes(4, "little") + payload)
    return digest


@dataclass
class GrafcetEngine:
    step: GrafcetStep = GrafcetStep.IDLE
    batch: list[Packet] = field(default_factory=list)
    matrix: ExistenceMatrix = field(default_factory=ExistenceMatrix)
    batch_size: int = 32
    knowledge_cache: dict[int, int] = field(default_factory=dict)
    processed: int = 0
    cache_hits: int = 0

    def _transition(self, event: str) -> None:
        table: dict[tuple[GrafcetStep, str], GrafcetStep] = {
            (GrafcetStep.IDLE, "packet"): GrafcetStep.ACCUMULATE,
            (GrafcetStep.ACCUMULATE, "full"): GrafcetStep.REDUNDANCY,
            (GrafcetStep.REDUNDANCY, "ok"): GrafcetStep.EXISTENCE,
            (GrafcetStep.EXISTENCE, "ok"): GrafcetStep.PROCESS,
            (GrafcetStep.PROCESS, "done"): GrafcetStep.EMIT,
            (GrafcetStep.EMIT, "flush"): GrafcetStep.IDLE,
        }
        self.step = table.get((self.step, event), self.step)

    def _flush_batch(self) -> list[int]:
        outputs: list[int] = []
        if not self.batch:
            return outputs

        self._transition("full")
        replica_rows: list[list[int]] = []
        for pkt in self.batch:
            replicas = redundant_hash(pkt.payload)
            if len(set(replicas)) > 1:
                replicas[2] = majority_vote(replicas)
            replica_rows.append(replicas)
            self.matrix.append_row(replicas)
        self._transition("ok")

        for idx, pkt in enumerate(self.batch):
            if not self.matrix.row_exists(idx):
                continue
            self._transition("ok")
            key = k3_hash(pkt.payload)
            if key in self.knowledge_cache:
                self.cache_hits += 1
                digest = self.knowledge_cache[key]
            else:
                digest = _heavy_hash(pkt.payload)
                self.knowledge_cache[key] = digest
            outputs.append(digest)
            self.processed += 1

        self._transition("done")
        self._transition("flush")
        self.batch.clear()
        return outputs

    def ingest(self, packet: Packet) -> list[int]:
        if self.step == GrafcetStep.IDLE:
            self.batch.clear()
            self.matrix = ExistenceMatrix()
            self._transition("packet")

        self.batch.append(packet)
        if len(self.batch) < self.batch_size:
            return []

        return self._flush_batch()

    def finalize(self) -> list[int]:
        """Process trailing packets smaller than batch_size."""
        return self._flush_batch()