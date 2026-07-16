"""Motor MMO: Grafcet (existencia + heavy K3) + mundo autoritativo."""

from __future__ import annotations

from dataclasses import dataclass, field

from grafcet import GrafcetEngine, Packet, redundant_hash
from hash_engine import k3_hash_buffer

try:
    from native_engine import k3_heavy_hash as _heavy_hash
except ImportError:
    from grafcet import _heavy_hash  # type: ignore[attr-defined]

from .game_protocol import unpack_state, verify_k3_tag
from .game_world import WorldState


@dataclass
class GameBatchStats:
    validated: int = 0
    rejected_existence: int = 0
    rejected_wire: int = 0


@dataclass
class MmoGrafcetEngine(GrafcetEngine):
    """Grafcet D extendido: tras validar existencia, aplica estado al mundo MMO."""

    world: WorldState = field(default_factory=WorldState)
    num_shards: int = 1
    shard_id: int | None = None
    rejected_existence: int = 0
    remote_heavy: object | None = None  # GameClusterFabric | None
    remote_validations: int = 0
    last_batch: GameBatchStats = field(default_factory=GameBatchStats)

    def _flush_batch(self) -> list[int]:
        outputs: list[int] = []
        if not self.batch:
            return outputs

        stats = GameBatchStats()
        batch_copy = list(self.batch)

        self._transition("full")
        for pkt in batch_copy:
            replicas = redundant_hash(pkt.payload)
            if len(set(replicas)) > 1:
                from grafcet import majority_vote

                replicas[2] = majority_vote(replicas)
            self.matrix.append_row(replicas)
        self._transition("ok")

        pending_remote: list[tuple[int, Packet, int]] = []

        for idx, pkt in enumerate(batch_copy):
            if not self.matrix.row_exists(idx):
                stats.rejected_existence += 1
                self.rejected_existence += 1
                continue

            if not verify_k3_tag(pkt.payload):
                stats.rejected_wire += 1
                self.world.rejected_wire += 1
                continue

            self._transition("ok")
            key = k3_hash_buffer(pkt.payload)
            if key in self.knowledge_cache:
                self.cache_hits += 1
                digest = self.knowledge_cache[key]
                game_pkt = unpack_state(pkt.payload)
                if game_pkt is not None:
                    self.world.apply_validated(
                        game_pkt,
                        digest,
                        shard_filter=self.shard_id,
                        num_shards=self.num_shards,
                    )
                    stats.validated += 1
                outputs.append(digest)
                self.processed += 1
            elif self.remote_heavy is not None:
                pending_remote.append((idx, pkt, key))
            else:
                digest = _heavy_hash(pkt.payload)
                self.knowledge_cache[key] = digest
                game_pkt = unpack_state(pkt.payload)
                if game_pkt is not None:
                    self.world.apply_validated(
                        game_pkt,
                        digest,
                        shard_filter=self.shard_id,
                        num_shards=self.num_shards,
                    )
                    stats.validated += 1
                outputs.append(digest)
                self.processed += 1

        if pending_remote and self.remote_heavy is not None:
            payloads = [item[1].payload for item in pending_remote]
            digests = self.remote_heavy.validate_heavy_batch(
                payloads, num_shards=self.num_shards
            )
            self.remote_validations += len(pending_remote)
            for (_, pkt, key), digest in zip(pending_remote, digests):
                self.knowledge_cache[key] = digest
                game_pkt = unpack_state(pkt.payload)
                if game_pkt is not None:
                    self.world.apply_validated(
                        game_pkt,
                        digest,
                        shard_filter=self.shard_id,
                        num_shards=self.num_shards,
                    )
                    stats.validated += 1
                outputs.append(digest)
                self.processed += 1

        self.last_batch = stats
        self._transition("done")
        self._transition("flush")
        self.batch.clear()
        return outputs


@dataclass
class GameGrafcetEngine:
    """Fachada MMO sobre MmoGrafcetEngine + contadores de ingestión."""

    batch_size: int = 32
    num_shards: int = 1
    shard_id: int | None = None
    remote_heavy: object | None = None
    engine: MmoGrafcetEngine = field(default_factory=MmoGrafcetEngine)
    packets_ingested: int = 0
    batches_flushed: int = 0

    def __post_init__(self) -> None:
        self.engine.batch_size = self.batch_size
        self.engine.num_shards = self.num_shards
        self.engine.shard_id = self.shard_id
        self.engine.remote_heavy = self.remote_heavy

    @property
    def world(self) -> WorldState:
        return self.engine.world

    def ingest_raw(self, payload: bytes, seq: int = 0) -> list[int]:
        self.packets_ingested += 1
        out = self.engine.ingest(Packet(seq=seq, payload=payload[:64]))
        if out:
            self.batches_flushed += 1
        return out

    def finalize(self) -> list[int]:
        out = self.engine.finalize()
        if out:
            self.batches_flushed += 1
        return out