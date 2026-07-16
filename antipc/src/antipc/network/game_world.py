"""Estado del mundo MMO — entidades, shards, dedup de estado."""

from __future__ import annotations

from dataclasses import dataclass, field

from .game_protocol import GameStatePacket


@dataclass
class EntityState:
    entity_id: int
    shard_id: int
    tick: int
    px: float
    py: float
    pz: float
    vx: float
    vy: float
    vz: float
    digest: int
    k3_tag: int
    updates: int = 0


@dataclass
class WorldState:
    """Mundo autoritativo post-validación Grafcet."""

    entities: dict[int, EntityState] = field(default_factory=dict)
    moves_applied: int = 0
    state_cache_hits: int = 0
    rejected_wire: int = 0
    rejected_stale: int = 0
    max_tick: int = 0

    def apply_validated(
        self,
        pkt: GameStatePacket,
        digest: int,
        *,
        shard_filter: int | None = None,
        num_shards: int = 1,
    ) -> bool:
        if shard_filter is not None and num_shards > 1:
            if pkt.entity_id % num_shards != shard_filter:
                return False

        prev = self.entities.get(pkt.entity_id)
        if prev is not None:
            if pkt.k3_tag == prev.k3_tag and pkt.tick <= prev.tick:
                self.state_cache_hits += 1
                return False
            if pkt.tick < prev.tick:
                self.rejected_stale += 1
                return False

        self.entities[pkt.entity_id] = EntityState(
            entity_id=pkt.entity_id,
            shard_id=pkt.shard_id,
            tick=pkt.tick,
            px=pkt.px,
            py=pkt.py,
            pz=pkt.pz,
            vx=pkt.vx,
            vy=pkt.vy,
            vz=pkt.vz,
            digest=digest,
            k3_tag=pkt.k3_tag,
            updates=(prev.updates + 1) if prev else 1,
        )
        self.moves_applied += 1
        if pkt.tick > self.max_tick:
            self.max_tick = pkt.tick
        return True

    def snapshot(self, limit: int = 5) -> list[dict]:
        rows = []
        for eid in sorted(self.entities.keys())[:limit]:
            e = self.entities[eid]
            rows.append(
                {
                    "entity": eid,
                    "shard": e.shard_id,
                    "tick": e.tick,
                    "pos": (round(e.px, 2), round(e.py, 2), round(e.pz, 2)),
                    "digest": f"{e.digest:08X}",
                }
            )
        return rows