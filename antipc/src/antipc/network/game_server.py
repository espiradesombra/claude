"""
Servidor MMO AntiPC — arquitectura B (SlotRing recvinto) → D (Grafcet MMO).

Un proceso autoritativo ingiere estados UDP de miles de jugadores simulados,
valida existencia, deduplica por K3 y mantiene WorldState por shard.
"""

from __future__ import annotations

import json
import socket
import sys
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path

_ANTIPC_SRC = Path(__file__).resolve().parents[1]
if str(_ANTIPC_SRC) not in sys.path:
    sys.path.insert(0, str(_ANTIPC_SRC))

from zypyzape.slot_ring import SlotRing  # noqa: E402

from .game_bots import spawn_bots  # noqa: E402
from .game_engine import GameGrafcetEngine  # noqa: E402
from .game_protocol import GAME_PKT_SIZE, unpack_state  # noqa: E402

GAME_PORT = 3344


@dataclass
class GameMmoMetrics:
    architecture: str
    port: int
    players: int
    duration_s: float
    packets_in: int
    packets_ingested: int
    moves_applied: int
    entities_active: int
    drops: int
    memcpy_user_copies: int
    grafcet_cache_hits: int
    state_cache_hits: int
    rejected_existence: int
    rejected_wire: int
    rejected_stale: int
    batches_flushed: int
    max_tick: int
    elapsed_s: float
    batch_size: int = 32
    num_shards: int = 1
    hubs: int = 0
    hubs_authenticated: int = 0
    remote_validations: int = 0
    remote_validations_ok: int = 0
    remote_fallback: int = 0

    @property
    def throughput_pps(self) -> float:
        return self.packets_in / self.elapsed_s if self.elapsed_s else 0.0

    @property
    def moves_per_sec(self) -> float:
        return self.moves_applied / self.elapsed_s if self.elapsed_s else 0.0


def _bind_game_server(port: int) -> socket.socket:
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind(("0.0.0.0", port))
    sock.setblocking(False)
    return sock


def run_game_server(
    sock: socket.socket,
    duration_s: float,
    *,
    batch_size: int = 32,
    num_shards: int = 1,
    shard_id: int | None = None,
) -> GameMmoMetrics:
    game = GameGrafcetEngine(
        batch_size=batch_size,
        num_shards=num_shards,
        shard_id=shard_id,
    )
    return _run_game_server_with_engine(sock, duration_s, game=game)


def run_game_demo(
    *,
    port: int = GAME_PORT,
    duration_s: float = 5.0,
    players: int = 128,
    tick_hz: float = 20.0,
    batch_size: int = 32,
    num_shards: int = 4,
    duplicate_ratio: float = 0.15,
    spawn_bots_flag: bool = True,
    out_path: Path | None = None,
) -> tuple[GameMmoMetrics, list[dict]]:
    sock = _bind_game_server(port)
    bot_threads: list[threading.Thread] = []
    stop = threading.Event()
    game = GameGrafcetEngine(batch_size=batch_size, num_shards=num_shards)

    if spawn_bots_flag:
        ticks = max(1, int(duration_s * tick_hz))
        bot_threads, _, stop = spawn_bots(
            master="127.0.0.1",
            port=port,
            players=players,
            ticks=ticks,
            tick_hz=tick_hz,
            duplicate_ratio=duplicate_ratio,
        )
        time.sleep(0.2)
        for th in bot_threads:
            th.start()

    metrics = _run_game_server_with_engine(
        sock, duration_s, game=game,
    )
    metrics.port = port
    metrics.players = players
    metrics.duration_s = duration_s
    world_sample = game.world.snapshot(limit=8)

    stop.set()
    for th in bot_threads:
        th.join(timeout=2.0)
    sock.close()

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "pipeline": "B_slot_ring → D_mmo_grafcet",
            "mmo": True,
            "metrics": asdict(metrics),
            "throughput_pps": metrics.throughput_pps,
            "moves_per_sec": metrics.moves_per_sec,
            "world_sample": world_sample,
        }
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return metrics, world_sample


def run_game_cluster_demo(
    *,
    port: int = GAME_PORT,
    duration_s: float = 5.0,
    players: int = 128,
    tick_hz: float = 20.0,
    batch_size: int = 32,
    num_shards: int = 4,
    hubs: int = 4,
    hub_hosts: list[str] | None = None,
    hub_base_port: int = 19701,
    duplicate_ratio: float = 0.15,
    spawn_bots_flag: bool = True,
    out_path: Path | None = None,
) -> tuple[GameMmoMetrics, list[dict]]:
    """MMO con validación heavy K3 en hubs L3 WORK/RESULT por shard."""
    from .game_cluster import GameClusterFabric

    cluster = GameClusterFabric(
        num_shards=num_shards,
        hubs=hubs,
        hub_hosts=hub_hosts,
        hub_base_port=hub_base_port,
    )
    hub_info = cluster.start()

    sock = _bind_game_server(port)
    bot_threads: list[threading.Thread] = []
    stop = threading.Event()
    game = GameGrafcetEngine(
        batch_size=batch_size,
        num_shards=num_shards,
        remote_heavy=cluster,
    )

    if spawn_bots_flag:
        ticks = max(1, int(duration_s * tick_hz))
        bot_threads, _, stop = spawn_bots(
            master="127.0.0.1",
            port=port,
            players=players,
            ticks=ticks,
            tick_hz=tick_hz,
            duplicate_ratio=duplicate_ratio,
        )
        time.sleep(0.35)
        for th in bot_threads:
            th.start()

    metrics = _run_game_server_with_engine(sock, duration_s, game=game)
    metrics.architecture = "B_slot_ring_L3_hubs_D_mmo_grafcet"
    metrics.port = port
    metrics.players = players
    metrics.duration_s = duration_s
    metrics.hubs = hub_info["hubs"]
    metrics.hubs_authenticated = hub_info["authenticated"]
    metrics.remote_validations = cluster.validations_sent
    metrics.remote_validations_ok = cluster.validations_ok
    metrics.remote_fallback = cluster.validations_fallback
    world_sample = game.world.snapshot(limit=8)

    stop.set()
    for th in bot_threads:
        th.join(timeout=2.0)
    sock.close()
    cluster.shutdown()

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "pipeline": "B_slot_ring → L3_hubs WORK/RESULT → D_mmo_grafcet",
            "mmo_cluster": True,
            "hub_info": hub_info,
            "metrics": asdict(metrics),
            "throughput_pps": metrics.throughput_pps,
            "moves_per_sec": metrics.moves_per_sec,
            "world_sample": world_sample,
        }
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return metrics, world_sample


def _run_game_server_with_engine(
    sock: socket.socket,
    duration_s: float,
    *,
    game: GameGrafcetEngine,
) -> GameMmoMetrics:
    ring = SlotRing(capacity=4096, data_size=GAME_PKT_SIZE)
    running = threading.Event()
    running.set()
    received = 0
    drops = 0
    copies = 0

    def worker() -> None:
        while running.is_set() or len(ring) > 0:
            idx = ring.read_index()
            if idx is None:
                time.sleep(0)
                continue
            view = ring.slot_view(idx)
            payload = bytes(view)
            pkt = unpack_state(payload)
            seq = pkt.seq if pkt else 0
            game.ingest_raw(payload, seq=seq)
            ring.commit_read()

    t0 = time.perf_counter()
    th = threading.Thread(target=worker, daemon=True, name="mmo-grafcet-worker")
    th.start()

    while time.perf_counter() - t0 < duration_s:
        widx = ring.write_index()
        if widx is None:
            drops += 1
            try:
                sock.recv(2048)
            except BlockingIOError:
                pass
            time.sleep(0.00002)
            continue
        view = ring.slot_view(widx)
        try:
            n, _ = sock.recvfrom_into(view)
        except BlockingIOError:
            time.sleep(0.0005)
            continue
        if n < GAME_PKT_SIZE:
            continue
        received += 1
        copies += 1
        ring.commit_write()

    running.clear()
    th.join(timeout=3.0)
    game.finalize()
    elapsed = time.perf_counter() - t0
    eng = game.engine

    return GameMmoMetrics(
        architecture="B_slot_ring_to_D_mmo_grafcet",
        port=0,
        players=0,
        duration_s=duration_s,
        packets_in=received,
        packets_ingested=game.packets_ingested,
        moves_applied=game.world.moves_applied,
        entities_active=len(game.world.entities),
        drops=drops,
        memcpy_user_copies=copies,
        grafcet_cache_hits=eng.cache_hits,
        state_cache_hits=game.world.state_cache_hits,
        rejected_existence=eng.rejected_existence,
        rejected_wire=game.world.rejected_wire,
        rejected_stale=game.world.rejected_stale,
        batches_flushed=game.batches_flushed,
        max_tick=game.world.max_tick,
        elapsed_s=elapsed,
        batch_size=game.batch_size,
        num_shards=game.num_shards,
    )