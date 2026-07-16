"""Bots MMO — simulan jugadores enviando estados UDP 64B."""

from __future__ import annotations

import math
import random
import socket
import threading
import time
from dataclasses import dataclass

from .game_protocol import GameMsgType, pack_state

DEFAULT_PORT = 3344


@dataclass
class BotMetrics:
    sent: int = 0
    errors: int = 0


def _bot_loop(
    *,
    master: str,
    port: int,
    bot_id: int,
    num_bots: int,
    ticks: int,
    tick_hz: float,
    duplicate_ratio: float,
    stop: threading.Event,
    metrics: BotMetrics,
) -> None:
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    target = (master, port)
    rng = random.Random(bot_id * 7919)
    entity_id = 1000 + bot_id
    shard_id = bot_id % max(1, num_bots // 4 + 1)
    px = rng.uniform(-500, 500)
    py = 0.0
    pz = rng.uniform(-500, 500)
    angle = rng.uniform(0, math.tau)
    speed = rng.uniform(2.0, 8.0)
    last_payload: bytes | None = None
    interval = 1.0 / max(tick_hz, 1.0)

    for tick in range(ticks):
        if stop.is_set():
            break
        angle += rng.uniform(-0.3, 0.3)
        vx = math.cos(angle) * speed
        vz = math.sin(angle) * speed
        px += vx * interval
        pz += vz * interval

        if rng.random() < duplicate_ratio and last_payload is not None:
            payload = last_payload
        else:
            payload = pack_state(
                entity_id=entity_id,
                tick=tick,
                shard_id=shard_id,
                px=px,
                py=py,
                pz=pz,
                vx=vx,
                vy=0.0,
                vz=vz,
                seq=tick,
                yaw=int(math.degrees(angle)) % 360,
                msg_type=GameMsgType.MOVE,
            )
            last_payload = payload

        try:
            sock.sendto(payload, target)
            metrics.sent += 1
        except OSError:
            metrics.errors += 1

        if tick % 64 == 0:
            hb = pack_state(
                entity_id=entity_id,
                tick=tick,
                shard_id=shard_id,
                px=px,
                py=py,
                pz=pz,
                msg_type=GameMsgType.HEARTBEAT,
            )
            try:
                sock.sendto(hb, target)
                metrics.sent += 1
            except OSError:
                metrics.errors += 1

        time.sleep(interval * rng.uniform(0.8, 1.2))

    sock.close()


def spawn_bots(
    *,
    master: str = "127.0.0.1",
    port: int = DEFAULT_PORT,
    players: int = 64,
    ticks: int = 500,
    tick_hz: float = 20.0,
    duplicate_ratio: float = 0.12,
) -> tuple[list[threading.Thread], list[BotMetrics], threading.Event]:
    """Lanza N hilos bot; retorna threads, métricas y evento de parada."""
    stop = threading.Event()
    threads: list[threading.Thread] = []
    all_metrics: list[BotMetrics] = []

    for bot_id in range(players):
        m = BotMetrics()
        all_metrics.append(m)
        th = threading.Thread(
            target=_bot_loop,
            kwargs={
                "master": master,
                "port": port,
                "bot_id": bot_id,
                "num_bots": players,
                "ticks": ticks,
                "tick_hz": tick_hz,
                "duplicate_ratio": duplicate_ratio,
                "stop": stop,
                "metrics": m,
            },
            daemon=True,
            name=f"mmo-bot-{bot_id}",
        )
        threads.append(th)
    return threads, all_metrics, stop


def run_bot_swarm(
    *,
    master: str = "127.0.0.1",
    port: int = DEFAULT_PORT,
    players: int = 64,
    duration_s: float = 5.0,
    tick_hz: float = 20.0,
    duplicate_ratio: float = 0.12,
) -> int:
    """CLI standalone: bots durante duration_s."""
    ticks = max(1, int(duration_s * tick_hz))
    threads, metrics, stop = spawn_bots(
        master=master,
        port=port,
        players=players,
        ticks=ticks,
        tick_hz=tick_hz,
        duplicate_ratio=duplicate_ratio,
    )
    for th in threads:
        th.start()
    time.sleep(duration_s)
    stop.set()
    for th in threads:
        th.join(timeout=2.0)
    total = sum(m.sent for m in metrics)
    print(f"[bots] {players} jugadores, {total} paquetes enviados, errores={sum(m.errors for m in metrics)}")
    return total