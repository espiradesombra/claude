"""
B → D pipeline: UDP slot-ring ingest (ZypyZape) → GrafcetEngine (existencia + reuse).

Arquitectura B recibe con recvinto(slot fijo); el worker alimenta Grafcet (D) sin
segunda copia usuario entre red y motor de estados.
"""

from __future__ import annotations

import json
import socket
import struct
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from grafcet import GrafcetEngine, Packet

_ANTIPC_ROOT = Path(__file__).resolve().parents[3]
if str(_ANTIPC_ROOT) not in sys.path:
    sys.path.insert(0, str(_ANTIPC_ROOT))

from zypyzape.slot_ring import SlotRing  # noqa: E402

PORT = 3333
MAGIC = 0x5A59


@dataclass
class BdDemoMetrics:
    architecture: str
    packets_in: int
    packets_grafcet: int
    drops: int
    memcpy_user_copies: int
    cache_hits: int
    grafcet_step: str
    elapsed_s: float
    hubs: int = 1
    batch_size: int = 32

    @property
    def throughput_pps(self) -> float:
        return self.packets_grafcet / self.elapsed_s if self.elapsed_s else 0.0


def _parse_seq(payload: memoryview | bytes) -> int:
    data = bytes(payload[:8])
    if len(data) >= 8:
        magic, _hub, seq = struct.unpack("!HHI", data[:8])
        if magic == MAGIC:
            return seq
    return 0


def _bind_master(port: int) -> socket.socket:
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind(("0.0.0.0", port))
    sock.setblocking(False)
    return sock


def spawn_hubs(master: str, port: int, hubs: int, count: int) -> list[threading.Thread]:
    emitter = _ANTIPC_ROOT / "zypyzape" / "hub_emitter.py"
    per_hub = max(1, count // hubs)
    threads: list[threading.Thread] = []

    def launch(hid: int) -> None:
        subprocess.run(
            [
                sys.executable,
                str(emitter),
                "--master",
                master,
                "--port",
                str(port),
                "--hub-id",
                str(hid),
                "--count",
                str(per_hub),
                "--burst",
                "128",
            ],
            check=False,
            capture_output=True,
        )

    for hid in range(1, hubs + 1):
        threads.append(threading.Thread(target=launch, args=(hid,), daemon=True))
    return threads


def run_bd_pipeline(
    sock: socket.socket,
    duration_s: float,
    *,
    batch_size: int = 32,
) -> BdDemoMetrics:
    """Slot ring (B) ingest loop; worker pushes packets into GrafcetEngine (D)."""
    ring = SlotRing(capacity=1024)
    engine = GrafcetEngine(batch_size=batch_size)
    running = threading.Event()
    running.set()
    received = 0
    drops = 0
    copies = 0
    ingested = 0

    def worker() -> None:
        nonlocal ingested
        while running.is_set() or len(ring) > 0:
            idx = ring.read_index()
            if idx is None:
                time.sleep(0)
                continue
            view = ring.slot_view(idx)
            payload = bytes(view)
            seq = _parse_seq(payload)
            engine.ingest(Packet(seq=seq, payload=payload))
            ring.commit_read()
            ingested += 1

    t0 = time.perf_counter()
    th = threading.Thread(target=worker, daemon=True)
    th.start()

    while time.perf_counter() - t0 < duration_s:
        widx = ring.write_index()
        if widx is None:
            drops += 1
            try:
                sock.recv(2048)
            except BlockingIOError:
                pass
            time.sleep(0.00005)
            continue
        view = ring.slot_view(widx)
        try:
            n, _ = sock.recvfrom_into(view)
        except BlockingIOError:
            time.sleep(0.001)
            continue
        except TimeoutError:
            continue
        if n <= 0:
            continue
        received += 1
        copies += 1
        ring.commit_write()

    running.clear()
    th.join(timeout=2.0)
    engine.finalize()
    elapsed = time.perf_counter() - t0

    return BdDemoMetrics(
        architecture="B_slot_ring_to_D_grafcet",
        packets_in=received,
        packets_grafcet=engine.processed,
        drops=drops,
        memcpy_user_copies=copies,
        cache_hits=engine.cache_hits,
        grafcet_step=engine.step.name,
        elapsed_s=elapsed,
        batch_size=batch_size,
    )


def run_network_demo(
    *,
    port: int = PORT,
    duration_s: float = 3.0,
    hubs: int = 4,
    packets: int = 20_000,
    batch_size: int = 32,
    spawn: bool = True,
    out_path: Path | None = None,
) -> BdDemoMetrics:
    sock = _bind_master(port)
    if spawn:
        threads = spawn_hubs("127.0.0.1", port, hubs, packets)
        time.sleep(0.15)
        for th in threads:
            th.start()
    metrics = run_bd_pipeline(sock, duration_s, batch_size=batch_size)
    metrics.hubs = hubs
    sock.close()

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "port": port,
            "hubs": hubs,
            "batch_size": batch_size,
            "pipeline": "B_slot_ring → D_grafcet",
            "metrics": asdict(metrics),
        }
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return metrics