#!/usr/bin/env python3
"""
ZypyZape / AntiPC — Demostración de viabilidad UDP (Windows HP Envy).

Compara:
  A) Convencional: recvinto(temp) + memcpy a cola
  B) AntiPC:      recvinto(slot fijo del ring) — sin segunda copia en usuario

Incluye benchmark con hubs simulados y métricas exportables.
"""

from __future__ import annotations

import argparse
import json
import os
import queue
import socket
import struct
import sys
import threading
import time
from dataclasses import dataclass, asdict
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hash_engine import k3_hash_buffer
from zypyzape.slot_ring import DATA_SIZE, SlotRing

PORT = 3333
MAGIC = 0x5A59


@dataclass
class RunMetrics:
    mode: str
    packets_in: int
    packets_processed: int
    drops: int
    memcpy_user_copies: int
    elapsed_s: float
    hubs: int = 1

    @property
    def throughput_pps(self) -> float:
        return self.packets_processed / self.elapsed_s if self.elapsed_s else 0.0

    @property
    def alu_proxy_ms_per_pkt(self) -> float:
        return (self.elapsed_s * 1000.0) / max(self.packets_processed, 1)


def _light_process(payload: memoryview | bytes) -> int:
    """Procesamiento tipo Grafcet/K3 — ALU acotada por paquete."""
    data = bytes(payload[: min(len(payload), 64)])
    d = k3_hash_buffer(data)
    return d & 0xFFFF


def run_arch_a(sock: socket.socket, duration_s: float, qsize: int = 4096) -> RunMetrics:
    """Arquitectura A: temp buffer + copia explícita a Queue."""
    q: queue.Queue[bytes] = queue.Queue(maxsize=qsize)
    temp = bytearray(2048)
    running = threading.Event()
    running.set()
    processed = 0
    received = 0
    drops = 0
    copies = 0

    def worker() -> None:
        nonlocal processed
        while running.is_set() or not q.empty():
            try:
                data = q.get(timeout=0.01)
            except queue.Empty:
                continue
            _light_process(data)
            processed += 1
            q.task_done()

    t0 = time.perf_counter()
    th = threading.Thread(target=worker, daemon=True)
    th.start()

    while time.perf_counter() - t0 < duration_s:
        try:
            n, _ = sock.recvfrom_into(temp)
        except BlockingIOError:
            time.sleep(0.001)
            continue
        except TimeoutError:
            continue
        if n < 8:
            continue
        received += 1
        copies += 1  # memcpy usuario: temp → objeto en cola
        try:
            q.put_nowait(bytes(temp[:n]))
        except queue.Full:
            drops += 1

    running.clear()
    q.join()
    th.join(timeout=2.0)
    elapsed = time.perf_counter() - t0
    copies += received  # recvinto kernel→user cuenta como 1 copia por paquete
    return RunMetrics("A_convencional_memcpy_cola", received, processed, drops, copies, elapsed)


def run_arch_b(sock: socket.socket, duration_s: float) -> RunMetrics:
    """Arquitectura B: recvinto directo al slot del ring (AntiPC/ZypyZape)."""
    ring = SlotRing(capacity=1024)
    running = threading.Event()
    running.set()
    processed = 0
    received = 0
    drops = 0
    copies = 0

    def worker() -> None:
        nonlocal processed
        while running.is_set() or len(ring) > 0:
            idx = ring.read_index()
            if idx is None:
                time.sleep(0)
                continue
            view = ring.slot_view(idx)
            _light_process(view)
            ring.commit_read()
            processed += 1

    t0 = time.perf_counter()
    th = threading.Thread(target=worker, daemon=True)
    th.start()

    while time.perf_counter() - t0 < duration_s:
        widx = ring.write_index()
        if widx is None:
            drops += 1
            try:
                sock.recv(2048)  # descartar si lleno
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
        copies += 1  # solo kernel→slot (sin segunda copia usuario)
        ring.commit_write()

    running.clear()
    th.join(timeout=2.0)
    elapsed = time.perf_counter() - t0
    return RunMetrics("B_antipc_slot_ring", received, processed, drops, copies, elapsed)


def _bind_master(port: int) -> socket.socket:
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind(("0.0.0.0", port))
    sock.setblocking(False)
    return sock


def spawn_hubs(master: str, port: int, hubs: int, count: int) -> list[threading.Thread]:
    import subprocess

    root = Path(__file__).resolve().parent
    emitter = root / "hub_emitter.py"
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
        th = threading.Thread(target=launch, args=(hid,), daemon=True)
        threads.append(th)
    return threads


def main() -> int:
    p = argparse.ArgumentParser(description="Viabilidad UDP ZypyZape/AntiPC")
    p.add_argument("--port", type=int, default=PORT)
    p.add_argument("--duration", type=float, default=3.0, help="segundos escucha por modo")
    p.add_argument("--hubs", type=int, default=4)
    p.add_argument("--packets", type=int, default=20_000, help="total paquetes de hubs")
    p.add_argument("--out", default="", help="JSON resultados")
    p.add_argument("--no-hubs", action="store_true", help="solo escucha (lanza hubs aparte)")
    args = p.parse_args()

    print()
    print("=" * 72)
    print("  ZypyZape / AntiPC — Viabilidad UDP (inyección a slots)")
    print("=" * 72)
    print(f"  Puerto {args.port} | hubs {args.hubs} | {args.packets} pkts | {args.duration}s/modo")
    print()

    results: list[RunMetrics] = []

    for mode_name, runner in [("A", run_arch_a), ("B", run_arch_b)]:
        sock = _bind_master(args.port)
        print(f"  [{mode_name}] Escuchando...")
        if not args.no_hubs:
            threads = spawn_hubs("127.0.0.1", args.port, args.hubs, args.packets)
            time.sleep(0.15)
            for th in threads:
                th.start()
        else:
            print("       (esperando hubs externos — python hub_emitter.py ...)")

        m = runner(sock, args.duration)
        m.hubs = args.hubs
        results.append(m)
        sock.close()
        time.sleep(0.3)
        print(
            f"       in={m.packets_in} proc={m.packets_processed} "
            f"drop={m.drops} copies={m.memcpy_user_copies} "
            f"{m.throughput_pps:.0f} pkt/s  proxy_ALU={m.alu_proxy_ms_per_pkt:.3f} ms/pkt"
        )

    if len(results) == 2:
        a, b = results
        if a.throughput_pps > 0:
            gain = (b.throughput_pps - a.throughput_pps) / a.throughput_pps * 100
            copy_save = (
                (1 - b.memcpy_user_copies / max(a.memcpy_user_copies, 1)) * 100
            )
            print()
            print("-" * 72)
            print("  COMPARATIVA")
            print("-" * 72)
            print(f"  Throughput B vs A     : {gain:+.1f}%")
            print(f"  Copias usuario ahorro : {copy_save:.1f}% (B evita memcpy→cola)")
            if b.alu_proxy_ms_per_pkt < a.alu_proxy_ms_per_pkt:
                print(
                    f"  Proxy carga ALU       : B {b.alu_proxy_ms_per_pkt:.3f} ms "
                    f"< A {a.alu_proxy_ms_per_pkt:.3f} ms"
                )

    out_path = Path(args.out) if args.out else Path(__file__).parent.parent / "output" / "zypyzape_viability.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "port": args.port,
        "hubs": args.hubs,
        "runs": [asdict(r) for r in results],
    }
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\n  JSON: {out_path}")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())