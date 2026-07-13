#!/usr/bin/env python3
"""
AntiPC Architecture E — real UDP network benchmark.

Spawns hub workers on localhost, floods WORK packets, measures throughput.

Usage:
    python udp_benchmark.py
    python udp_benchmark.py --packets 2000 --hubs 4
"""

from __future__ import annotations

import argparse
import os
import socket
import subprocess
import sys
import time
from dataclasses import dataclass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from architectures import BenchmarkResult, _make_packets
from udp_protocol import MsgType, pack_done, pack_ping, pack_work, unpack_result, RESULT


BASE_HUB_PORT = 19701
MASTER_PORT = 19700


@dataclass
class HubProc:
    port: int
    process: subprocess.Popen[bytes]


def _spawn_hubs(count: int, script_dir: str) -> list[HubProc]:
    hubs: list[HubProc] = []
    hub_script = os.path.join(script_dir, "hub_node.py")
    for i in range(count):
        port = BASE_HUB_PORT + i
        proc = subprocess.Popen(
            [sys.executable, hub_script, "--port", str(port),
             "--master-host", "127.0.0.1", "--master-port", str(MASTER_PORT)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        hubs.append(HubProc(port, proc))
    time.sleep(0.3)
    return hubs


def _stop_hubs(hubs: list[HubProc], sock: socket.socket) -> None:
    for hub in hubs:
        try:
            sock.sendto(pack_done(), ("127.0.0.1", hub.port))
        except OSError:
            pass
    time.sleep(0.2)
    for hub in hubs:
        hub.process.terminate()
        try:
            hub.process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            hub.process.kill()


def _wait_pong(sock: socket.socket, hub_port: int, timeout: float = 2.0) -> bool:
    deadline = time.perf_counter() + timeout
    sock.sendto(pack_ping(), ("127.0.0.1", hub_port))
    while time.perf_counter() < deadline:
        try:
            data, _ = sock.recvfrom(1024)
            if data and data[0] == MsgType.PONG:
                return True
        except (TimeoutError, BlockingIOError, OSError):
            pass
    return False


def run_udp_network(packets_count: int, payload_size: int, hubs: int) -> BenchmarkResult:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    packets = _make_packets(packets_count, payload_size)
    hub_procs = _spawn_hubs(hubs, script_dir)

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind(("127.0.0.1", MASTER_PORT))
    sock.settimeout(0.005)

    hub_addrs = [("127.0.0.1", BASE_HUB_PORT + i) for i in range(hubs)]
    for port in [h.port for h in hub_procs]:
        if not _wait_pong(sock, port):
            _stop_hubs(hub_procs, sock)
            sock.close()
            raise RuntimeError(f"hub on port {port} did not respond to PING")

    pending = packets_count
    received: dict[int, int] = {}
    next_send = 0

    t0 = time.perf_counter()
    while pending > 0 or next_send < len(packets):
        while next_send < len(packets):
            hub = hub_addrs[next_send % hubs]
            pkt = packets[next_send]
            sock.sendto(pack_work(pkt.seq, pkt.payload), hub)
            next_send += 1
            if next_send % 64 == 0:
                break

        try:
            data, _ = sock.recvfrom(1024)
        except (TimeoutError, BlockingIOError, OSError):
            continue

        if len(data) >= RESULT.size and data[0] == MsgType.RESULT:
            seq, _digest = unpack_result(data)
            if seq not in received:
                received[seq] = _digest
                pending -= 1

    elapsed = time.perf_counter() - t0
    _stop_hubs(hub_procs, sock)
    sock.close()

    return BenchmarkResult(
        f"E — UDP real ({hubs} hubs)",
        packets_count,
        elapsed,
        len(received),
        hubs=hubs,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="AntiPC UDP network benchmark")
    parser.add_argument("--packets", type=int, default=2000)
    parser.add_argument("--payload", type=int, default=64)
    parser.add_argument("--hubs", type=int, default=4)
    args = parser.parse_args()

    print()
    print("=" * 72)
    print("  AntiPC — Arquitectura E: UDP real (red local como cómputo)")
    print("=" * 72)
    print()
    print(f"  Hubs: {args.hubs}  |  Paquetes: {args.packets}  |  Payload: {args.payload} B")
    print()

    try:
        result = run_udp_network(args.packets, args.payload, args.hubs)
    except RuntimeError as exc:
        print(f"  ERROR: {exc}")
        return 1

    print(f"  Procesados : {result.processed} / {result.packets}")
    print(f"  Tiempo     : {result.elapsed_s:.4f} s")
    print(f"  Throughput : {result.throughput:.0f} pkt/s")
    print(f"  Hubs       : {result.hubs}")
    print()
    print("  La red deja de ser solo transporte: cada hub ejecuta HASH→INDEX")
    print("  con buffer lock-free; el maestro solo coordina por UDP.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())