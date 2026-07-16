#!/usr/bin/env python3
"""Hub esclavo ZypyZape — emite ráfagas UDP 64B al maestro (simula periferia)."""

from __future__ import annotations

import argparse
import os
import socket
import struct
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PORT = 3333
MAGIC = 0x5A59  # "ZY"


def make_payload(hub_id: int, seq: int, extra: bytes = b"") -> bytes:
    """64 bytes: magic(2) hub(2) seq(4) + payload hash-like."""
    head = struct.pack("!HHI", MAGIC, hub_id & 0xFFFF, seq & 0xFFFFFFFF)
    body = (head + extra)[:64]
    return body + b"\x00" * (64 - len(body))


def main() -> int:
    p = argparse.ArgumentParser(description="Hub UDP ZypyZape")
    p.add_argument("--master", default="127.0.0.1")
    p.add_argument("--port", type=int, default=PORT)
    p.add_argument("--hub-id", type=int, default=1)
    p.add_argument("--count", type=int, default=10_000)
    p.add_argument("--burst", type=int, default=64, help="paquetes por ráfaga")
    p.add_argument("--delay-ms", type=float, default=0.0)
    args = p.parse_args()

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    target = (args.master, args.port)
    sent = 0
    t0 = time.perf_counter()

    for seq in range(args.count):
        payload = make_payload(args.hub_id, seq, struct.pack("!Q", seq * 0x9E3779B1))
        sock.sendto(payload, target)
        sent += 1
        if args.burst and (seq + 1) % args.burst == 0 and args.delay_ms > 0:
            time.sleep(args.delay_ms / 1000.0)

    elapsed = time.perf_counter() - t0
    print(f"[hub {args.hub_id}] enviados {sent} pkts → {target[0]}:{target[1]} en {elapsed:.3f}s")
    sock.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())