#!/usr/bin/env python3
"""
AntiPC v0.1.0 — Flow Kernel demo
Implements gptcomputing.txt roadmap: UPS + Event Bus + Knowledge Buffer.
"""

from __future__ import annotations

import argparse
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plugins.k3_plugin import K3HashPlugin
from runtime.kernel import FlowKernel
from runtime.modes import BatteryMode, BenchmarkMode, LocalMode, ResearchMode


def run_demo(count: int, repeat: float, mode_name: str) -> None:
    modes = {
        "local": LocalMode(),
        "battery": BatteryMode(),
        "research": ResearchMode(),
        "benchmark": BenchmarkMode(),
    }
    mode = modes[mode_name]

    kernel = FlowKernel(mode=mode)
    kernel.register_plugin(K3HashPlugin())
    if mode_name == "benchmark":
        kernel.network_permission = True

    rng = random.Random(42)
    pool = [rng.randbytes(64) for _ in range(20)]

    t0 = time.perf_counter()
    for i in range(count):
        payload = pool[i % len(pool)] if rng.random() < repeat else rng.randbytes(64)
        raw = kernel.ingest_raw(payload, label=f"pkt_{i}")
        kernel.submit_operation("K3_HASH", [raw])
        kernel.run_until_idle(max_steps=4)

    elapsed = time.perf_counter() - t0
    stats = kernel.stats()

    print()
    print("=" * 68)
    print("  AntiPC v0.1.0 — Flow Kernel")
    print("=" * 68)
    print(f"  Modo              : {stats['mode']}")
    print(f"  Operaciones       : {count}")
    print(f"  Ejecutadas ALU    : {stats['executed']}")
    print(f"  Canceladas (reuse): {stats['cancelled']}")
    print(f"  Knowledge hits    : {stats['knowledge_hits']} / {stats['knowledge_queries']}")
    print(f"  Event coalescing  : {stats['events_coalesced']}")
    print(f"  Tiempo            : {elapsed:.4f} s")
    if stats['executed']:
        print(f"  ALU ahorrada      : {(1 - stats['executed'] / count) * 100:.1f}%")
    print()
    print("  P_util(N) = N·E(N) + K(N)  —  K(N) = knowledge hits")
    print("=" * 68)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=500)
    parser.add_argument("--repeat", type=float, default=0.3, help="ratio duplicados")
    parser.add_argument("--mode", choices=["local", "battery", "research", "benchmark"],
                        default="battery")
    args = parser.parse_args()
    run_demo(args.count, args.repeat, args.mode)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())