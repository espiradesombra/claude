#!/usr/bin/env python3
"""
AntiPC — Adaptive Network Through Parallel Computing
Benchmark reproducible de 4 arquitecturas.

Uso:
    python benchmark.py
    python benchmark.py --packets 8000 --payload 256 --hubs 4
    python benchmark.py --udp          # incluye arquitectura E (UDP real)
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from architectures import (
    BenchmarkResult,
    _make_packets,
    run_conventional,
    run_distributed,
    run_grafcet,
    run_lockfree,
)


def _print_header() -> None:
    print()
    print("=" * 72)
    print("  AntiPC — Adaptive Network Through Parallel Computing")
    print('  Lema: "Transforming communication into computation."')
    print("=" * 72)
    print()


def _print_result(r: BenchmarkResult, baseline_s: float) -> None:
    rel = (r.elapsed_s / baseline_s) * 100.0 if baseline_s else 100.0
    bar = "#" * max(1, int(50 * baseline_s / max(r.elapsed_s, 1e-9)))
    print(f"  {r.name}")
    print(f"    Paquetes procesados : {r.processed:>8} / {r.packets}")
    print(f"    Tiempo              : {r.elapsed_s:>8.4f} s")
    print(f"    Throughput          : {r.throughput:>8.0f} pkt/s")
    print(f"    Carga ALU relativa  : {rel:>7.1f}%  (menor = mejor)")
    if r.cache_hits:
        print(f"    Cache hits (K(N))   : {r.cache_hits:>8}")
    if r.hubs > 1:
        print(f"    Hubs activos        : {r.hubs:>8}")
    print(f"    [{bar}]")
    print()


def run_benchmark(packets: int, payload: int, hubs: int, include_udp: bool = False) -> list[BenchmarkResult]:
    data = _make_packets(packets, payload)
    results: list[BenchmarkResult] = []

    print(f"  Generando {packets} paquetes sintéticos ({payload} bytes c/u)...")
    print(f"  Ratio de duplicados ~15% (simula conocimiento reutilizable K(N))")
    print()

    archs: list[tuple[str, object]] = [
        ("A", lambda: run_conventional(data)),
        ("B", lambda: run_lockfree(data)),
        ("C", lambda: run_distributed(data, hubs=hubs)),
        ("D", lambda: run_grafcet(data)),
    ]
    if include_udp:
        from udp_benchmark import run_udp_network
        archs.append(("E", lambda: run_udp_network(packets, payload, hubs)))

    for label, fn in archs:
        print(f"  Ejecutando arquitectura {label}...", end="", flush=True)
        t0 = time.perf_counter()
        result = fn()
        print(f" {time.perf_counter() - t0:.2f}s")
        results.append(result)

    return results


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="AntiPC benchmark suite")
    parser.add_argument("--packets", type=int, default=4000, help="número de paquetes")
    parser.add_argument("--payload", type=int, default=128, help="tamaño de payload en bytes")
    parser.add_argument("--hubs", type=int, default=4, help="hubs para arquitectura C/E")
    parser.add_argument("--udp", action="store_true", help="incluir arquitectura E (UDP real)")
    args = parser.parse_args(argv)

    _print_header()
    results = run_benchmark(args.packets, args.payload, args.hubs, include_udp=args.udp)

    baseline = results[0].elapsed_s
    print("-" * 72)
    print("  RESULTADOS")
    print("-" * 72)
    print()

    for r in results:
        _print_result(r, baseline)

    best = min(results, key=lambda r: r.elapsed_s)
    print("-" * 72)
    print(f"  Mejor arquitectura: {best.name}")
    print(f"  Mejora vs convencional: {(1 - best.elapsed_s / baseline) * 100:.1f}% menos tiempo ALU")
    print()
    print("  Fórmula AntiPC: P_util(N) = N · E(N) + K(N)")
    print("    N  = capacidad de cómputo agregada (hubs / workers)")
    print("    E(N) = eficiencia de coordinación (0..1)")
    print("    K(N) = valor del conocimiento reutilizable (cache hits)")
    print()
    print("  Ver explicacion_cientifica.txt para el marco formal completo.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())