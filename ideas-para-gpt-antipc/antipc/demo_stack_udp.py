#!/usr/bin/env python3
"""Demo: AntiPCStack enlazado a hubs UDP con HMAC."""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from architecture.antipc_stack import AntiPCStack
from hash_engine import k3_hash_buffer
from runtime.profile import ExecutionProfile


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(Path(__file__).parent.parent / "sale-it"))
    parser.add_argument("--hubs", type=int, default=4)
    parser.add_argument("--max-files", type=int, default=30)
    args = parser.parse_args()

    stack = AntiPCStack(profile=ExecutionProfile.cluster())

    print()
    print("=" * 72)
    print("  AntiPC — STACK + UDP HUBS (HMAC enlazado)")
    print("=" * 72)

    print("\n  [1] Permiso local HMAC...")
    if not stack.request_network_permission():
        print("  DENEGADO")
        return 1
    print("  OK")

    print(f"\n  [2] Arrancando {args.hubs} hubs UDP + autenticación...")
    try:
        net = stack.start_network(hubs=args.hubs)
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return 1
    print(f"  Hubs autenticados: {net['authenticated']}/{net['total']}")
    if net["authenticated"] == 0:
        stack.stop_network()
        return 1

    root = Path(args.root)
    files = [str(p) for p in root.rglob("*") if p.is_file()][: args.max_files]

    print(f"\n  [3] HASH remoto en hubs ({len(files)} ficheros)...")
    t0 = time.perf_counter()
    remote = stack.pipeline_remote_files(files)
    t_remote = time.perf_counter() - t0

    print(f"\n  [4] HASH local monolítico (misma carga)...")
    t0 = time.perf_counter()
    for path in files:
        with open(path, "rb") as f:
            k3_hash_buffer(f.read())
    t_local = time.perf_counter() - t0

    print("\n" + "-" * 72)
    print("  RESULTADOS")
    print("-" * 72)
    print(f"  Ficheros       : {len(remote)}")
    print(f"  Tiempo remoto  : {t_remote:.4f} s  (ALU en hubs)")
    print(f"  Tiempo local   : {t_local:.4f} s  (ALU en maestro)")
    if t_remote > 0:
        print(f"  Descarga maestro: {(1 - t_remote / max(t_local, 1e-9)) * 100:.1f}%")
    print(f"\n  Muestra:")
    for r in remote[:3]:
        print(f"    {r['digest_hex']}  {r['path']}")

    out = Path(__file__).parent / "output"
    stack.runtime.telemetry.export_csv(out / "telemetria_stack_udp.csv")
    print(f"\n  Telemetría: {out / 'telemetria_stack_udp.csv'}")

    stack.stop_network()
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())