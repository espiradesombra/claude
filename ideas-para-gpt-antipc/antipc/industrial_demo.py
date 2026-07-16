#!/usr/bin/env python3
"""
Demo industrial AntiPC — inferido de gptcomputing.txt

Caso real: auditoría de ficheros (huella K3 + dedup).
Compara:
  1. Monolítico (PC clásico: hash todo siempre)
  2. Industrial cluster (resolver + knowledge + workers)
  3. Segunda pasada (demuestra K(N) — conocimiento reutilizable)
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hash_engine import k3_hash_buffer, get_backend
from industrial.runtime import IndustrialRuntime
from plugins.k3_dedup_plugin import K3DedupPlugin
from plugins.k3_file_plugin import K3FilePlugin
from runtime.profile import ExecutionProfile


def collect_files(root: Path, max_files: int) -> list[str]:
    files: list[str] = []
    for p in root.rglob("*"):
        if p.is_file() and not p.name.startswith("."):
            files.append(str(p))
            if len(files) >= max_files:
                break
    return files


def run_monolithic(paths: list[str]) -> dict:
    """PC clásico: memcpy + hash secuencial, sin knowledge buffer."""
    t0 = time.perf_counter()
    refs_meta = []
    for path in paths:
        try:
            with open(path, "rb") as f:
                data = f.read()
            digest = k3_hash_buffer(data)
            size = len(data)
            refs_meta.append((path, size, digest))
        except OSError:
            pass
    elapsed = time.perf_counter() - t0

    groups: dict[str, int] = {}
    for path, size, digest in refs_meta:
        key = f"{size}:{digest:08X}"
        groups[key] = groups.get(key, 0) + 1
    dup_groups = sum(1 for c in groups.values() if c > 1)

    return {
        "name": "Monolítico (PC clásico)",
        "files": len(refs_meta),
        "elapsed_s": elapsed,
        "alu_executions": len(refs_meta),
        "duplicate_groups": dup_groups,
        "knowledge_hits": 0,
    }


def run_industrial(paths: list[str], profile: ExecutionProfile,
                   workers: int, pass_label: str) -> dict:
    rt = IndustrialRuntime(profile=profile)
    rt.bootstrap(K3FilePlugin(), K3DedupPlugin())

    t0 = time.perf_counter()
    refs = rt.scan_files(paths, workers=workers)
    dedup_ref = rt.run_dedup_report(refs)
    elapsed = time.perf_counter() - t0
    stats = rt.report()

    dup_groups = 0
    recoverable = 0
    if dedup_ref:
        dup_groups = dedup_ref.metadata.get("duplicate_groups", 0)
        recoverable = dedup_ref.metadata.get("bytes_recoverable", 0)

    out_dir = Path(__file__).parent / "output"
    rt.telemetry.export_csv(out_dir / f"telemetria_{pass_label}.csv")
    rt.telemetry.export_json(out_dir / f"informe_{pass_label}.json", extra=stats)

    return {
        "name": f"Industrial ({pass_label})",
        "files": len(refs),
        "elapsed_s": elapsed,
        "alu_executions": stats["executed"],
        "duplicate_groups": dup_groups,
        "bytes_recoverable": recoverable,
        "knowledge_hits": stats["knowledge_hits"],
        "plugin_stats": stats.get("plugin_stats", {}),
    }


def print_result(r: dict, baseline_s: float | None = None) -> None:
    print(f"  {r['name']}")
    print(f"    Ficheros          : {r['files']}")
    print(f"    Tiempo            : {r['elapsed_s']:.4f} s")
    print(f"    Ejecuciones ALU   : {r['alu_executions']}")
    print(f"    Knowledge hits    : {r['knowledge_hits']}")
    print(f"    Grupos duplicados : {r['duplicate_groups']}")
    if r.get("bytes_recoverable"):
        print(f"    Bytes recuperables: {r['bytes_recoverable']:,}")
    if baseline_s and r["elapsed_s"]:
        saving = (1 - r["elapsed_s"] / baseline_s) * 100
        print(f"    Ahorro tiempo     : {saving:.1f}%")
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description="AntiPC demo industrial")
    parser.add_argument("--root", type=str,
                        default=str(Path(__file__).parent.parent / "sale-it"))
    parser.add_argument("--max-files", type=int, default=80)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    root = Path(args.root)
    if not root.is_dir():
        print(f"  ERROR: no existe {root}")
        return 1

    paths = collect_files(root, args.max_files)
    if not paths:
        print("  ERROR: sin ficheros")
        return 1

    print()
    print("=" * 72)
    print("  AntiPC — DEMO INDUSTRIAL")
    print("  Caso: auditoría de integridad + deduplicación (k3dedup)")
    print(f"  Hash backend: {get_backend()}")
    print(f"  Carpeta: {root}  ({len(paths)} ficheros)")
    print("=" * 72)
    print()

    mono = run_monolithic(paths)
    print_result(mono)
    baseline = mono["elapsed_s"]

    profile = ExecutionProfile.cluster()
    rt = IndustrialRuntime(profile=profile)
    rt.bootstrap(K3FilePlugin(), K3DedupPlugin())

    t0 = time.perf_counter()
    refs1 = rt.scan_files(paths, workers=args.workers)
    dedup1 = rt.run_dedup_report(refs1)
    stats1 = rt.report()
    pass1 = {
        "name": "Industrial (pasada1 — construye K(N))",
        "files": len(refs1),
        "elapsed_s": time.perf_counter() - t0,
        "alu_executions": stats1["executed"],
        "duplicate_groups": dedup1.metadata.get("duplicate_groups", 0) if dedup1 else 0,
        "knowledge_hits": stats1["knowledge_hits"],
    }
    out = Path(__file__).parent / "output"
    rt.telemetry.export_csv(out / "telemetria_pasada1.csv")
    print_result(pass1, baseline)

    exec_before = rt.report()["executed"]
    t0 = time.perf_counter()
    refs2 = rt.scan_files(paths, workers=args.workers)
    dedup2 = rt.run_dedup_report(refs2)
    stats2 = rt.report()
    pass2 = {
        "name": "Industrial (pasada2 — reutiliza K(N))",
        "files": len(refs2),
        "elapsed_s": time.perf_counter() - t0,
        "alu_executions": stats2["executed"] - exec_before,
        "knowledge_hits": stats2["knowledge_hits"],
        "duplicate_groups": dedup2.metadata.get("duplicate_groups", 0) if dedup2 else 0,
        "bytes_recoverable": dedup2.metadata.get("bytes_recoverable", 0) if dedup2 else 0,
    }
    rt.telemetry.export_csv(out / "telemetria_pasada2.csv")
    rt.telemetry.export_json(out / "informe_industrial.json", extra=stats2)
    print_result(pass2, baseline)

    alu_saved = 0
    if pass1["alu_executions"]:
        hits_pass2 = pass2["knowledge_hits"] - pass1["knowledge_hits"]
        alu_saved = (hits_pass2 / max(len(paths), 1)) * 100

    print("-" * 72)
    print("  CONCLUSIÓN INDUSTRIAL")
    print("-" * 72)
    print(f"  Reutilización pasada 2: {alu_saved:.1f}% ficheros sin recalcular (K(N))")
    print(f"  Total knowledge hits: {pass2['knowledge_hits']} (objetivo del chat)")
    print(f"  Telemetría CSV : antipc/output/telemetria_pasada*.csv")
    print(f"  Informe JSON   : antipc/output/informe_pasada*.json")
    print()
    print("  Ley AntiPC: P_util(N) = N·E(N) + K(N)")
    print("  La 2ª auditoría del mismo lote reutiliza conocimiento → menos ALU.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())