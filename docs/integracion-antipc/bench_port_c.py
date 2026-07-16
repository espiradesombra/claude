#!/usr/bin/env python3
"""
Benchmark Python vs C (repo antipc_native + prototipos locales).
No mueve archivos; solo lee rutas existentes.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO_ANTIPC = Path(r"C:\Users\cuent\Desktop\repo\antipc\src\antipc")
REPO_VMA = Path(r"C:\Users\cuent\Desktop\repo\vma-methods")
DEMO_EXE = ROOT / "port_demo.exe"


def _add_paths() -> None:
    for p in (REPO_ANTIPC, REPO_VMA):
        if p.is_dir() and str(p) not in sys.path:
            sys.path.insert(0, str(p))


def bench(label: str, fn, repeats: int = 3) -> tuple[float, object]:
    best = None
    out = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = fn()
        ms = (time.perf_counter() - t0) * 1000
        if best is None or ms < best:
            best = ms
    return best or 0.0, out


def main() -> int:
    _add_paths()
    limit = 50_000
    n_mdc = 1147

    print("=" * 60)
    print("  bench_port_c — Escritorio (sin mover archivos)")
    print("=" * 60)

    # --- Cribas Python ---
    from vma_methods.cribas import CribaModular6k, CribaHibrida

    ms, c1 = bench("CribaModular6k Py", lambda: len(CribaModular6k(limit).run()))
    print(f"  CribaModular6k Python : {ms:6.1f} ms  count={c1}")

    ms, c2 = bench("CribaHibrida Py", lambda: len(CribaHibrida(limit).run()))
    print(f"  CribaHibrida Python   : {ms:6.1f} ms  count={c2}")

    # --- antipc_native DLL si existe ---
    try:
        from native_engine import (
            geo_converge,
            is_full_native,
            mdc_factor,
            sieve_count,
            status_report,
        )

        print()
        print(status_report())
        if is_full_native():
            ms, sc = bench("sieve_count C", lambda: sieve_count(limit))
            print(f"  Eratóstenes C (DLL)   : {ms:6.1f} ms  count={sc}")
            ms, f = bench("mdc_factor C", lambda: mdc_factor(1_047_029))
            print(f"  MDC factor C (DLL)    : {ms:6.1f} ms  factor={f}")
            ms, g = bench(
                "geo C",
                lambda: geo_converge(
                    "0101101011000010", [3, 5, 8, 13, 21], [6, 12, 18]
                ),
            )
            print(f"  Geo converge C (DLL)  : {ms:6.1f} ms  raw={g}")
    except Exception as e:
        print(f"  [DLL] {e}")

    # --- MDC analyze Python ---
    try:
        from mdc_lib.mdc_analyze import analyze

        ms, r = bench("mdc analyze Py", lambda: analyze(n_mdc))
        print(f"  MDC analyze Python    : {ms:6.1f} ms  union={len(r.union_both)}")
    except Exception as e:
        print(f"  [MDC] {e}")

    # --- Prototipo local port_demo.exe ---
    print()
    if DEMO_EXE.is_file():
        t0 = time.perf_counter()
        proc = subprocess.run([str(DEMO_EXE)], capture_output=True, text=True, check=False)
        ms = (time.perf_counter() - t0) * 1000
        print(f"  port_demo.exe (local C): {ms:6.1f} ms")
        print(proc.stdout)
        if proc.returncode != 0:
            print(proc.stderr)
    else:
        print("  port_demo.exe no compilado — ejecuta build_desktop.bat")

    print("=" * 60)
    print("  Ver 01_MAPA_PORT_A_C.md para siguiente oleada")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())