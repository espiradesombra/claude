#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparador / benchmark de cribas VMA + referencias clásicas.

Autor material: Víctor Manzanares Alberola (VMA)
Pack: VMA_mates_rescat_2026/01_cribas
Uso civil · monorepo espiradesombra/claude

Teoría: TEORIA_CRIBAS.md
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

# Asegura import desde la misma carpeta
sys.path.insert(0, str(Path(__file__).resolve().parent))

from cribas import (  # noqa: E402
    CribaDesmemoriada,
    CribaHibrida,
    CribaModular6k,
    _is_prime_td,
)


# ─────────────────────────────────────────────
# Referencias clásicas (para comparación justa)
# ─────────────────────────────────────────────

def eratostenes(limit: int) -> List[int]:
    """Criba de Eratóstenes clásica (boolean array)."""
    if limit < 2:
        return []
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if is_prime[p]:
            start = p * p
            step = p
            is_prime[start : limit + 1 : step] = b"\x00" * (((limit - start) // step) + 1)
    return [i for i in range(2, limit + 1) if is_prime[i]]


def trial_division_list(limit: int) -> List[int]:
    """Lista de primos por trial division (solo límites pequeños)."""
    return [i for i in range(2, limit + 1) if _is_prime_td(i)]


def wheel_6k_pm1_erat(limit: int) -> List[int]:
    """
    Eratóstenes solo sobre candidatos 6k±1 (rueda mod 6).
    Referencia intermedia entre clásico y modular VMA.
    """
    if limit < 2:
        return []
    primes = [2, 3] if limit >= 3 else ([2] if limit >= 2 else [])
    if limit < 5:
        return [p for p in primes if p <= limit]

    # Mapa índice → valor 6k±1
    # Construimos lista de candidatos y un dict valor→índice sería pesado;
    # usamos array booleano completo pero saltamos pares y múltiplos de 3 en el bucle de p.
    is_comp = bytearray(limit + 1)
    is_comp[0:2] = b"\x01\x01"
    for i in range(4, limit + 1, 2):
        is_comp[i] = 1
    for i in range(9, limit + 1, 6):
        is_comp[i] = 1

    p = 5
    skip = 2
    while p * p <= limit:
        if not is_comp[p]:
            # fase 2p/4p (misma corrección que CribaModular6k)
            j = p * p
            use_2p = p % 6 == 5
            while j <= limit:
                is_comp[j] = 1
                j += (2 * p) if use_2p else (4 * p)
                use_2p = not use_2p
        p += skip
        skip = 6 - skip
    return [i for i in range(2, limit + 1) if not is_comp[i]]


# ─────────────────────────────────────────────
# Registro de métodos
# ─────────────────────────────────────────────

@dataclass
class Metodo:
    id: str
    nombre: str
    familia: str  # VMA | clasica | referencia
    fn: Callable[[int], List[int]]
    max_limit: Optional[int] = None  # None = sin tope artificial
    nota: str = ""


def _registry() -> List[Metodo]:
    return [
        Metodo("desmemoriada", "Criba Desmemoriada", "VMA",
               lambda L: CribaDesmemoriada(L).run(),
               nota="Patrones periodo 6p + AND cíclico (implementación didáctica)"),
        Metodo("modular6k", "Criba Modular 6k±1", "VMA",
               lambda L: CribaModular6k(L).run(),
               nota="Saltos 2p/4p con fase p≡±1 (mod 6) — fix 2026-07-23"),
        Metodo("hibrida", "Criba Híbrida asc+desc", "VMA",
               lambda L: CribaHibrida(L).run(),
               nota="Ascenso hasta mid + descenso por residuos"),
        Metodo("hibrida_seg", "Criba Híbrida segmentada", "VMA",
               lambda L: CribaHibrida(L).segmented_run(seg_size=max(500, L // 20 or 500)),
               nota="Bloques independientes + primos base √L"),
        Metodo("eratostenes", "Eratóstenes clásico", "clasica",
               eratostenes,
               nota="Referencia O(n log log n) en práctica"),
        Metodo("rueda6k", "Eratóstenes rueda 6k±1", "clasica",
               wheel_6k_pm1_erat,
               nota="Solo marca en clases 6k±1 (misma fase 2p/4p)"),
        Metodo("trial", "Trial division", "referencia",
               trial_division_list,
               max_limit=50_000,
               nota="Referencia de corrección; no usarla para L grandes"),
    ]


# ─────────────────────────────────────────────
# Benchmark
# ─────────────────────────────────────────────

@dataclass
class FilaResultado:
    limit: int
    method_id: str
    nombre: str
    familia: str
    n_primes: int
    time_ms: float
    ok: bool
    pi_ref: int
    speedup_vs_erat: float  # >1 = más rápido que Eratóstenes
    nota: str


def _time_call(fn: Callable[[int], List[int]], limit: int, repeats: int) -> Tuple[List[int], float]:
    """Devuelve (resultado última corrida, mejor tiempo en ms)."""
    best = float("inf")
    primes: List[int] = []
    # warm-up
    primes = fn(limit)
    for _ in range(max(1, repeats)):
        t0 = time.perf_counter()
        primes = fn(limit)
        dt = (time.perf_counter() - t0) * 1000.0
        if dt < best:
            best = dt
    return primes, best


def run_benchmark(
    limits: Sequence[int],
    methods: Optional[Sequence[str]] = None,
    repeats: int = 3,
    include_trial: bool = True,
) -> List[FilaResultado]:
    reg = {m.id: m for m in _registry()}
    if methods:
        ids = list(methods)
    else:
        ids = [m.id for m in _registry() if m.id != "trial" or include_trial]

    filas: List[FilaResultado] = []
    for L in limits:
        ref = eratostenes(L)
        pi_ref = len(ref)
        t_erat = None
        for mid in ids:
            m = reg[mid]
            if m.max_limit is not None and L > m.max_limit:
                filas.append(FilaResultado(
                    limit=L, method_id=m.id, nombre=m.nombre, familia=m.familia,
                    n_primes=-1, time_ms=float("nan"), ok=False, pi_ref=pi_ref,
                    speedup_vs_erat=float("nan"), nota=f"omitido L>{m.max_limit}",
                ))
                continue
            try:
                primes, ms = _time_call(m.fn, L, repeats)
            except Exception as e:
                filas.append(FilaResultado(
                    limit=L, method_id=m.id, nombre=m.nombre, familia=m.familia,
                    n_primes=-1, time_ms=float("nan"), ok=False, pi_ref=pi_ref,
                    speedup_vs_erat=float("nan"), nota=f"error: {e}",
                ))
                continue
            ok = primes == ref
            if m.id == "eratostenes":
                t_erat = ms
            speed = (t_erat / ms) if (t_erat and ms > 0) else float("nan")
            # si erat aún no medido en este L, rellena después
            filas.append(FilaResultado(
                limit=L, method_id=m.id, nombre=m.nombre, familia=m.familia,
                n_primes=len(primes), time_ms=ms, ok=ok, pi_ref=pi_ref,
                speedup_vs_erat=speed, nota=m.nota if ok else "DISCREPANCIA vs Eratóstenes",
            ))
        # segunda pasada: speedup si erat se midió
        t_erat_row = next((f.time_ms for f in filas if f.limit == L and f.method_id == "eratostenes" and f.ok), None)
        if t_erat_row:
            for f in filas:
                if f.limit == L and f.ok and f.time_ms > 0:
                    f.speedup_vs_erat = t_erat_row / f.time_ms
    return filas


def print_table(filas: List[FilaResultado]) -> None:
    print()
    print("=" * 88)
    print("  COMPARADOR DE CRIBAS VMA + referencias")
    print("  Víctor Manzanares Alberola · pack 01_cribas")
    print("=" * 88)
    cur = None
    for f in filas:
        if f.limit != cur:
            cur = f.limit
            print()
            print(f"  Límite L = {cur:,}   ·   π(L) ref Eratóstenes = {f.pi_ref}")
            print(f"  {'Método':28s} {'familia':10s} {'#π':>8s} {'ms':>10s} {'vs Erat':>8s}  ok")
            print("  " + "-" * 80)
        if math.isnan(f.time_ms):
            print(f"  {f.nombre:28s} {f.familia:10s} {'—':>8s} {'—':>10s} {'—':>8s}  ✗  {f.nota}")
            continue
        vs = f"{f.speedup_vs_erat:7.2f}x" if f.speedup_vs_erat == f.speedup_vs_erat else "     — "
        mark = "✓" if f.ok else "✗"
        print(f"  {f.nombre:28s} {f.familia:10s} {f.n_primes:8d} {f.time_ms:10.3f} {vs:>8s}  {mark}")
    print()
    print("  Notas:")
    print("  · 'vs Erat' > 1  ⇒ más rápido que Eratóstenes en esta máquina")
    print("  · Corrección modular 2p/4p (fase) es crítica: p≡5(mod6) arranca +2p")
    print("  · Teoría: TEORIA_CRIBAS.md  ·  Pseudocodi: ../00_sunraman_eratostenes/")
    print("=" * 88)


def write_csv(filas: List[FilaResultado], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(filas[0]).keys()) if filas else [])
        if filas:
            w.writeheader()
            for f in filas:
                w.writerow(asdict(f))
    print(f"  CSV → {path}")


def write_report_txt(filas: List[FilaResultado], path: Path) -> None:
    lines = [
        "COMPARADOR CRIBAS VMA — informe",
        f"Generado: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
    ]
    cur = None
    for f in filas:
        if f.limit != cur:
            cur = f.limit
            lines.append(f"--- L={cur}  π_ref={f.pi_ref} ---")
        lines.append(
            f"  {f.method_id:14s}  n={f.n_primes:8d}  ms={f.time_ms:10.3f}  "
            f"ok={f.ok}  speedup_erat={f.speedup_vs_erat:.3f}  {f.nota}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  TXT → {path}")


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────

def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Comparador de cribas VMA (benchmark + validación + CSV)",
    )
    ap.add_argument(
        "--limits",
        type=str,
        default="1000,5000,20000,100000",
        help="Límites separados por coma (default: 1000,5000,20000,100000)",
    )
    ap.add_argument(
        "--methods",
        type=str,
        default="",
        help="ids separados por coma (default: todos salvo trial en L grandes)",
    )
    ap.add_argument("--repeats", type=int, default=3, help="Repeticiones; se toma el mejor ms")
    ap.add_argument("--no-trial", action="store_true", help="No incluir trial division")
    ap.add_argument(
        "--quick",
        action="store_true",
        help="Atajo: L=500,2000,10000 y 2 repeats",
    )
    ap.add_argument(
        "--csv",
        type=str,
        default="",
        help="Ruta CSV de salida (default: results/benchmark_cribas_<ts>.csv)",
    )
    ap.add_argument(
        "--report",
        type=str,
        default="",
        help="Ruta informe txt (default: results/benchmark_cribas_<ts>.txt)",
    )
    ap.add_argument("--list", action="store_true", help="Lista métodos y sale")
    args = ap.parse_args(argv)

    if args.list:
        for m in _registry():
            print(f"  {m.id:14s}  [{m.familia:10s}]  {m.nombre}")
            print(f"                 {m.nota}")
        return 0

    if args.quick:
        limits = [500, 2000, 10_000]
        repeats = 2
    else:
        limits = [int(x.strip()) for x in args.limits.split(",") if x.strip()]
        repeats = args.repeats

    methods = [x.strip() for x in args.methods.split(",") if x.strip()] or None
    filas = run_benchmark(
        limits,
        methods=methods,
        repeats=repeats,
        include_trial=not args.no_trial,
    )
    print_table(filas)

    out_dir = Path(__file__).resolve().parent / "results"
    ts = time.strftime("%Y%m%d_%H%M%S")
    csv_path = Path(args.csv) if args.csv else out_dir / f"benchmark_cribas_{ts}.csv"
    report_path = Path(args.report) if args.report else out_dir / f"benchmark_cribas_{ts}.txt"
    write_csv(filas, csv_path)
    write_report_txt(filas, report_path)

    # resumen ganadores VMA por límite (entre ok)
    print("  Resumen (más rápido con ok=True, por L):")
    for L in limits:
        cands = [f for f in filas if f.limit == L and f.ok and not math.isnan(f.time_ms)]
        if not cands:
            continue
        best = min(cands, key=lambda f: f.time_ms)
        vma = [f for f in cands if f.familia == "VMA"]
        best_vma = min(vma, key=lambda f: f.time_ms) if vma else None
        print(f"    L={L:>8,}: global={best.nombre} ({best.time_ms:.2f} ms)", end="")
        if best_vma:
            print(f"  |  mejor VMA={best_vma.nombre} ({best_vma.time_ms:.2f} ms)")
        else:
            print()
    print()
    return 0 if all(f.ok or math.isnan(f.time_ms) or "omitido" in f.nota for f in filas) else 1


if __name__ == "__main__":
    raise SystemExit(main())
