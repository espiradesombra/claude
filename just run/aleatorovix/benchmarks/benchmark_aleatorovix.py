#!/usr/bin/env python3
"""
Benchmark Aleatorovix v1.0 — cribas, MDC, organismo.
Genera benchmarks/results.csv y benchmarks/RESULTADOS_BENCHMARK.txt
"""
from __future__ import annotations

import csv
import math
import random
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from aleatorovix import AleatorovixOrganismo
from nucleo.criba import criba_desmemoriada, criba_hibrida, es_primo
from nucleo.mdc_memoria import AleatorovixMemory, diente_frac

OUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUT_DIR / "results.csv"
TXT_PATH = OUT_DIR / "RESULTADOS_BENCHMARK.txt"

def _mk_semiprimes() -> list[tuple[int, int, int]]:
    pairs = [(11, 13), (29, 31), (59, 61), (101, 103), (311, 329), (1009, 1013)]
    out = []
    for a, b in pairs:
        n = a * b
        out.append((n, a, b))
    return out


SEMIPRIMOS_FULL = _mk_semiprimes()


@dataclass
class BenchRow:
    suite: str
    caso: str
    metrica: str
    valor: float
    unidad: str
    notas: str = ""


def _timer(fn, *args, **kwargs) -> tuple[float, object]:
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    ms = (time.perf_counter() - t0) * 1000.0
    return ms, result


def bench_cribas(rows: list[BenchRow]) -> None:
    limites = [1_000, 5_000, 10_000, 50_000]
    for lim in limites:
        ms_d, primos_d = _timer(criba_desmemoriada, lim)
        ms_h, primos_h = _timer(criba_hibrida, lim)
        rows.append(BenchRow("criba", f"limite_{lim}", "tiempo_desmemoriada_ms", ms_d, "ms"))
        rows.append(BenchRow("criba", f"limite_{lim}", "tiempo_hibrida_ms", ms_h, "ms"))
        rows.append(BenchRow("criba", f"limite_{lim}", "primos_desmemoriada", len(primos_d), "count"))
        rows.append(BenchRow("criba", f"limite_{lim}", "primos_hibrida", len(primos_h), "count"))
        ratio = ms_d / ms_h if ms_h > 0 else 0
        rows.append(
            BenchRow(
                "criba",
                f"limite_{lim}",
                "ratio_desmem_vs_hibrida",
                ratio,
                "x",
                f"desmemoriada {'mas rapida' if ratio < 1 else 'mas lenta'}",
            )
        )


def mdc_lineal(n: int, max_m: int | None = None) -> tuple[int, float, int]:
    """Escaneo lineal hasta max_m — cuenta evaluaciones d(m)."""
    if max_m is None:
        max_m = math.isqrt(n) + 1
    mejor_m, mejor_d, evals = 1, 1.0, 0
    for m in range(1, max_m + 1):
        d = diente_frac(n, m)
        evals += 1
        if abs(d - 0.5) < abs(mejor_d - 0.5):
            mejor_m, mejor_d = m, d
    return mejor_m, mejor_d, evals


def mdc_adaptativo(n: int, pasos: int | None = None) -> tuple[int, float, int]:
    if pasos is None:
        pasos = min(24, max(6, math.isqrt(n) // 4))
    mem = AleatorovixMemory()
    tray = mem.explorar_mdc(n, pasos=pasos)
    mejor = min(tray, key=lambda t: abs(t[1] - 0.5))
    return mejor[0], mejor[1], len(tray)


def bench_mdc(rows: list[BenchRow]) -> None:
    for n, fa, fb in SEMIPRIMOS_FULL:
        caso = f"N_{n}"
        ms_lin, (m_lin, d_lin, evals_lin) = _timer(mdc_lineal, n)
        ms_adp, (m_adp, d_adp, evals_adp) = _timer(mdc_adaptativo, n)

        rows.append(BenchRow("mdc", caso, "evals_lineal", evals_lin, "count", f"factores {fa}x{fb}"))
        rows.append(BenchRow("mdc", caso, "evals_adaptativo", evals_adp, "count"))
        rows.append(BenchRow("mdc", caso, "tiempo_lineal_ms", ms_lin, "ms"))
        rows.append(BenchRow("mdc", caso, "tiempo_adaptativo_ms", ms_adp, "ms"))
        rows.append(BenchRow("mdc", caso, "error_lineal_abs_d-0.5", abs(d_lin - 0.5), "frac"))
        rows.append(BenchRow("mdc", caso, "error_adaptativo_abs_d-0.5", abs(d_adp - 0.5), "frac"))
        ahorro = (1.0 - evals_adp / evals_lin) * 100 if evals_lin else 0
        rows.append(BenchRow("mdc", caso, "ahorro_evals_pct", ahorro, "%"))


def bench_organismo(rows: list[BenchRow], ciclos: int = 200) -> None:
    for ping in (False, True):
        label = "con_ping" if ping else "sin_ping"
        org = AleatorovixOrganismo(usar_ping=ping)
        acciones_salto = 0
        t0 = time.perf_counter()
        for _ in range(ciclos):
            r = org.ciclo(verbose=False)
            if r["decision"] in (1, 2):
                acciones_salto += 1
        total_ms = (time.perf_counter() - t0) * 1000.0
        primos = len(set(org.primos_encontrados))
        rows.append(BenchRow("organismo", label, "ciclos", ciclos, "count"))
        rows.append(BenchRow("organismo", label, "tiempo_total_ms", total_ms, "ms"))
        rows.append(BenchRow("organismo", label, "ms_por_ciclo", total_ms / ciclos, "ms"))
        rows.append(BenchRow("organismo", label, "saltos_6k", acciones_salto, "count"))
        rows.append(BenchRow("organismo", label, "primos_unicos", primos, "count"))
        if acciones_salto:
            rows.append(
                BenchRow(
                    "organismo",
                    label,
                    "tasa_primo_en_salto_pct",
                    100.0 * primos / acciones_salto,
                    "%",
                )
            )


def _entropia_bits(contadores: Counter[int]) -> float:
    total = sum(contadores.values())
    if total == 0:
        return 0.0
    h = 0.0
    for c in contadores.values():
        p = c / total
        if p > 0:
            h -= p * math.log2(p)
    return h


def bench_azar_masivo(rows: list[BenchRow], ciclos: int = 10_000) -> None:
    """
    Sesión larga RTL: muchas ejecuciones acumulan azar físico masivo.
    No comparar con 300 intentos PRNG — comparar cobertura y entropía a igual presupuesto.
    """
    for direccion in ("rtl", "ltr"):
        org = AleatorovixOrganismo(usar_ping=False, direccion=direccion)
        t0 = time.perf_counter()
        for _ in range(ciclos):
            org.ciclo(verbose=False)
        ms = (time.perf_counter() - t0) * 1000.0

        primos = len(set(org.primos_encontrados))
        candidatos = len(org.candidatos_vistos)
        ks = len(org.k_vistos)
        saltos = sum(1 for d in org.decisiones if d in (1, 2))
        ent_k = _entropia_bits(Counter(org.k_historial))
        ent_dec = _entropia_bits(Counter(org.decisiones))
        cobertura_k_pct = 100.0 * ks / 10_000

        label = f"ciclos_{ciclos}_{direccion}"
        rows.append(BenchRow("azar_masivo", label, "ciclos", ciclos, "count"))
        rows.append(BenchRow("azar_masivo", label, "tiempo_ms", ms, "ms"))
        rows.append(BenchRow("azar_masivo", label, "primos_unicos", primos, "count"))
        rows.append(BenchRow("azar_masivo", label, "candidatos_unicos", candidatos, "count"))
        rows.append(BenchRow("azar_masivo", label, "k_unicos", ks, "count"))
        rows.append(BenchRow("azar_masivo", label, "cobertura_k_pct", cobertura_k_pct, "%"))
        rows.append(BenchRow("azar_masivo", label, "saltos_6k", saltos, "count"))
        rows.append(BenchRow("azar_masivo", label, "entropia_k_bits", ent_k, "bits"))
        rows.append(BenchRow("azar_masivo", label, "entropia_decision_bits", ent_dec, "bits"))

    # PRNG mismo presupuesto: 10000 k aleatorios
    random.seed(42)
    primos_prng: set[int] = set()
    ks_prng: set[int] = set()
    t1 = time.perf_counter()
    for _ in range(ciclos):
        k = random.randint(0, 9_999)
        ks_prng.add(k)
        for cand in (6 * k + 1, 6 * k - 1):
            if cand > 1 and es_primo(cand):
                primos_prng.add(cand)
    ms_prng = (time.perf_counter() - t1) * 1000.0

    label = f"ciclos_{ciclos}_prng"
    rows.append(BenchRow("azar_masivo", label, "ciclos", ciclos, "count", "random.randint"))
    rows.append(BenchRow("azar_masivo", label, "tiempo_ms", ms_prng, "ms"))
    rows.append(BenchRow("azar_masivo", label, "primos_unicos", len(primos_prng), "count"))
    rows.append(BenchRow("azar_masivo", label, "k_unicos", len(ks_prng), "count"))
    rows.append(BenchRow("azar_masivo", label, "cobertura_k_pct", 100.0 * len(ks_prng) / 10_000, "%"))
    ks_prng_hist = []
    random.seed(42)
    for _ in range(ciclos):
        ks_prng_hist.append(random.randint(0, 9_999))
    rows.append(BenchRow("azar_masivo", label, "entropia_k_bits", _entropia_bits(Counter(ks_prng_hist)), "bits"))


def bench_sintonia_vs_azar(rows: list[BenchRow], intentos: int = 500) -> None:
    """Compara primos hallados por organismo vs k aleatorio 6k±1."""
    org = AleatorovixOrganismo(usar_ping=False)
    t0 = time.perf_counter()
    org.respirar(n=intentos, verbose=False)
    ms_org = (time.perf_counter() - t0) * 1000.0
    primos_org = len(set(org.primos_encontrados))

    random.seed(42)
    primos_azar: set[int] = set()
    t1 = time.perf_counter()
    for _ in range(intentos):
        k = random.randint(1, 10_000)
        for cand in (6 * k + 1, 6 * k - 1):
            if cand > 1 and es_primo(cand):
                primos_azar.add(cand)
    ms_azar = (time.perf_counter() - t1) * 1000.0

    rows.append(BenchRow("sintonia", f"intentos_{intentos}", "primos_organismo", primos_org, "count"))
    rows.append(BenchRow("sintonia", f"intentos_{intentos}", "primos_azar", len(primos_azar), "count"))
    rows.append(BenchRow("sintonia", f"intentos_{intentos}", "tiempo_organismo_ms", ms_org, "ms"))
    rows.append(BenchRow("sintonia", f"intentos_{intentos}", "tiempo_azar_ms", ms_azar, "ms"))


def escribir_csv(rows: list[BenchRow]) -> None:
    with CSV_PATH.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["suite", "caso", "metrica", "valor", "unidad", "notas"])
        for r in rows:
            w.writerow([r.suite, r.caso, r.metrica, f"{r.valor:.6f}", r.unidad, r.notas])


def escribir_txt(rows: list[BenchRow]) -> None:
    lines = [
        "=" * 72,
        "ALEATOROVIX v1.0 — RESULTADOS BENCHMARK",
        f"Fecha: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"Maquina: Python {sys.version.split()[0]} · Windows",
        "=" * 72,
        "",
    ]

    by_suite: dict[str, list[BenchRow]] = {}
    for r in rows:
        by_suite.setdefault(r.suite, []).append(r)

    for suite, srows in by_suite.items():
        lines.append(f"## {suite.upper()}")
        lines.append("-" * 40)
        current_caso = None
        for r in srows:
            if r.caso != current_caso:
                current_caso = r.caso
                lines.append(f"  [{r.caso}]")
            lines.append(f"    {r.metrica}: {r.valor:.4f} {r.unidad}  {r.notas}")
        lines.append("")

    lines.extend(
        [
            "=" * 72,
            "INTERPRETACION",
            "=" * 72,
            "",
            "CRIBA: compara criba_desmemoriada (wheel 6k) vs criba_hibrida (saltos 2p/4p).",
            "MDC:   evals_adaptativo << evals_lineal indica memoria v6 reduce trabajo.",
            "ORGANISMO: ms_por_ciclo mide overhead entropia+mascara+desmemoria.",
            "AZAR_MASIVO: muchas ejecuciones RTL acumulan entropia fisica.",
            "  No es PRNG: cada ciclo roba nanos+pila+basura; k=9999-(nanos%%10000).",
            "  La sesion larga ES el metodo — comparar cobertura_k y entropia_k a 10k ciclos.",
            "SINTONIA (300): solo referencia corta; NO invalida azar masivo acumulado.",
            "",
            "Regenerar:",
            "  python benchmarks/benchmark_aleatorovix.py",
            "",
        ]
    )
    TXT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    print("[*] Benchmark Aleatorovix v1.0")
    rows: list[BenchRow] = []

    print("  [1/4] Cribas...")
    bench_cribas(rows)

    print("  [2/4] MDC lineal vs adaptativo...")
    bench_mdc(rows)

    print("  [3/4] Organismo (ciclos)...")
    bench_organismo(rows, ciclos=100)

    print("  [4/5] Azar masivo (sesion larga RTL)...")
    bench_azar_masivo(rows, ciclos=10_000)

    print("  [5/5] Sintonia vs azar (corto — referencia)...")
    bench_sintonia_vs_azar(rows, intentos=300)

    escribir_csv(rows)
    escribir_txt(rows)
    print(f"\nOK: {len(rows)} metricas")
    print(f"  CSV: {CSV_PATH}")
    print(f"  TXT: {TXT_PATH}")

    # Resumen rapido en consola
    print("\n--- RESUMEN ---")
    for r in rows:
        if r.metrica in (
            "ahorro_evals_pct",
            "ms_por_ciclo",
            "primos_unicos",
            "ratio_desmem_vs_hibrida",
            "cobertura_k_pct",
            "entropia_k_bits",
        ):
            print(f"  {r.suite}/{r.caso} {r.metrica} = {r.valor:.2f} {r.unidad}")


if __name__ == "__main__":
    main()