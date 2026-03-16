#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mini_comparativa.py — Benchmark mínimo (precisión y tiempo) para:
  • TU_ALGO (doble ajuste simultáneo) para log_b(a) y raíz a^(1/n)
  • Referencias: math / mpmath (si disponible) / combinaciones exp/log

Uso:
  python3 mini_comparativa.py --N 200000 --seed 123 --modo both

Salidas:
  • Imprime métricas: ns/call, ULP@95%, ULPmax, errores
  • (opcional) --csv results.csv

Notas:
  • Si tienes mpmath instalado, se usa como oráculo de alta precisión.
  • Si no, se usa math como referencia (menos estricto en ULP).
"""

from __future__ import annotations
import argparse, math, random, time, sys, csv
from dataclasses import dataclass
from typing import Callable, List, Tuple

# --- Oráculo de alta precisión (si mpmath está disponible) ---
try:
    import mpmath as mp
    MP_OK = True
    mp.mp.dps = 200
except Exception:
    MP_OK = False

# --- Utilidades de ULP ---
import struct

def _to_bits(x: float) -> int:
    return struct.unpack('>q', struct.pack('>d', x))[0]

def ulp_distance(a: float, b: float) -> float:
    if math.isnan(a) or math.isnan(b):
        return float('inf')
    if math.isinf(a) or math.isinf(b):
        return 0.0 if a == b else float('inf')
    ia = _to_bits(a)
    ib = _to_bits(b)
    # Mapeo ordenado (negativos decrecen):
    if ia < 0: ia = 0x8000000000000000 - ia
    if ib < 0: ib = 0x8000000000000000 - ib
    return float(abs(ia - ib))

# --- Tu algoritmo: LOG y RAÍZ ---

def tu_log_b_a(a: float, b: float) -> float:
    """Log_b(a) por doble ajuste simultáneo (multiplicativo + exponencial).
    Parámetros optimizados (según tus docs): f1=1.2, f2=3.4, f3=-0.8, f4=1.9.
    """
    if not (a > 0.0) or not (b > 0.0) or b == 1.0:
        return float('nan')
    f1, f2, f3, f4 = 1.2, 3.4, -0.8, 1.9
    j = 1.0
    d2 = 0.0
    d1 = 0.0
    precision = 1e-15
    for _ in range(15):
        try:
            bj = b ** j
        except OverflowError:
            bj = float('inf')
        if bj != 0.0 and a / bj >= bj:
            j *= f1
            continue
        j1 = j
        d3 = d2
        d2 = d1
        # Cálculo base
        try:
            j = (a + (a / (b ** j)) - 1.0) / a
        except OverflowError:
            j = 1.0
        # Ajuste multiplicativo
        d1 = j - j1
        j *= (1.0 + f2 * abs(d1))
        # Ajuste exponencial (limitamos exponente por estabilidad)
        if abs(d3 - d1) > 1e-12:
            factor = (f3 + f4 * (d3 - d2) / (d3 - d1))
            expo = min(abs(factor), 2.0)
            try:
                j = (abs(j) ** expo)
            except OverflowError:
                j = 1e308
        if abs(j - j1) < precision:
            break
    return j


def tu_root(a: float, n: float) -> float:
    """Raíz n-ésima de a con misma filosofía: trabajar en log para estabilidad.
    Evita complejos (dominio real).
    """
    if not (n > 0.0) or not (a >= 0.0):
        return float('nan')
    if a == 0.0:
        return 0.0
    f2, f3, f4 = 3.4, -0.8, 1.9
    x = max(math.exp(math.log(a) / n), 1e-300)
    precision = 1e-15
    for _ in range(20):
        # error en dominio logarítmico: n*log(x) - log(a)
        ex = n * math.log(max(x, 1e-300))
        err = ex - math.log(a)
        x1 = x
        # corrección base tipo multiplicativa exponenciada
        corr = -err / n if err != 0.0 else 0.0
        try:
            x *= math.exp(corr)
        except OverflowError:
            x = 1e308
        d = x - x1
        # ajuste multiplicativo proporcional al desplazamiento
        x *= (1.0 + f2 * abs(d) / max(abs(x1), 1e-300))
        # ajuste exponencial suave
        trend = (d / x) if abs(x) > 0 else 0.0
        trend = max(min(trend, 1.0), -1.0)
        expo = min(abs(f3 + f4 * trend), 2.0)
        try:
            x = (abs(x) ** expo)
        except OverflowError:
            x = 1e308
        if abs(x - x1) < precision * max(1.0, abs(x1)):
            break
        x = max(x, 1e-300)
    return x

# --- Casos de prueba ---
@dataclass
class BenchResult:
    name: str
    time_ns: float
    ulp95: float
    ulpmax: float
    bad: int


def run_bench_log(N: int, seed: int) -> List[BenchResult]:
    rng = random.Random(seed)
    A = []
    B = []
    for _ in range(N):
        u = rng.random(); v = rng.random()
        # a ≈ 2^e con e en [-500, 500], b en (1, 16]
        e = -500 + u * 1000
        a = math.ldexp(1.0, int(e)) * (1.0 + 1e-12 * v)
        b = 1.0 + v * 15.0
        A.append(a); B.append(b)

    def measure(name: str, fn: Callable[[float, float], float], ref: Callable[[float, float], float]) -> BenchResult:
        ulps: List[float] = []
        bad = 0
        t0 = time.perf_counter_ns()
        for a, b in zip(A, B):
            try:
                got = fn(a, b)
            except Exception:
                got = float('nan')
            try:
                r = ref(a, b)
            except Exception:
                r = float('nan')
            u = ulp_distance(got, r)
            ulps.append(u)
            if math.isnan(got):
                bad += 1
        t1 = time.perf_counter_ns()
        ulps.sort()
        ulp95 = ulps[int(0.95 * (len(ulps) - 1))]
        ulpmax = max(ulps)
        time_ns = (t1 - t0) / len(A)
        return BenchResult(name, time_ns, ulp95, ulpmax, bad)

    # Referencias
    if MP_OK:
        def ref_logb(a, b):
            return float(mp.log(mp.mpf(a)) / mp.log(mp.mpf(b)))
    else:
        def ref_logb(a, b):
            return math.log(a) / math.log(b)

    results = []
    results.append(measure('TU_LOG', tu_log_b_a, ref_logb))
    results.append(measure('LIBM_LOG', lambda a, b: math.log(a) / math.log(b), ref_logb))
    # SLEEF/CRlibm no están disponibles en Python puro; para Docker/CPP usamos el otro paquete.
    return results


def run_bench_root(N: int, seed: int) -> List[BenchResult]:
    rng = random.Random(seed + 1)
    A = []
    Nth = []
    choices = [0.5, 2.0, 3.0, 10.0, 100.0]
    for _ in range(N):
        u = rng.random(); v = rng.random()
        e = -500 + u * 1000
        a = math.ldexp(1.0, int(e)) * (1.0 + 1e-12 * v)
        n = choices[int(v * len(choices)) % len(choices)]
        if a < 0.0 and n != 2.0:
            a = abs(a)
        A.append(a); Nth.append(n)

    def measure(name: str, fn: Callable[[float, float], float], ref: Callable[[float, float], float]) -> BenchResult:
        ulps: List[float] = []
        bad = 0
        t0 = time.perf_counter_ns()
        for a, n in zip(A, Nth):
            try:
                got = fn(a, n)
            except Exception:
                got = float('nan')
            try:
                r = ref(a, n)
            except Exception:
                r = float('nan')
            u = ulp_distance(got, r)
            ulps.append(u)
            if math.isnan(got):
                bad += 1
        t1 = time.perf_counter_ns()
        ulps.sort()
        ulp95 = ulps[int(0.95 * (len(ulps) - 1))]
        ulpmax = max(ulps)
        time_ns = (t1 - t0) / len(A)
        return BenchResult(name, time_ns, ulp95, ulpmax, bad)

    if MP_OK:
        def ref_root(a, n):
            return float(mp.power(mp.mpf(a), mp.mpf(1.0) / mp.mpf(n)))
    else:
        def ref_root(a, n):
            return math.pow(a, 1.0 / n)

    results = []
    results.append(measure('TU_ROOT', tu_root, ref_root))
    results.append(measure('LIBM_ROOT', lambda a, n: math.pow(a, 1.0 / n), ref_root))
    return results


# --- CLI ---

def main():
    p = argparse.ArgumentParser(description='Mini comparativa: TU_ALGO vs referencias (log y raíz).')
    p.add_argument('--N', type=int, default=200000, help='nº de casos aleatorios por prueba')
    p.add_argument('--seed', type=int, default=123, help='semilla RNG')
    p.add_argument('--modo', choices=['log','root','both'], default='both', help='qué probar')
    p.add_argument('--csv', type=str, default='', help='ruta CSV opcional para volcar resultados')
    args = p.parse_args()

    print(f"mpmath disponible: {MP_OK}")
    rows: List[BenchResult] = []

    if args.modo in ('log','both'):
        print("\n[LOG_b(a)] ejecutando…")
        r = run_bench_log(args.N, args.seed)
        rows.extend(r)
    if args.modo in ('root','both'):
        print("\n[ROOT a^(1/n)] ejecutando…")
        r = run_bench_root(args.N, args.seed)
        rows.extend(r)

    # Mostrar tabla
    print("\nResultados (ns/call, ULP@95%, ULPmax, bad):")
    for br in rows:
        print(f"{br.name:>10s}  {br.time_ns:10.2f}  {br.ulp95:10.3f}  {br.ulpmax:10.3f}  bad={br.bad}")

    # CSV opcional
    if args.csv:
        with open(args.csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['name','time_ns','ulp95','ulpmax','bad'])
            for br in rows:
                w.writerow([br.name, f"{br.time_ns:.2f}", f"{br.ulp95:.3f}", f"{br.ulpmax:.3f}", br.bad])
        print(f"\nCSV escrito en {args.csv}")

if __name__ == '__main__':
    main()
