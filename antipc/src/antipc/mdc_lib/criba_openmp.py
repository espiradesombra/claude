"""
Criba 6k±1 estilo anexoF (OpenMP / VMA).

Fuente canónica:
  archivos-vma/codigo/anexoF_criba6kpm1_openmp.c
  diamante/anexoF_criba6kpm1_openmp.c

Notas:
  - El .c original usa OpenMP solo para wall-time (omp_get_wtime), no para
    paralelizar el marcado. Aquí el algoritmo es secuencial en Python.
  - Corrección: tras el barrido, los candidatos no marcados (estado 1) son
    primos > √N. El anexoF solo contaba estado 2 y subcontaba π(N).

Estado: [ALG] — criba exacta sobre la rejilla 6k±1.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass


def _is_cand(v: int) -> bool:
    r = v % 6
    return r == 1 or r == 5


def _value(i: int) -> int:
    """Mapa índice → candidato 6k±1: 5,7,11,13,… (2i+3)."""
    return 2 * i + 3


@dataclass
class CribaOpenMPResult:
    limit: int
    count: int
    ms: float
    corrected_count: bool
    primes_sample: list[int]


def sieve_6k_openmp(limit: int, *, sample: int = 0) -> CribaOpenMPResult:
    """
    Criba de primos ≤ limit en la rejilla 6k±1 + {2,3}.

    Representación P[v]:
      0 = compuesto / no candidato
      1 = candidato vivo (primo si sobrevive)
      2 = primo confirmado (≤ √N, semilla del marcado)
    """
    if limit < 2:
        return CribaOpenMPResult(limit, 0, 0.0, True, [])

    t0 = time.perf_counter()
    p = bytearray(limit + 1)
    if limit >= 2:
        p[2] = 2
    if limit >= 3:
        p[3] = 2
    for v in range(5, limit + 1):
        if _is_cand(v):
            p[v] = 1

    root = int(math.isqrt(limit))
    n = 0
    while True:
        pn = _value(n)
        if pn > root or pn > limit:
            break
        if p[pn] == 0:
            n += 1
            continue
        p[pn] = 2
        m = n
        while True:
            qm = _value(m)
            if qm == 0:
                break
            if pn > limit // qm:
                break
            if qm % 3 == 0 or not _is_cand(qm):
                m += 1
                continue
            j = pn * qm
            salto1 = 2 * pn
            salto2 = 4 * pn
            toggle = True
            while j <= limit:
                if p[j] == 1:
                    p[j] = 0
                j += salto1 if toggle else salto2
                toggle = not toggle
            m += 1
        n += 1

    # Corrección anexoF: supervivientes (1) y confirmados (2) son primos.
    count = 0
    sample_list: list[int] = []
    for v in range(2, limit + 1):
        if p[v] >= 1:
            count += 1
            if sample and len(sample_list) < sample:
                sample_list.append(v)

    ms = (time.perf_counter() - t0) * 1000.0
    return CribaOpenMPResult(
        limit=limit,
        count=count,
        ms=ms,
        corrected_count=True,
        primes_sample=sample_list,
    )


def count_primes(limit: int) -> int:
    return sieve_6k_openmp(limit).count


def format_result(r: CribaOpenMPResult) -> str:
    lines = [
        f"Criba 6k±1 anexoF (Python) ≤ {r.limit}",
        f"  Primos   : {r.count}",
        f"  Conteo   : corregido (1∪2; anexoF solo contaba 2)",
        f"  Tiempo   : {r.ms:.1f} ms",
    ]
    if r.primes_sample:
        lines.append(f"  Muestra  : {r.primes_sample}")
    return "\n".join(lines)
