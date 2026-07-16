"""
Salto máximo VMA — ventana [n−⌊√(n+3)⌋−3, n+3] con ≥2 primos.

Filtrado desde filesclaude 6-5/salto_maximo.py y teorema 03 (Libro 1).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

E_MINUS_2 = math.e - 2
DENSITY_FACTOR = 1 - E_MINUS_2


def ventana(n: int) -> tuple[int, int]:
    radio = math.isqrt(n + 3)
    return n - radio - 3, n + 3


def cota_minima_primos(n: int) -> float:
    return DENSITY_FACTOR * math.sqrt(n + 3)


def n_minimo_teorema() -> float:
    return (2 / DENSITY_FACTOR) ** 2 - 3


def _sieve(limit: int) -> list[bool]:
    if limit < 2:
        return [False] * (limit + 1)
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = False
    return is_prime


def contar_primos_en_ventana(n: int, is_prime: list[bool]) -> tuple[int, list[int]]:
    lo, hi = ventana(n)
    lo = max(lo, 2)
    hi = min(hi, len(is_prime) - 1)
    primos = [p for p in range(lo, hi + 1) if is_prime[p]]
    return len(primos), primos


def verificar_conjetura(n_start: int = 50, n_max: int = 10_000) -> tuple[bool, list[tuple[int, int, list[int]]]]:
    is_prime = _sieve(n_max + 10)
    fallos: list[tuple[int, int, list[int]]] = []
    for n in range(n_start, n_max + 1):
        count, primos = contar_primos_en_ventana(n, is_prime)
        if count < 2:
            fallos.append((n, count, primos))
    return len(fallos) == 0, fallos


def L(n: int) -> float:
    return math.isqrt(n + 3) + 7


def m_accum(n: int, k_max: int | None = None) -> float:
    if k_max is None:
        k_max = max(2, int(math.isqrt(n) * 9 // 24))
    sqrt_n3 = math.sqrt(n + 3)
    total = 0.0
    factorial = 1.0
    for i in range(1, k_max + 1):
        factorial *= i
        if i >= 2:
            total += sqrt_n3 / factorial
    return total


@dataclass
class VentanaReport:
    n: int
    lo: int
    hi: int
    cota: float
    count: int
    primos: list[int]


def ventana_report(n: int) -> VentanaReport:
    is_prime = _sieve(n + 20)
    lo, hi = ventana(n)
    count, primos = contar_primos_en_ventana(n, is_prime)
    return VentanaReport(n=n, lo=lo, hi=hi, cota=cota_minima_primos(n), count=count, primos=primos)


def format_ventana(r: VentanaReport) -> str:
    ok = "SÍ" if r.count >= 2 else "NO"
    primos = str(r.primos[:8]) + ("..." if len(r.primos) > 8 else "")
    return "\n".join(
        [
            f"Salto máximo VMA — n={r.n}",
            f"  Ventana   : [{r.lo}, {r.hi}]",
            f"  Cota dens.: {r.cota:.4f}  (teórico ≥2 primos si n>{n_minimo_teorema():.0f})",
            f"  Primos    : {r.count}  {primos}",
            f"  ≥2 primos : {ok}",
        ]
    )


def format_verify(n_start: int, n_max: int, ok_all: bool, fallos: list[tuple[int, int, list[int]]]) -> str:
    lines = [
        f"Verificación salto máximo — n ∈ [{n_start}, {n_max}]",
        f"  (1−(e−2)) = {DENSITY_FACTOR:.5f}",
        f"  n_min teórico ≈ {n_minimo_teorema():.1f}",
    ]
    if ok_all:
        lines.append(f"  Resultado : OK — siempre ≥2 primos en ventana")
    else:
        lines.append(f"  Resultado : FALLO — {len(fallos)} contraejemplos")
        for n, c, p in fallos[:5]:
            lines.append(f"    n={n}: {c} primo(s) {p}")
    return "\n".join(lines)