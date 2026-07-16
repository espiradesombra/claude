"""
Método discriminante VMA — Δ(S) = S² + 6S − (N−9) = k².

Filtrado desde filesclaude 6-5/discriminant.py y teoremas 09/11 (Libro 5).
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def delta_s(S: int, M: int) -> int:
    return S * S + 6 * S - M


def is_perfect_square(n: int) -> tuple[bool, int]:
    if n < 0:
        return False, 0
    k = math.isqrt(n)
    return k * k == n, k


def stop_condition(S: int, delta: int) -> bool:
    if delta < 0:
        return True
    sqrt_delta = math.isqrt(delta)
    return (4 * S + 16) > (2 * sqrt_delta + 1)


@dataclass
class DiscriminantResult:
    n: int
    factor1: int | None
    factor2: int | None
    steps: int
    final_s: int | None
    reason: str

    @property
    def composite(self) -> bool:
        return self.factor1 is not None


def discriminant_factor(n: int) -> DiscriminantResult:
    steps = 0
    final_s: int | None = None
    reason = ""

    for p in (2, 3, 5, 7, 11, 13):
        if n % p == 0 and n > p:
            return DiscriminantResult(n, p, n // p, 0, None, f"small factor {p}")
        if n == p:
            return DiscriminantResult(n, None, None, 0, None, "small prime")

    if n % 6 not in (1, 5):
        return DiscriminantResult(n, None, None, 0, None, "not in 6k±1 form")

    m = n - 9
    sqrt_n = math.isqrt(n)
    v0 = max(0, (sqrt_n - 3) // 2)
    s = 2 * v0

    while s >= 0:
        delta = delta_s(s, m)
        steps += 1

        if delta >= 0:
            is_sq, k = is_perfect_square(delta)
            if is_sq:
                a = s + 3 - k
                b_val = s + 3 + k
                if a > 1 and b_val > 1 and a * b_val == n:
                    return DiscriminantResult(n, a, b_val, steps, s, "discriminant square")
                v = s // 2
                b = (k - s - 3) // 2
                f1 = 2 * v + 3
                f2 = 2 * b + 3
                if f1 > 1 and f2 > 1 and f1 * f2 == n:
                    return DiscriminantResult(n, f1, f2, steps, s, "discriminant square (2v+3)")

        if stop_condition(s, delta):
            final_s = s
            reason = "deterministic stop"
            break

        s -= 2

    if not reason:
        reason = "exhausted"
    return DiscriminantResult(n, None, None, steps, final_s, reason)


@dataclass
class TrajectoryRow:
    s: int
    delta: int
    perfect: bool
    sqrt_k: str
    stop: bool


def delta_trajectory(n: int, *, show_all: bool = False) -> list[TrajectoryRow]:
    m = n - 9
    sqrt_n = math.isqrt(n)
    v0 = max(0, (sqrt_n - 3) // 2)
    s = 2 * v0
    rows: list[TrajectoryRow] = []

    while s >= 0:
        delta = delta_s(s, m)
        perfect = False
        k_str = ""
        if delta >= 0:
            perfect, k = is_perfect_square(delta)
            if perfect:
                k_str = str(k)
        stop = stop_condition(s, delta)
        if show_all or perfect or stop:
            rows.append(TrajectoryRow(s, delta, perfect, k_str, stop))
        if stop:
            break
        s -= 2

    return rows


def format_result(r: DiscriminantResult) -> str:
    lines = [
        f"Discriminante VMA — N={r.n}",
        f"  Pasos     : {r.steps}",
        f"  S final   : {r.final_s}",
        f"  Motivo    : {r.reason}",
    ]
    if r.composite:
        lines.append(f"  Factores  : {r.factor1} × {r.factor2}")
        lines.append(f"  Producto  : {r.factor1 * r.factor2}")
    else:
        lines.append("  Resultado : candidato primo (sin cuadrado perfecto en ventana)")
    return "\n".join(lines)


def format_trajectory(n: int, rows: list[TrajectoryRow]) -> str:
    lines = [
        f"Trayectoria Δ(S) — N={n}  (M=N-9={n - 9})",
        f"{'S':>6}  {'Δ(S)':>12}  {'√Δ':>8}  {'cuadrado':>9}  stop",
        "-" * 48,
    ]
    for row in rows:
        mark = ""
        if row.perfect:
            mark = "  ← FACTOR"
        elif row.stop:
            mark = "  ← STOP"
        lines.append(
            f"{row.s:>6}  {row.delta:>12}  {row.sqrt_k:>8}  "
            f"{'SÍ' if row.perfect else 'no':>9}  {'SÍ' if row.stop else 'no'}{mark}"
        )
    return "\n".join(lines)