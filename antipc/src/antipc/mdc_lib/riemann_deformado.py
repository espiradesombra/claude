"""
Riemann deformado VMA — R̂, R̃ y densidad por capas D₀=4/9 (Libro 3).

Filtrado desde filesclaude 6-5/riemann_deformado.py.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

D0 = 4 / 9
EULER_GAMMA = 0.5772156649015328606


def mobius(n: int) -> int:
    if n == 1:
        return 1
    factors = 0
    d = 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            factors += 1
            if n % d == 0:
                return 0
        d += 1
    if n > 1:
        factors += 1
    return (-1) ** factors


def Li(x: float, terms: int = 50) -> float:
    if x <= 0:
        return float("nan")
    if x == 1:
        return float("-inf")
    t = math.log(x)
    ei = EULER_GAMMA + math.log(abs(t)) if t != 0 else float("-inf")
    term = 1.0
    factorial = 1.0
    for k in range(1, terms + 1):
        factorial *= k
        term = t**k / (k * factorial)
        ei += term
        if abs(term) < 1e-15 * abs(ei):
            break
    return ei


def Li_safe(x: float) -> float:
    if x > 1:
        return Li(x)
    if x > 0:
        return Li(x)
    if x < 0:
        return -Li(-x) if -x > 1 else 0.0
    return 0.0


def R_clasico(n: float, k_max: int = 50) -> float:
    total = 0.0
    for k in range(1, k_max + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        total += (mu / k) * Li_safe(n ** (1.0 / k))
    return total


def R_hat(n: float, k_max: int = 50) -> float:
    total = 0.0
    for k in range(1, k_max + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        total += Li_safe(mu * (n ** (1.0 / k)))
    return total


def R_tilde(n: float, k_max: int = 50) -> float:
    total = 0.0
    for k in range(1, k_max + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        total += Li_safe(mu * (n ** (1.0 / (k + 1))))
    return total


def sigma_transition(x: float, center: float, s: float = 1e4) -> float:
    z = (x - center) / s
    z = max(-500, min(500, z))
    return 1.0 / (1.0 + math.exp(z))


def pi_density_layers(x: float, s: float | None = None) -> float:
    if x < 2:
        return 0.0
    if s is None:
        s = max(x / 10, 1.0)
    total = 0.0
    n = 0
    while 10**n <= x:
        layer = D0 / (2**n)
        center = 10.0**n
        total += layer * sigma_transition(x, center, s)
        n += 1
    return x * total


def pi_sieve(n: int) -> int:
    if n < 2:
        return 0
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0] = sieve[1] = 0
    p = 2
    while p * p <= n:
        if sieve[p]:
            sieve[p * p : n + 1 : p] = b"\x00" * len(sieve[p * p : n + 1 : p])
        p += 1
    return sum(sieve)


@dataclass
class RiemannCompareRow:
    n: int
    pi_n: int
    r_classic: float
    r_hat: float
    r_tilde: float
    layers: float


def compare_at(n: int, k_max: int = 50) -> RiemannCompareRow:
    return RiemannCompareRow(
        n=n,
        pi_n=pi_sieve(n),
        r_classic=R_clasico(n, k_max),
        r_hat=R_hat(n, k_max),
        r_tilde=R_tilde(n, k_max),
        layers=pi_density_layers(n),
    )


def compare_table(n_values: list[int], k_max: int = 50) -> list[RiemannCompareRow]:
    return [compare_at(n, k_max) for n in n_values]


def format_compare(rows: list[RiemannCompareRow], k_max: int) -> str:
    lines = [
        f"Riemann deformado VMA — K={k_max}  (D₀=4/9={D0:.6f})",
        f"  {'n':>10}  {'π(n)':>8}  {'R cl':>10}  {'R̂':>10}  {'R̃':>10}  {'capas':>10}",
        "-" * 68,
    ]
    for row in rows:
        lines.append(
            f"  {row.n:>10}  {row.pi_n:>8}  "
            f"{row.r_classic:>10.1f}  {row.r_hat:>10.1f}  "
            f"{row.r_tilde:>10.1f}  {row.layers:>10.1f}"
        )
    return "\n".join(lines)


def format_single(row: RiemannCompareRow, k_max: int) -> str:
    return "\n".join(
        [
            f"Riemann deformado — n={row.n}  K={k_max}",
            f"  π(n) real : {row.pi_n}",
            f"  R clásico : {row.r_classic:.2f}  (Δ {row.r_classic - row.pi_n:+.2f})",
            f"  R̂ hat     : {row.r_hat:.2f}  (Δ {row.r_hat - row.pi_n:+.2f})",
            f"  R̃ tilde   : {row.r_tilde:.2f}  (Δ {row.r_tilde - row.pi_n:+.2f})",
            f"  Capas D₀  : {row.layers:.2f}  (Δ {row.layers - row.pi_n:+.2f})",
        ]
    )