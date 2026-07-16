"""
MRAUV — densidad cinemática de primos (Libro 2 / teoremas 14 y 23).

Extraído y unificado desde archivos históricos (filtrado 2026-07-16):
  Desktop/archivos/mrauv.py
  Desktop/archivos/anexoF_mrauv_calibrador.py
  repo/py/mrauv_goldbach.py
  graficas y explicaciones/generar_todo.py (calibrar_mrauv)

Estado: [HEUR] — modelo predictivo; no teorema analítico.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def L(n: int) -> int:
    """Longitud del corredor: ⌊√(n+3)⌋ + 7."""
    return math.isqrt(n + 3) + 7


def K(n: int) -> int:
    """Ventana factorial: ⌊⌊√n⌋ · 9/24⌋."""
    return max(2, int(math.isqrt(n) * 9 // 24))


def m_sum(n: int, k_max: int = 2000) -> float:
    """Sobreconteo m(n) vía Σ_{i=2}^{K} √(n+3)/i!."""
    k = min(K(n), k_max)
    sqrt_n3 = math.sqrt(n + 3)
    total = 0.0
    factorial = 1.0
    for i in range(1, k + 1):
        factorial *= i
        if i >= 2:
            total += sqrt_n3 / factorial
    return total


def corridor_margin(n: int) -> float:
    """Espacio libre L(n) − m(n); ≥2 ⇒ criterio dos primos (heurístico)."""
    return L(n) - m_sum(n)


def D(n: int) -> float:
    """Densidad MRAUV en [n, 2n]: (L(n) − m(n)) / (2n)."""
    if n <= 0:
        return 0.0
    return corridor_margin(n) / (2 * n)


def sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i * i : limit + 1 : i] = b"\x00" * len(is_prime[i * i : limit + 1 : i])
    return [i for i, v in enumerate(is_prime) if v]


def F_eff(n: int, small_primes: list[int]) -> float:
    """Fallo asimétrico efectivo (heurística Goldbach, ficha 14)."""
    if n <= 0:
        return 0.0
    bound = int(math.sqrt(2 * n)) + 1
    pi_2n = sum(1 for p in small_primes if p <= 2 * n)
    density_2n = pi_2n / (2 * n)
    fault = 0.0
    for p in small_primes:
        if p > bound:
            break
        fault += (2 * n // p) * density_2n
    return fault


def goldbach_margin(n: int, small_primes: list[int]) -> float:
    return D(n) - F_eff(n, small_primes) / (2 * n)


@dataclass
class CalibrarReport:
    n0: int
    dn: int
    points: list[int]
    med: list[float]
    V0: float
    a0: float
    j: float
    D0: float
    D_pred: float
    primes_est: int
    two_prime_ok: bool


def calibrar(n0: int, dn: int | None = None) -> CalibrarReport:
    """
    Calibración cinemática en 3 puntos (separación dn ≈ 2√n₀).
    V₀ = 1/med₀, a₀ = V₀/med₁, j = a₀/med₂.
    """
    if n0 < 3:
        raise ValueError("n0 debe ser >= 3")
    if dn is None:
        dn = int(2 * math.sqrt(n0))
    if dn < 1:
        dn = 1

    pts = [n0, n0 + dn, n0 + 2 * dn]
    med = [L(x) - m_sum(x) for x in pts]
    med0, med1, med2 = med

    v0 = 1.0 / med0 if med0 else 0.0
    a0 = v0 / med1 if med1 else 0.0
    j = a0 / med2 if med2 else 0.0

    d0 = 1.0 / max(1, int(math.log(max(3, n0))))
    d_pred = d0 + v0 * dn + 0.5 * a0 * (dn**2)
    primes_est = int(round(d_pred * dn))

    return CalibrarReport(
        n0=n0,
        dn=dn,
        points=pts,
        med=med,
        V0=v0,
        a0=a0,
        j=j,
        D0=d0,
        D_pred=d_pred,
        primes_est=primes_est,
        two_prime_ok=corridor_margin(n0) >= 2,
    )


@dataclass
class DensidadReport:
    n: int
    L_val: int
    m_val: float
    margin: float
    density: float
    two_prime_ok: bool


def densidad_report(n: int) -> DensidadReport:
    l_val = L(n)
    m_val = m_sum(n)
    margin = l_val - m_val
    return DensidadReport(
        n=n,
        L_val=l_val,
        m_val=m_val,
        margin=margin,
        density=D(n),
        two_prime_ok=margin >= 2,
    )


def mrauv_fit(n0: int, n1: int, n2: int) -> tuple[float, float, float]:
    d0, d1, d2 = D(n0), D(n1), D(n2)
    h = n1 - n0
    if h == 0:
        return d0, 0.0, 0.0
    v0 = (d1 - d0) / h
    a0 = (d2 - 2 * d1 + d0) / (h * h)
    return d0, v0, a0


def mrauv_predict(d0: float, v0: float, a0: float, delta_n: int) -> float:
    return d0 + v0 * delta_n + 0.5 * a0 * delta_n * delta_n


def count_goldbach(two_n: int, primes_set: set[int]) -> int:
    n = two_n // 2
    count = 0
    for p in sorted(primes_set):
        if p > n:
            break
        if (two_n - p) in primes_set:
            count += 1
    return count


@dataclass
class GoldbachScanRow:
    n: int
    density: float
    f_norm: float
    margin: float
    ok: bool


def scan_goldbach(
    n_start: int,
    n_max: int,
    delta: int,
) -> tuple[list[GoldbachScanRow], bool]:
    small = sieve_primes(int(math.sqrt(2 * n_max)) + 2)
    rows: list[GoldbachScanRow] = []
    all_ok = True
    n = n_start
    while n <= n_max:
        d = D(n)
        f_norm = F_eff(n, small) / (2 * n) if n > 0 else 0.0
        margin = d - f_norm
        ok = margin > 0
        if not ok:
            all_ok = False
        rows.append(GoldbachScanRow(n=n, density=d, f_norm=f_norm, margin=margin, ok=ok))
        n += delta
    return rows, all_ok


def format_calibrar(r: CalibrarReport) -> str:
    lines = [
        "MRAUV — calibración por tramos (3 puntos)",
        f"  n0        : {r.n0}",
        f"  dn        : {r.dn}",
        f"  puntos    : {r.points}",
        f"  med       : {[round(x, 6) for x in r.med]}",
        f"  V0        : {r.V0:.6e}",
        f"  a0        : {r.a0:.6e}",
        f"  j         : {r.j:.6e}",
        f"  D0        : {r.D0:.6f}",
        f"  D_pred    : {r.D_pred:.6f}",
        f"  primos est: {r.primes_est} en [{r.n0}, {r.n0 + r.dn}]",
        f"  2-primos  : {'OK' if r.two_prime_ok else 'no'} (L-m >= 2 en n0)",
    ]
    return "\n".join(lines)


def format_densidad(r: DensidadReport) -> str:
    return "\n".join(
        [
            f"MRAUV — densidad en n={r.n}",
            f"  L(n)      : {r.L_val}",
            f"  m(n)      : {r.m_val:.6f}",
            f"  L - m     : {r.margin:.6f}",
            f"  D(n)      : {r.density:.8f}",
            f"  2-primos  : {'OK' if r.two_prime_ok else 'no'} (heurístico)",
        ]
    )


def format_goldbach_scan(rows: list[GoldbachScanRow], all_ok: bool) -> str:
    lines = [
        f"{'n':>8}  {'D(n)':>10}  {'F/2n':>10}  {'margin':>10}  ok",
        "-" * 48,
    ]
    for row in rows:
        mark = "OK" if row.ok else "ALERT"
        lines.append(
            f"{row.n:>8}  {row.density:>10.6f}  {row.f_norm:>10.6f}  {row.margin:>+10.6f}  {mark}"
        )
    lines.append("")
    lines.append(
        "Criterio MRAUV-Goldbach satisfecho en rango"
        if all_ok
        else "Criterio falló en al menos un punto (heurístico)"
    )
    return "\n".join(lines)