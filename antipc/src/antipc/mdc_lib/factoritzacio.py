"""Factorización MDC — solo toy-numbers (≤10 dígitos, uso didáctico)."""

from __future__ import annotations

import math

MAX_TOY_N = 9_999_999_999  # 10 dígitos
PETITS = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
E_MINUS_1_OVER_2 = 0.8591409142295227


def k_sweep_mdc(n: int, m_ini: int, m_fi: int) -> int:
    if m_ini < 1:
        m_ini = 1
    if m_fi < m_ini:
        return 0
    pos_fi = 2 * m_fi + 3
    pos_ini = 2 * m_ini + 3
    k_lo = max(1, n // pos_fi) if pos_fi else 1
    k_hi = n // pos_ini
    for k in range(k_lo, k_hi + 1):
        if k == 0:
            continue
        candidat = n // k
        if candidat < 3 or (candidat & 1) == 0:
            continue
        if n % candidat == 0 and candidat > 1 and candidat < n:
            return candidat
    return 0


def factorizar_mdc_toy(n: int) -> int:
    try:
        from native_engine import is_full_native, mdc_factor as native_mdc

        if is_full_native():
            return native_mdc(n)
    except Exception:
        pass

    if n > MAX_TOY_N:
        raise ValueError(f"MDC toy limit: N must be <= {MAX_TOY_N}")
    if n < 4:
        return 0
    if (n & 1) == 0:
        return 2
    for p in PETITS:
        if n % p == 0 and p < n:
            return p
    sq = int(math.isqrt(n))
    while (sq + 1) * (sq + 1) <= n:
        sq += 1
    if sq * sq == n:
        return sq
    m_max = (sq - 3) // 2
    if m_max < 1:
        return 0
    lim_criba = min(m_max, 500_000)
    d, idx = 11, 1
    salts = [2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8,
             4, 2, 4, 2, 4, 8, 6, 4, 2, 4, 6, 2, 6, 4, 6, 2, 6, 6, 4, 2,
             4, 6, 2, 6, 4, 2, 4, 2]  # simplificación roda p210
    while d <= 2 * lim_criba + 3:
        if n % d == 0 and d < n:
            return d
        d += salts[idx % len(salts)]
        idx += 1
    if lim_criba >= m_max:
        return 0
    m_conv = max(lim_criba, int(E_MINUS_1_OVER_2 * m_max))
    f = k_sweep_mdc(n, m_conv, m_max)
    if f:
        return min(f, n // f)
    f = k_sweep_mdc(n, lim_criba, m_conv)
    if f:
        return min(f, n // f)
    return 0