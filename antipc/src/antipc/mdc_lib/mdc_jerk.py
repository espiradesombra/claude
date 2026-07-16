"""
MDC-U — pinça doble 4+4 i invariante jerk (Libro 5 / DeepSeek 6 2026).

Port compacte de mdc_v17.py (filestot l5) per a CLI AntiPC.
"""

from __future__ import annotations

from fractions import Fraction
import math
from dataclasses import dataclass


def d_frac(m: int, N: int) -> Fraction:
    """Distància signada al zero de fase: R/D − ½."""
    D = 2 * (2 * m + 3)
    R = N % D
    return Fraction(R, D) - Fraction(1, 2)


def pinca_doble(m: int, N: int) -> dict:
    """
    Pinça doble 4+4: sèries endavant i enrere sobre d_frac.
    Retorna V, A, J per cada direcció i salt mestre m − V/A.
    """
    d = [d_frac(m + i, N) for i in range(4)]
    Vf = d[1] - d[0]
    Af = d[2] - 2 * d[1] + d[0]
    Jf = d[3] - 3 * d[2] + 3 * d[1] - d[0]

    dr = [d_frac(m - i, N) for i in range(4)]
    Vr = dr[1] - dr[0]
    Ar = dr[2] - 2 * dr[1] + dr[0]
    Jr = dr[3] - 3 * dr[2] + 3 * dr[1] - dr[0]

    delta_J = abs(Jf - Jr)

    salt_f = None
    if Af != 0:
        despl = Vf / Af
        despl_int = (
            int(despl + Fraction(1, 2))
            if despl >= 0
            else -int(-despl + Fraction(1, 2))
        )
        m_next = m - despl_int
        if m_next >= 1:
            salt_f = m_next

    salt_r = None
    if Ar != 0:
        despl = Vr / Ar
        despl_int = (
            int(despl + Fraction(1, 2))
            if despl >= 0
            else -int(-despl + Fraction(1, 2))
        )
        m_next = m + despl_int
        if m_next >= 1:
            salt_r = m_next

    return {
        "m": m,
        "d0": d[0],
        "Vf": Vf,
        "Af": Af,
        "Jf": Jf,
        "salt_f": salt_f,
        "Vr": Vr,
        "Ar": Ar,
        "Jr": Jr,
        "salt_r": salt_r,
        "delta_J": delta_J,
        "V_dir": Vf - abs(Vf),
    }


@dataclass
class JerkReport:
    n: int
    m: int
    pos: int
    pinza: dict
    macro_envolvente: bool
    salto_balistico: int | None
    factor_D: int | None


def _default_m(n: int) -> int:
    m = (math.isqrt(n) - 3) // 2
    return max(1, m)


def analyze_jerk(n: int, m: int | None = None, *, eps: float = 1e-12) -> JerkReport:
    """
    Analitza cinemàtica local: |J| petit → dent interior; ΔJ gran → macro-envolvente.
    Salto balístic: m* ≈ m − V/A (sèrie endavant).
    """
    if n < 4:
        raise ValueError("n debe ser >= 4")
    m0 = m if m is not None else _default_m(n)
    p = pinca_doble(m0, n)

    jf = float(p["Jf"])
    af = float(p["Af"])
    macro = float(p["delta_J"]) > eps

    salto = p["salt_f"]
    factor = None
    if salto is not None:
        D = 2 * salto + 3
        if 1 < D < n and n % D == 0:
            factor = D

    return JerkReport(
        n=n,
        m=m0,
        pos=2 * m0 + 3,
        pinza=p,
        macro_envolvente=macro,
        salto_balistico=salto,
        factor_D=factor,
    )


def format_jerk_report(r: JerkReport) -> str:
    p = r.pinza
    lines = [
        f"MDC-U jerk — n={r.n:,}  m={r.m:,}  pos=2m+3={r.pos:,}",
        f"  Macro-envolvente (ΔJ>ε) : {'sí' if r.macro_envolvente else 'no'}",
        f"  Endavant  V={p['Vf']}  A={p['Af']}  J={p['Jf']}",
        f"  Enrere    V={p['Vr']}  A={p['Ar']}  J={p['Jr']}",
        f"  ΔJ        {p['delta_J']}",
        f"  Salto balístico m−V/A     : {r.salto_balistico}",
    ]
    if r.factor_D:
        lines.append(f"  Factor trobat D=2m+3      : {r.factor_D}  (q={r.n // r.factor_D})")
    return "\n".join(lines)