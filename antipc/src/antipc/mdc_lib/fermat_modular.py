"""
Fermat modular VMA — base B_n, residuo privilegiado, alineación.

Filtrado desde diamante/fermat_modular.py y teorema 15 (Libro 6).
"""

from __future__ import annotations

from dataclasses import dataclass


def fermat(n: int) -> int:
    return 2 ** (2**n) + 1


def base_b(n: int) -> int:
    prod = 1
    for k in range(n + 1):
        prod *= fermat(k)
    return 2 * prod


def residuo_privilegiado(n: int) -> int:
    mods = [fermat(k) for k in range(n + 1)]
    residues = [2] * len(mods)
    try:
        from sympy.ntheory.modular import crt

        r, _ = crt(mods, residues)
        return int(r)
    except ImportError:
        r = 2
        for m in mods:
            r = r % m
        return r


@dataclass
class FermatAlignRow:
    m: int
    aligned: bool
    fermat_mod: int
    target_mod: int


@dataclass
class FermatAlignReport:
    n: int
    base: int
    residuo: int
    rows: list[FermatAlignRow]


def alineacion_modular(n: int, span: int = 5) -> FermatAlignReport:
    bn = base_b(n)
    rn = residuo_privilegiado(n)
    rows: list[FermatAlignRow] = []
    for m in range(n + 1, n + 1 + span):
        fm = fermat(m) % bn
        aligned = fm == rn % bn
        rows.append(FermatAlignRow(m=m, aligned=aligned, fermat_mod=fm, target_mod=rn % bn))
    return FermatAlignReport(n=n, base=bn, residuo=rn, rows=rows)


def format_fermat_align(r: FermatAlignReport) -> str:
    lines = [
        f"Fermat modular VMA — n={r.n}",
        f"  B_n       : {r.base}",
        f"  r_n (CRT) : {r.residuo}",
        f"  {'m':>4}  {'F_m mod B_n':>14}  alineado",
        "  " + "-" * 32,
    ]
    for row in r.rows:
        mark = "SÍ" if row.aligned else "no"
        lines.append(f"  {row.m:>4}  {row.fermat_mod:>14}  {mark}")
    return "\n".join(lines)