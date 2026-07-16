"""
Criterio dos primos VMA — L(n), m(n), I(n) (Libro 2 / teorema 01).

Filtrado desde anexoE_L_m_script.py + tablas L_m v3 (OneDrive pack 2025-09).
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

E_MINUS_2 = math.e - 2.0

_REPO = Path(__file__).resolve().parents[4]
_DATOS = _REPO / "mrauv" / "datos"
_DATASET = _DATOS / "Lm_dataset_v3.csv"
_TABLA = _DATOS / "Tabla_L_m_intervalos_v3.csv"

_FAC: list[float] = []


def _ensure_fac() -> list[float]:
    global _FAC
    if not _FAC:
        fac = [1.0]
        for i in range(1, 201):
            fac.append(fac[-1] / i)
        _FAC = fac
    return _FAC


def L(n: int) -> int:
    return math.isqrt(n + 3) + 7


def K(n: int) -> int:
    return max(2, int(math.isqrt(n) * 9 // 24))


def m_sum(n: int) -> float:
    """m(n) = Σ_{i=2}^{K(n)} √(n+3)/i!  (anexoE, inversos factorial)."""
    fac = _ensure_fac()
    k_max = min(max(2, K(n)), 200)
    root = math.sqrt(n + 3)
    return sum(root * fac[i] for i in range(2, k_max + 1))


def m_approx(n: int) -> float:
    return E_MINUS_2 * math.sqrt(n + 3)


def interval_I(n: int) -> tuple[int, int]:
    """I(n) = [n − ⌊√(n+3)⌋ − 3, n + 3]."""
    r = math.isqrt(n + 3)
    return n - r - 3, n + 3


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


def count_primes_in_interval(n: int) -> tuple[int, list[int]]:
    lo, hi = interval_I(n)
    lo = max(lo, 2)
    is_prime = _sieve(hi + 2)
    hi = min(hi, len(is_prime) - 1)
    primos = [p for p in range(lo, hi + 1) if is_prime[p]]
    return len(primos), primos


@dataclass
class DosPrimosReport:
    n: int
    L_val: int
    K_val: int
    m_val: float
    m_aprox: float
    margin: float
    interval: tuple[int, int]
    prime_count: int
    primes_sample: list[int]
    criterio_ok: bool
    interval_ok: bool
    csv_margin: float | None = None
    csv_match: bool | None = None


def verify_n(n: int, *, check_csv: bool = True, tol: float = 1e-4) -> DosPrimosReport:
    l_val = L(n)
    k_val = K(n)
    m_val = m_sum(n)
    margin = l_val - m_val
    interval = interval_I(n)
    count, primos = count_primes_in_interval(n)

    csv_margin = None
    csv_match = None
    if check_csv and _DATASET.is_file():
        for row in load_dataset():
            if row.n == n:
                csv_margin = row.l_minus_m
                csv_match = abs(margin - csv_margin) <= tol
                break

    return DosPrimosReport(
        n=n,
        L_val=l_val,
        K_val=k_val,
        m_val=m_val,
        m_aprox=m_approx(n),
        margin=margin,
        interval=interval,
        prime_count=count,
        primes_sample=primos[:10],
        criterio_ok=margin >= 2,
        interval_ok=count >= 2,
        csv_margin=csv_margin,
        csv_match=csv_match,
    )


@dataclass
class DatasetRow:
    n: int
    L_val: int
    K_val: int
    m_val: float
    l_minus_m: float


@dataclass
class TablaRow:
    n: int
    sqrt_n3: float
    L_val: int
    K_val: int
    m_sum: float
    m_aprox: float
    l_minus_m_sum: float
    l_minus_m_aprox: float
    ge_2: bool
    interval: str


def load_dataset() -> list[DatasetRow]:
    if not _DATASET.is_file():
        return []
    rows: list[DatasetRow] = []
    with _DATASET.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            rows.append(
                DatasetRow(
                    n=int(raw["n"]),
                    L_val=int(raw["L(n)"]),
                    K_val=int(raw["K"]),
                    m_val=float(raw["m_sum(n)"]),
                    l_minus_m=float(raw["L_minus_m"]),
                )
            )
    return rows


def load_tabla() -> list[TablaRow]:
    if not _TABLA.is_file():
        return []
    rows: list[TablaRow] = []
    with _TABLA.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            rows.append(
                TablaRow(
                    n=int(raw["n"]),
                    sqrt_n3=float(raw["sqrt(n+3)"]),
                    L_val=int(raw["L(n)"]),
                    K_val=int(raw["K"]),
                    m_sum=float(raw["m_sum"]),
                    m_aprox=float(raw["m_aprox"]),
                    l_minus_m_sum=float(raw["L_minus_m_sum"]),
                    l_minus_m_aprox=float(raw["L_minus_m_aprox"]),
                    ge_2=raw["ge_2_sum"].strip().lower().startswith("s"),
                    interval=raw["I(n)"],
                )
            )
    return rows


def audit_dataset(tol: float = 1e-4) -> tuple[bool, list[str]]:
    errors: list[str] = []
    for row in load_dataset():
        m_calc = m_sum(row.n)
        margin = L(row.n) - m_calc
        if abs(m_calc - row.m_val) > tol:
            errors.append(f"n={row.n}: m_calc={m_calc} != csv={row.m_val}")
        if abs(margin - row.l_minus_m) > tol:
            errors.append(f"n={row.n}: L-m={margin} != csv={row.l_minus_m}")
    return len(errors) == 0, errors


def audit_tabla_primes() -> tuple[bool, list[str]]:
    errors: list[str] = []
    for row in load_tabla():
        count, _ = count_primes_in_interval(row.n)
        if row.ge_2 and count < 2:
            errors.append(f"n={row.n}: CSV ge_2 pero solo {count} primos en I(n)")
        if not row.ge_2 and count >= 2:
            errors.append(f"n={row.n}: CSV no ge_2 pero {count} primos en I(n)")
    return len(errors) == 0, errors


def format_verify(r: DosPrimosReport) -> str:
    lo, hi = r.interval
    primos = str(r.primes_sample) + ("..." if r.prime_count > len(r.primes_sample) else "")
    lines = [
        f"Dos primos VMA — n={r.n}",
        f"  L(n)      : {r.L_val}",
        f"  K(n)      : {r.K_val}",
        f"  m(n)      : {r.m_val:.6f}  (aprox (e−2)√n+3 = {r.m_aprox:.6f})",
        f"  L−m       : {r.margin:.6f}  (criterio ≥2: {'SÍ' if r.criterio_ok else 'NO'})",
        f"  I(n)      : [{lo}, {hi}]",
        f"  #primos   : {r.prime_count}  {primos}",
        f"  ≥2 primos : {'SÍ' if r.interval_ok else 'NO'}",
    ]
    if r.csv_margin is not None:
        match = "SÍ" if r.csv_match else "NO"
        lines.append(f"  CSV L−m   : {r.csv_margin:.6f}  (coincide: {match})")
    return "\n".join(lines)


def format_tabla(rows: list[TablaRow]) -> str:
    lines = [
        "Tabla L_m intervalos v3 (anexoE / OneDrive)",
        f"  {'n':>10}  {'L':>5}  {'K':>4}  {'L−m':>10}  {'≥2':>4}  I(n)",
        "-" * 62,
    ]
    for row in rows:
        lines.append(
            f"  {row.n:>10}  {row.L_val:>5}  {row.K_val:>4}  "
            f"{row.l_minus_m_sum:>10.3f}  {'SÍ' if row.ge_2 else 'no':>4}  {row.interval}"
        )
    return "\n".join(lines)


def format_audit_dataset(ok: bool, errors: list[str]) -> str:
    if ok:
        return f"Auditoría Lm_dataset_v3.csv: OK ({len(load_dataset())} filas)"
    lines = [f"Auditoría Lm_dataset_v3.csv: FALLO ({len(errors)} errores)"]
    lines.extend(f"  {e}" for e in errors[:8])
    return "\n".join(lines)


def format_audit_tabla(ok: bool, errors: list[str]) -> str:
    if ok:
        return f"Auditoría Tabla_L_m + I(n): OK ({len(load_tabla())} filas)"
    lines = [f"Auditoría Tabla_L_m + I(n): FALLO ({len(errors)} errores)"]
    lines.extend(f"  {e}" for e in errors[:8])
    return "\n".join(lines)