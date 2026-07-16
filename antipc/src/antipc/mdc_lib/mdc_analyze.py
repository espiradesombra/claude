"""
MDC analyze — dos trenes (variación x / variación y), colisiones enteras.

Modelo:
  Trayectoria diofántica: (2x+3)(2y+3) = n  <=>  x·y = n en parametrización MDC
  (con S=2x+3, T=2y+3 factores reales).

  Tren X: x entero avanza con paso 1 (velocidad constante en x).
          En cada x, pendiente local Δy/Δx → candidato y = (n-6x-9)/(4x+6).
          Colisión si y es entero y S·T = n.

  Tren Y: y entero, simétrico con x = (n-6y-9)/(4y+6).

  Unión: pares (S,T) que coinciden en ambos sentidos → factorización confirmada.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction


def diofantic_y(n: int, x: Fraction) -> Fraction | None:
    den = 4 * x + 6
    if den == 0:
        return None
    return Fraction(n - 6 * x - 9, den)


def diofantic_x(n: int, y: Fraction) -> Fraction | None:
    den = 4 * y + 6
    if den == 0:
        return None
    return Fraction(n - 6 * y - 9, den)


def factors_from_xy(x: Fraction, y: Fraction, n: int) -> tuple[int, int] | None:
    if x.denominator != 1 or y.denominator != 1:
        return None
    s, t = 2 * int(x) + 3, 2 * int(y) + 3
    if s > 0 and t > 0 and s * t == n:
        return (min(s, t), max(s, t))
    return None


def slope_delta_y(n: int, x: int) -> Fraction | None:
    """Pendiente discreta MDC en x: Δy/Δx con Δx=1."""
    y0 = diofantic_y(n, Fraction(x))
    y1 = diofantic_y(n, Fraction(x + 1))
    if y0 is None or y1 is None:
        return None
    return y1 - y0


def slope_delta_x(n: int, y: int) -> Fraction | None:
    """Pendiente discreta MDC en y: Δx/Δy con Δy=1."""
    x0 = diofantic_x(n, Fraction(y))
    x1 = diofantic_x(n, Fraction(y + 1))
    if x0 is None or x1 is None:
        return None
    return x1 - x0


@dataclass
class Collision:
    x: int | str
    y: int | str
    s: int
    t: int
    k: int
    source: str
    slope: str = ""

    @property
    def pair(self) -> tuple[int, int]:
        return (self.s, self.t)


@dataclass
class TrainReport:
    axis: str
    steps: int
    collisions: list[Collision] = field(default_factory=list)
    avg_slope: Fraction | None = None


@dataclass
class MdcAnalyzeResult:
    n: int
    train_x: TrainReport
    train_y: TrainReport
    union_both: list[tuple[int, int]]
    union_x_only: list[tuple[int, int]]
    union_y_only: list[tuple[int, int]]
    toy_factor: int
    aux_step: int
    gcd_core: list[dict]


def _x_range(n: int) -> range:
    return range(-3, max(2, n // 2) + 1)


def _y_range(n: int) -> range:
    return range(-3, max(2, n // 2) + 1)


def _train_report_from_native(axis: str, raw) -> TrainReport:
    hits: list[Collision] = []
    for i in range(raw.n_collisions):
        h = raw.hits[i]
        hits.append(
            Collision(
                x=h.x,
                y=h.y,
                s=h.s,
                t=h.t,
                k=h.k,
                source=f"tren-{axis.upper()}",
                slope="",
            )
        )
    return TrainReport(axis=axis, steps=raw.steps, collisions=hits, avg_slope=None)


def scan_train_x(n: int) -> TrainReport:
    hits: list[Collision] = []
    slopes: list[Fraction] = []
    seen: set[tuple[int, int]] = set()

    for xi in _x_range(n):
        y = diofantic_y(n, Fraction(xi))
        if y is None:
            continue
        sl = slope_delta_y(n, xi)
        if sl is not None:
            slopes.append(sl)
        if y.denominator != 1:
            continue
        yi = int(y)
        pair = factors_from_xy(Fraction(xi), Fraction(yi), n)
        if pair is None or pair in seen:
            continue
        seen.add(pair)
        s, t = pair
        hits.append(
            Collision(
                x=xi,
                y=yi,
                s=s,
                t=t,
                k=(t - s) // 2,
                source="tren-X",
                slope=str(sl) if sl else "",
            )
        )

    avg = None
    if slopes:
        avg = sum(slopes, Fraction(0)) / len(slopes)

    return TrainReport(axis="x", steps=len(_x_range(n)), collisions=hits, avg_slope=avg)


def scan_train_y(n: int) -> TrainReport:
    hits: list[Collision] = []
    slopes: list[Fraction] = []
    seen: set[tuple[int, int]] = set()

    for yi in _y_range(n):
        x = diofantic_x(n, Fraction(yi))
        if x is None:
            continue
        sl = slope_delta_x(n, yi)
        if sl is not None:
            slopes.append(sl)
        if x.denominator != 1:
            continue
        xi = int(x)
        pair = factors_from_xy(Fraction(xi), Fraction(yi), n)
        if pair is None or pair in seen:
            continue
        seen.add(pair)
        s, t = pair
        hits.append(
            Collision(
                x=xi,
                y=yi,
                s=s,
                t=t,
                k=(t - s) // 2,
                source="tren-Y",
                slope=str(sl) if sl else "",
            )
        )

    avg = None
    if slopes:
        avg = sum(slopes, Fraction(0)) / len(slopes)

    return TrainReport(axis="y", steps=len(_y_range(n)), collisions=hits, avg_slope=avg)


def gcd_candidates(n: int, m_max: int | None = None) -> list[dict]:
    if m_max is None:
        m_max = math.isqrt(n)
    out: list[dict] = []
    for m in range(0, m_max + 1):
        s = 2 * m + 3
        g = math.gcd(n, s)
        if g > 1:
            out.append(
                {
                    "m": m,
                    "S": s,
                    "gcd": g,
                    "T": n // s if n % s == 0 else None,
                }
            )
    return out


def aux_train_step(n: int, d: int | None = None) -> int:
    """Paso constante Δx en recta auxiliar y = -(n/d)x + n (trenes en cruces enteros)."""
    if d is None:
        d = min(2 * math.isqrt(n) + 3, max(5, n // 6))
    g = math.gcd(n, d)
    return d // g if g else d


def union_collisions(
    tx: TrainReport, ty: TrainReport
) -> tuple[list[tuple[int, int]], list[tuple[int, int]], list[tuple[int, int]]]:
    px = {c.pair for c in tx.collisions}
    py = {c.pair for c in ty.collisions}
    both = sorted(px & py)
    x_only = sorted(px - py)
    y_only = sorted(py - px)
    return both, x_only, y_only


def analyze(n: int, d_aux: int | None = None) -> MdcAnalyzeResult:
    from .factoritzacio import factorizar_mdc_toy

    try:
        from native_engine import mdc_scan_trains

        tx_raw, ty_raw = mdc_scan_trains(n)
        tx = _train_report_from_native("x", tx_raw)
        ty = _train_report_from_native("y", ty_raw)
    except Exception:
        tx = scan_train_x(n)
        ty = scan_train_y(n)
    both, x_only, y_only = union_collisions(tx, ty)
    toy = 0
    try:
        toy = factorizar_mdc_toy(n)
    except ValueError:
        pass

    return MdcAnalyzeResult(
        n=n,
        train_x=tx,
        train_y=ty,
        union_both=both,
        union_x_only=x_only,
        union_y_only=y_only,
        toy_factor=toy,
        aux_step=aux_train_step(n, d_aux),
        gcd_core=gcd_candidates(n),
    )


def _is_proper_pair(s: int, t: int, n: int) -> bool:
    return s > 1 and t > 1 and t < n


def format_report(
    r: MdcAnalyzeResult, d_aux: int | None = None, *, proper_only: bool = False
) -> str:
    lines: list[str] = []
    sn = math.isqrt(r.n)
    lines.append("=" * 64)
    lines.append(f"MDC ANALYZE — n = {r.n}")
    lines.append(f"Ecuación: (2x+3)(2y+3) = {r.n}   |   √n ≈ {sn}")
    lines.append("")

    lines.append("--- Tren X (variación x entero, paso Δx=1) ---")
    lines.append(f"  Pasos explorados : {r.train_x.steps}")
    if r.train_x.avg_slope is not None:
        lines.append(f"  Pendiente media Δy/Δx : {r.train_x.avg_slope}")
    if not r.train_x.collisions:
        lines.append("  Colisiones enteras : (ninguna)")
    for c in r.train_x.collisions:
        if proper_only and not _is_proper_pair(c.s, c.t, r.n):
            continue
        lines.append(
            f"  x={c.x}, y={c.y}  →  S={c.s}, T={c.t}, k={c.k}  "
            f"pendiente_local={c.slope or '—'}"
        )
    lines.append("")

    lines.append("--- Tren Y (variación y entero, paso Δy=1) ---")
    lines.append(f"  Pasos explorados : {r.train_y.steps}")
    if r.train_y.avg_slope is not None:
        lines.append(f"  Pendiente media Δx/Δy : {r.train_y.avg_slope}")
    if not r.train_y.collisions:
        lines.append("  Colisiones enteras : (ninguna)")
    for c in r.train_y.collisions:
        if proper_only and not _is_proper_pair(c.s, c.t, r.n):
            continue
        lines.append(
            f"  y={c.y}, x={c.x}  →  S={c.s}, T={c.t}, k={c.k}  "
            f"pendiente_local={c.slope or '—'}"
        )
    lines.append("")

    lines.append("--- Unión (resultados iguales X ∩ Y) ---")
    union = (
        [p for p in r.union_both if _is_proper_pair(p[0], p[1], r.n)]
        if proper_only
        else r.union_both
    )
    if union:
        for s, t in union:
            lines.append(f"  CONFIRMADO  {s} × {t} = {r.n}  (ambos trenes)")
    else:
        lines.append("  (sin intersección — revisar rango o n primo)")
    if r.union_x_only:
        lines.append(f"  Solo tren-X : {r.union_x_only}")
    if r.union_y_only:
        lines.append(f"  Solo tren-Y : {r.union_y_only}")
    lines.append("")

    if d_aux is None:
        d_aux = min(2 * sn + 3, max(5, r.n // 6))
    g = math.gcd(r.n, d_aux)
    lines.append(f"--- Recta auxiliar (velocidad constante Δx = {r.aux_step}) ---")
    lines.append(f"  y = -(n/d)x + n   con d={d_aux}, gcd(n,d)={g}")
    lines.append(f"  MCM(n,d) = {r.n // g * d_aux if g else '—'}")
    lines.append("")

    lines.append("--- Candidatos gcd(n, 2m+3) ---")
    for c in r.gcd_core[:12]:
        extra = f" => {c['S']}×{c['T']}" if c["T"] else ""
        lines.append(f"  m={c['m']}, S={c['S']}, gcd={c['gcd']}{extra}")
    if len(r.gcd_core) > 12:
        lines.append(f"  ... +{len(r.gcd_core) - 12} más")
    lines.append("")

    if r.toy_factor:
        lines.append(f"Factor toy (mdc_lib): {r.toy_factor}")
    lines.append("=" * 64)
    return "\n".join(lines)