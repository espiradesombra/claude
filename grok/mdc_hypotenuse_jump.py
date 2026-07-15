#!/usr/bin/env python3
"""
MDC / Hipotenusa — saltos MCM·gcd y cruces enteros
Víctor Manzanares Alberola · uso local educativo

Recta diofántica exacta:
    (2x+3)(2y+3) = n  <=>  y = (n - (6x - 9) / (4x + 6)

Recta auxiliar (pendiente -n/d):
    y = -(n/d)*x + n

Paso mínimo para y entero en la auxiliar:
    Δx = d / gcd(n, d) = lcm(n, d) / n
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction


def gcd(a: int, b: int) -> int:
    return math.gcd(abs(a), abs(b))


def lcm(a: int, b: int) -> int:
    g = gcd(a, b)
    if g == 0:
        return 0
    return abs(a // g * b)


def sqrt_int(n: int) -> int | None:
    r = math.isqrt(n)
    return r if r * r == n else None


def diofantic_y(n: int, x: Fraction) -> Fraction | None:
    """y tal que (2x+3)(2y+3)=n; None si denominador 0."""
    den = 4 * x + 6
    if den == 0:
        return None
    return Fraction(n - 6 * x - 9, den)


def factor_from_xy(x: Fraction, y: Fraction, n: int) -> tuple[int, int] | None:
    S = 2 * x + 3
    T = 2 * y + 3
    if S.denominator != 1 or T.denominator != 1:
        return None
    s, t = int(S), int(T)
    if s > 0 and t > 0 and s * t == n:
        return s, t
    return None


def scan_diofantic(n: int, x_max: int | None = None) -> list[dict]:
    """Cruces enteros (x,y) y semi-enteros útiles en la recta exacta."""
    if x_max is None:
        x_max = n // 2
    hits: list[dict] = []
    seen: set[tuple[int, int]] = set()

    for xi in range(-3, x_max + 1):
        x = Fraction(xi, 1)
        y = diofantic_y(n, x)
        if y is None:
            continue
        pair = factor_from_xy(x, y, n)
        if pair and pair not in seen:
            seen.add(pair)
            s, t = pair
            k = Fraction(t - s, 2)
            hits.append(
                {
                    "x": xi,
                    "y": y,
                    "S": s,
                    "T": t,
                    "k": k,
                    "kind": "entero" if y.denominator == 1 else "racional",
                }
            )
        elif y.denominator == 1 and int(y) >= 0:
            s = 2 * xi + 3
            t = 2 * int(y) + 3
            if s > 0 and t > 0 and s * t == n and (s, t) not in seen:
                seen.add((s, t))
                hits.append(
                    {
                        "x": xi,
                        "y": y,
                        "S": s,
                        "T": t,
                        "k": Fraction(t - s, 2),
                        "kind": "entero",
                    }
                )

    # Semi-enteros: x = m + 1/2  <=>  S = 2m+4, útil cuando n cuadrado
    for m in range(-2, x_max):
        x = Fraction(2 * m + 1, 2)
        y = diofantic_y(n, x)
        if y is None:
            continue
        pair = factor_from_xy(x, y, n)
        if pair and pair not in seen:
            seen.add(pair)
            s, t = pair
            hits.append(
                {
                    "x": str(x),
                    "y": y,
                    "S": s,
                    "T": t,
                    "k": Fraction(t - s, 2),
                    "kind": "semi-entero",
                }
            )

    hits.sort(key=lambda h: h["S"])
    return hits


def aux_line_hits(n: int, d: int, x_cap: int | None = None) -> list[tuple[int, int]]:
    """Puntos (x,y) enteros en y = -(n/d)x + n avanzando Δx = d/gcd(n,d)."""
    if d == 0:
        return []
    g = gcd(n, d)
    step = d // g
    if x_cap is None:
        x_cap = max(d, n // 2)
    pts: list[tuple[int, int]] = []
    x = 0
    while x <= x_cap:
        # y = n - (n/d)*x  con aritmética exacta
        y = n - (n * x) // d
        if (n * x) % d == 0:
            pts.append((x, y))
        x += step
    return pts


def mdc_candidates(n: int, m_max: int | None = None) -> list[dict]:
    """Candidatos S = 2m+3 con gcd(n,S) > 1 (núcleo MDC)."""
    if m_max is None:
        m_max = math.isqrt(n)
    out: list[dict] = []
    for m in range(0, m_max + 1):
        S = 2 * m + 3
        g = gcd(n, S)
        if g > 1:
            out.append({"m": m, "S": S, "gcd": g, "T": n // S if n % S == 0 else None})
    return out


def triangle_info(n: int) -> dict:
    sn = sqrt_int(n)
    m_vis = Fraction(2 - (sn if sn else math.isqrt(n)), 2 - Fraction(n, 2))
    return {
        "A": (Fraction(n, 2), 2),
        "B": (2, sn if sn else f"√{n} (no entero)"),
        "pendiente_visual": m_vis,
        "ecuacion_exacta": "y = (n - 6x - 9) / (4x + 6)",
    }


def fmt_frac(z: Fraction) -> str:
    return str(z) if z.denominator != 1 else str(z.numerator)


def report(n: int, d: int | None = None) -> None:
    print("=" * 60)
    print(f"n = {n}")
    tri = triangle_info(n)
    print(f"Triángulo visual: A={tri['A']}, B={tri['B']}")
    print(f"Pendiente A→B (aprox.): {tri['pendiente_visual']}")
    print(f"Recta exacta factores: {tri['ecuacion_exacta']}")
    print()

    print("--- Cruces diofánticos (factores) ---")
    hits = scan_diofantic(n)
    if not hits:
        print("  (ninguno en rango explorado)")
    for h in hits:
        k = h["k"]
        print(
            f"  x={h['x']}, y={fmt_frac(h['y'])}  "
            f"S={h['S']}, T={h['T']},  T-S=2k con k={fmt_frac(k)}  [{h['kind']}]"
        )
        if k.denominator == 1:
            ki = int(k)
            print(f"    S+k={h['S']+ki}, S-k={h['S']-ki}")
    print()

    if d is None:
        # d de ejemplo: denominador de pendiente -n/d en forma reducida para recta (0,n)-(d,0)
        d = min(2 * math.isqrt(n) + 3, max(5, n // 6))

    g = gcd(n, d)
    step = d // g
    print(f"--- Recta auxiliar y = -(n/d)x + n  con d={d} ---")
    print(f"  gcd(n,d)={g},  MCM(n,d)={lcm(n,d)}")
    print(f"  Paso Δx = d/gcd = {step}  (= MCM/d ... no, = MCM/n cuando aplica)")
    print(f"  Fracción {{n/d}} = {Fraction(n, d) - (n // d)}")
    aux = aux_line_hits(n, d)
    for x, y in aux[:12]:
        mark = ""
        if y >= 0:
            pair = factor_from_xy(Fraction(x), Fraction(y), n)
            if pair:
                mark = f"  --> factor {pair[0]}×{pair[1]}"
        print(f"  (x={x}, y={y}){mark}")
    if len(aux) > 12:
        print(f"  ... +{len(aux)-12} puntos más")
    print()

    print("--- Candidatos MDC (gcd(n, 2m+3) > 1) ---")
    for c in mdc_candidates(n)[:15]:
        t = c["T"]
        extra = f" => T={t}" if t else ""
        print(f"  m={c['m']}, S={c['S']}, gcd={c['gcd']}{extra}")
    print("=" * 60)


def main() -> None:
    p = argparse.ArgumentParser(description="Hipotenusa + saltos MCM/gcd + factores S,T,k")
    p.add_argument("n", type=int, nargs="?", default=36, help="entero a factorizar (default 36)")
    p.add_argument("-d", type=int, default=None, help="denominador d de pendiente -n/d")
    p.add_argument("--demo", action="store_true", help="ejemplos 36, 143, 35")
    args = p.parse_args()

    if args.demo:
        for val in (36, 35, 143):
            report(val, args.d)
            print()
    else:
        report(args.n, args.d)


if __name__ == "__main__":
    main()