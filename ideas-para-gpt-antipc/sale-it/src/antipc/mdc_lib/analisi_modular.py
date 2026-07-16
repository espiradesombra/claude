"""Análisis modular MDC — port parcial analisi_modular.hpp."""

from __future__ import annotations

import math

PI2 = 2.0 * math.pi


def fase_modular(ny: int, m: int) -> float:
    if m == 0:
        return 0.0
    return PI2 * (ny % m) / m


def patron_bits(ny: int, ventana: int = 8) -> int:
    if ventana == 0 or ventana > 64:
        return 0
    transiciones = 0
    mascara = (1 << ventana) - 1 if ventana < 64 else (1 << 64) - 1
    bits = ny & mascara
    ult = bits & 1
    bits >>= 1
    for _ in range(1, ventana):
        act = bits & 1
        if act != ult:
            transiciones += 1
            ult = act
        bits >>= 1
    return transiciones


def curvatura_modular(ny: int, m: int) -> float:
    if ny < 2 or m == 0:
        return 0.0
    f0 = fase_modular(ny - 1, m)
    f1 = fase_modular(ny, m)
    f2 = fase_modular(ny + 1, m)
    v1, v2 = f1 - f0, f2 - f1
    acc = v2 - v1
    denom = math.pow(1.0 + v1 * v1, 1.5)
    return 0.0 if abs(denom) < 1e-12 else abs(acc) / denom