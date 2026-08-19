"""Regla mecánica VMA — pendiente, choque, (2v+3)(2l+3). Sin DLL."""

from __future__ import annotations

import math


def f_desde_v(v: int) -> int:
    return 2 * v + 3


def f_desde_l(l: int) -> int:
    return 2 * l + 3


def encoger_en_choque(x: int, pendiente: int) -> int:
    if pendiente == 0:
        return x
    return x - (x + 1) // pendiente


def choque_mitad(x: int) -> int:
    return (x + 1) // 2


def pendiente_desde_v(n: int, v: int) -> int:
    f1 = f_desde_v(v)
    return n // f1 if f1 else 0


def delta_pendiente_v(n: int, v: int) -> int:
    k0 = pendiente_desde_v(n, v)
    k1 = pendiente_desde_v(n, v + 1)
    return max(0, k0 - k1)


def factorizar_par(n: int) -> tuple[bool, int, int, int, int]:
    if n < 9 or n % 2 == 0:
        return False, 0, 0, 0, 0
    sq = math.isqrt(n)
    v_max = (sq - 3) // 2
    v_ini = 0
    if v_max > 4:
        v_cota = choque_mitad(v_max)
        if 0 < v_cota < v_max:
            v_ini = max(v_ini, v_cota // 3)

    for v in range(v_ini, v_max + 1):
        f1 = f_desde_v(v)
        if n % f1:
            continue
        f2 = n // f1
        if f2 < 3 or f2 % 2 == 0 or (f2 - 3) % 2:
            continue
        l = (f2 - 3) // 2
        return True, v, l, min(f1, f2), max(f1, f2)
    return False, 0, 0, 0, 0


if __name__ == "__main__":
    print("Regla mecanica — demo")
    x = 100
    print(f"x={x} choque_mitad={choque_mitad(x)}")
    for p in range(2, 6):
        x = encoger_en_choque(x, p)
        print(f"  pendiente {p} -> x={x}")

    for n, esp in [(15, 3), (101 * 103, 101), (100003 * 100019, 100003)]:
        ok, v, l, f1, f2 = factorizar_par(n)
        mark = "OK" if ok and f1 == esp else "FAIL"
        print(f"N={n} -> ({f1})*({f2}) v={v} l={l} k_v={pendiente_desde_v(n,v)} [{mark}]")
        print(f"  delta v->v+1: {delta_pendiente_v(n, v)}")