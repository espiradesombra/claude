"""Criba 6k±1 — híbrida y desmemoriada."""
from __future__ import annotations

import math


def es_candidato_6k(v: int) -> bool:
    return v % 6 in (1, 5)


def criba_desmemoriada(limite: int) -> list[int]:
    """Criba 6k±1 clásica (wheel mod 6)."""
    if limite < 2:
        return []
    p6 = [True] * (limite + 1)
    p6[0] = p6[1] = False
    for i in range(2, limite + 1):
        if i % 2 == 0 or i % 3 == 0:
            p6[i] = False
    raiz = int(limite**0.5)
    for p in range(5, raiz + 1, 6):
        if not p6[p]:
            continue
        for multiple in range(p * p, limite + 1, 6 * p):
            p6[multiple] = False
            if multiple + 2 * p <= limite:
                p6[multiple + 2 * p] = False
    return [i for i, ok in enumerate(p6) if ok]


def criba_hibrida(n: int) -> list[int]:
    """Criba híbrida con saltos 2p/4p (archivos-vma)."""
    p = [0] * (n + 1)
    for v in range(2, n + 1):
        if v in (2, 3) or es_candidato_6k(v):
            p[v] = 1
    m = n // 2
    for primo in range(5, m + 1):
        if p[primo] != 1:
            continue
        j = primo * primo
        salto1 = 2 * primo
        salto2 = 4 * primo
        toggle = True
        while j <= m:
            p[j] = 0
            j += salto1 if toggle else salto2
            toggle = not toggle
    for primo in range(5, int(math.sqrt(n)) + 1):
        if p[primo] != 1:
            continue
        j = (n // primo) * primo
        while j > m:
            if es_candidato_6k(j):
                p[j] = 0
            j -= primo
    return [i for i, flag in enumerate(p) if flag == 1]


def es_primo(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True