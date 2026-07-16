"""Puente opcional a vma-methods cuando no hay DLL o para comparar Python."""

from __future__ import annotations

import math
import sys
from pathlib import Path

_VMA = Path(__file__).resolve().parents[3] / "vma-methods"
if _VMA.is_dir() and str(_VMA) not in sys.path:
    sys.path.insert(0, str(_VMA))


def sieve_desmemoriada_count(limit: int) -> int:
    try:
        from vma_methods.cribas import CribaDesmemoriada

        return len(CribaDesmemoriada(limit).run())
    except ImportError:
        return len(sieve_primes_python(limit))


def sieve_count_python(limit: int) -> int:
    try:
        from vma_methods.cribas import CribaModular6k

        return len(CribaModular6k(limit).run())
    except ImportError:
        return len(sieve_primes_python(limit))


def sieve_primes_python(limit: int) -> list[int]:
    if limit < 2:
        return []
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    root = int(limit**0.5)
    for p in range(2, root + 1):
        if not is_prime[p]:
            continue
        is_prime[p * p : limit + 1 : p] = b"\x00" * len(is_prime[p * p : limit + 1 : p])
    return [i for i in range(2, limit + 1) if is_prime[i]]


def newton_log_python(
    E: float,
    b: float,
    familia: str,
    *,
    n_exp: int = 2,
    k_known: float = 1.0,
) -> dict:
    try:
        from vma_methods.newton import log_con_oraculo

        kw = {}
        if familia == "potencia":
            kw["n"] = n_exp
        if familia == "kp":
            kw["k"] = k_known
        return log_con_oraculo(E, b=b, familia=familia, **kw)
    except ImportError:
        import math

        j = math.log(E) / math.log(b)
        return {
            "j": j,
            "j_exacto": j,
            "iteraciones": 0,
            "error": 0.0,
            "converged": True,
        }


def criva_python(x: float, layers: int, iterations: int) -> float:
    try:
        from vma_methods.criva import criva

        return float(criva(x, layers=layers, iterations=iterations))
    except ImportError:
        d = 0.66025
        for _ in range(iterations):
            d = (d + d * (1.0 - 1.0 / math.log(x + 1.0))) * 0.5
        return d