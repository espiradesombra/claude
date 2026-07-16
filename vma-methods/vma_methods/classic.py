"""Utilidades clásicas — referencia para validar métodos VMA."""

from __future__ import annotations

import math

from .criva import sieve_primes


def pi_count(x: int) -> int:
    """π(x) — número de primos ≤ x (Eratóstenes)."""
    if x < 2:
        return 0
    return len(sieve_primes(x))


def prime_density(x: float) -> float:
    """π(x)/x."""
    if x < 2:
        return 0.0
    return pi_count(int(x)) / x


def pnt_estimate(x: float) -> float:
    """Aproximación PNT: x / ln(x)."""
    if x <= 1:
        return 0.0
    return x / math.log(x)


def is_prime_trial(n: int) -> bool:
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True