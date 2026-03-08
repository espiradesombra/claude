"""
salto_maximo.py
===============
Maximum Prime Gap Conjecture / Theorem

Author: Víctor Manzanares Alberola (VMA)

Statement:
    In the interval  [n - floor(sqrt(n+3)) - 3,  n + 3]
    there are at least TWO primes.

Equivalently (via the density estimator):
    The minimum safe count of primes in a window of length sqrt(n+3) is:
        (1 - (e-2)) * sqrt(n+3)  >=  2

    This holds for n > ~100 (verified), making the conjecture a theorem
    in that range. Below n ~ 100 it must be checked case by case.

Key constants:
    e - 2  ≈  0.718...   (overcount fraction from factorial accumulator)
    1 - (e-2)  ≈  0.282  (minimum prime density factor)

    Condition: (1-(e-2)) * sqrt(n+3) >= 2
              => sqrt(n+3) >= 2 / (1-(e-2))  ≈  7.09
              => n+3 >= 50.3
              => n >= ~47  (conservative: n > ~100 for safety margin)

Connection to MRAUV:
    The accumulator m(n) = sum_{i=2}^{K} sqrt(n+3)/i! converges to
    (e-2)*sqrt(n+3), which is the overcount of composites.
    The residual L(n) - m(n) >= 2 is the MRAUV criterion.
"""

import math
from typing import List, Tuple


# ─────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────

E_MINUS_2 = math.e - 2          # ≈ 0.71828...
DENSITY_FACTOR = 1 - E_MINUS_2  # ≈ 0.28172...


# ─────────────────────────────────────────────
# Window definition
# ─────────────────────────────────────────────

def ventana(n: int) -> Tuple[int, int]:
    """
    Return the interval [lo, hi] = [n - floor(sqrt(n+3)) - 3,  n+3]
    This is the search window guaranteed to contain >= 2 primes.
    """
    radio = math.isqrt(n + 3)
    lo = n - radio - 3
    hi = n + 3
    return lo, hi


def longitud_ventana(n: int) -> int:
    """Length of the search window = floor(sqrt(n+3)) + 6."""
    return math.isqrt(n + 3) + 6


# ─────────────────────────────────────────────
# Density bound
# ─────────────────────────────────────────────

def cota_minima_primos(n: int) -> float:
    """
    Theoretical minimum count of primes in window of length sqrt(n+3):
        (1 - (e-2)) * sqrt(n+3)

    When this >= 2, the conjecture is guaranteed (for n > ~100).
    """
    return DENSITY_FACTOR * math.sqrt(n + 3)


def n_minimo_teorema() -> float:
    """
    Find the n above which cota_minima_primos(n) >= 2:
        (1-(e-2)) * sqrt(n+3) >= 2
        n+3 >= (2 / DENSITY_FACTOR)^2
        n   >= (2 / DENSITY_FACTOR)^2 - 3
    """
    return (2 / DENSITY_FACTOR) ** 2 - 3


# ─────────────────────────────────────────────
# Primality and counting
# ─────────────────────────────────────────────

def sieve(limit: int) -> List[bool]:
    """Sieve of Eratosthenes, returns is_prime list."""
    if limit < 2:
        return [False] * (limit + 1)
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return is_prime


def contar_primos_en_ventana(n: int, is_prime: List[bool]) -> Tuple[int, List[int]]:
    """
    Count and list primes in the window [lo, hi] for given n.
    Returns (count, list_of_primes).
    """
    lo, hi = ventana(n)
    lo = max(lo, 2)
    hi = min(hi, len(is_prime) - 1)
    primos = [p for p in range(lo, hi + 1) if is_prime[p]]
    return len(primos), primos


# ─────────────────────────────────────────────
# Verification
# ─────────────────────────────────────────────

def verificar_conjetura(n_max: int = 10_000, verbose: bool = True) -> bool:
    """
    Verify that every n in [50, n_max] has at least 2 primes in its window.

    Returns True if conjecture holds throughout, False on first failure.
    """
    is_prime = sieve(n_max + 10)
    fallos = []

    for n in range(50, n_max + 1):
        count, primos = contar_primos_en_ventana(n, is_prime)
        if count < 2:
            fallos.append((n, count, primos))

    if verbose:
        if fallos:
            print(f"⚠️  Conjecture FAILS at {len(fallos)} values:")
            for n, c, p in fallos[:10]:
                print(f"   n={n}: only {c} prime(s) in window: {p}")
        else:
            print(f"✅ Conjecture verified for all n in [50, {n_max}]: "
                  f"always >= 2 primes in window.")

    return len(fallos) == 0


def tabla_verificacion(n_values: List[int], is_prime: List[bool]) -> None:
    """
    Print verification table for selected n values.
    """
    n_min_teo = n_minimo_teorema()

    print(f"\nSaltoMáximo verification table")
    print(f"(1-(e-2)) = {DENSITY_FACTOR:.5f},  n_min_teorema ≈ {n_min_teo:.1f}")
    print()
    print(f"{'n':>8}  {'[lo, hi]':>18}  {'cota':>7}  "
          f"{'#primos':>8}  {'primos en ventana':>30}  {'ok':>4}")
    print("-" * 85)

    for n in n_values:
        lo, hi = ventana(n)
        cota = cota_minima_primos(n)
        count, primos = contar_primos_en_ventana(n, is_prime)
        ok = "✅" if count >= 2 else "❌"
        primos_str = str(primos[:6]) + ("..." if len(primos) > 6 else "")
        print(f"{n:>8}  [{lo:>5}, {hi:>5}]  {cota:>7.3f}  "
              f"{count:>8}  {primos_str:>30}  {ok:>4}")


# ─────────────────────────────────────────────
# MRAUV connection: L(n) - m(n) >= 2
# ─────────────────────────────────────────────

def L(n: int) -> float:
    """Search corridor length: floor(sqrt(n+3)) + 7"""
    return math.isqrt(n + 3) + 7


def m(n: int, K: int = None) -> float:
    """
    Factorial accumulator (overcount of composites).
    m(n) = sum_{i=2}^{K} sqrt(n+3)/i!
    Converges to (e-2)*sqrt(n+3).
    """
    if K is None:
        K = max(2, int(math.isqrt(n) * 9 // 24))
    sqrt_n3 = math.sqrt(n + 3)
    total = 0.0
    factorial = 1
    for i in range(1, K + 1):
        factorial *= i
        if i >= 2:
            total += sqrt_n3 / factorial
    return total


def mrauv_criterion(n: int) -> float:
    """
    MRAUV criterion value: L(n) - m(n)
    Must be >= 2 for the conjecture to hold at n.
    """
    return L(n) - m(n)


def convergencia_acumulador(n: int = 1000, K_max: int = 20) -> None:
    """
    Show how m(n)/sqrt(n+3) converges to (e-2) as K increases.
    """
    sqrt_n3 = math.sqrt(n + 3)
    print(f"\nConvergence of m(n)/sqrt(n+3) → (e-2) for n={n}")
    print(f"Target: e-2 = {E_MINUS_2:.8f}")
    print(f"{'K':>4}  {'m(n)':>12}  {'m/sqrt':>10}  {'error vs e-2':>14}")
    print("-" * 46)

    factorial = 1
    total = 0.0
    for i in range(1, K_max + 1):
        factorial *= i
        if i >= 2:
            total += sqrt_n3 / factorial
            ratio = total / sqrt_n3
            err = abs(ratio - E_MINUS_2)
            print(f"{i:>4}  {total:>12.6f}  {ratio:>10.6f}  {err:>14.2e}")


# ─────────────────────────────────────────────
# Gap statistics
# ─────────────────────────────────────────────

def estadisticas_salto(n_max: int = 100_000) -> None:
    """
    Compare actual prime gaps vs sqrt(n) bound.
    Checks SaltoMáximo(n) <= sqrt(n+3) empirically.
    """
    is_prime = sieve(n_max)
    primes = [i for i, v in enumerate(is_prime) if v]

    max_gap_real = 0
    max_gap_n = 0
    violaciones = []

    for i in range(1, len(primes)):
        gap = primes[i] - primes[i-1]
        n = primes[i-1]
        cota = math.sqrt(n + 3)
        if gap > max_gap_real:
            max_gap_real = gap
            max_gap_n = n
        if gap > cota:
            violaciones.append((n, primes[i], gap, cota))

    print(f"\nPrime gap statistics up to {n_max}")
    print(f"Max observed gap: {max_gap_real} at n={max_gap_n}")
    print(f"sqrt(n+3) at max gap: {math.sqrt(max_gap_n+3):.2f}")
    print(f"Ratio max_gap/sqrt(n+3): {max_gap_real/math.sqrt(max_gap_n+3):.4f}")

    if violaciones:
        print(f"\n⚠️  {len(violaciones)} gaps > sqrt(n+3):")
        for n, p2, gap, cota in violaciones[:5]:
            print(f"   gap({n},{p2})={gap} > sqrt({n}+3)={cota:.2f}")
    else:
        print(f"✅ No gaps > sqrt(n+3) found up to {n_max}.")
        print(f"   SaltoMáximo(n) ≤ sqrt(n+3) holds throughout.")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("SaltoMáximo — Maximum Prime Gap Conjecture")
    print("Author: Víctor Manzanares Alberola")
    print("=" * 65)

    # 1. Theoretical threshold
    n_teo = n_minimo_teorema()
    print(f"\nTheoretical threshold n_min = {n_teo:.1f}")
    print(f"  (1-(e-2)) * sqrt(n+3) >= 2  for n >= {n_teo:.0f})")
    print(f"  Conservative validity: n > 100")

    # 2. Convergence of accumulator to (e-2)
    convergencia_acumulador(n=500, K_max=15)

    # 3. Verification
    print()
    verificar_conjetura(n_max=10_000)

    # 4. Verification table
    print()
    is_prime_ref = sieve(20_000)
    sample = [50, 100, 200, 500, 1000, 2000, 5000, 10000]
    tabla_verificacion(sample, is_prime_ref)

    # 5. MRAUV criterion values
    print(f"\nMRAUV criterion L(n) - m(n) >= 2:")
    print(f"{'n':>8}  {'L(n)':>8}  {'m(n)':>10}  {'L-m':>8}  {'>=2':>5}")
    print("-" * 45)
    for n in sample:
        ln = L(n)
        mn = m(n)
        diff = ln - mn
        ok = "✅" if diff >= 2 else "❌"
        print(f"{n:>8}  {ln:>8.3f}  {mn:>10.5f}  {diff:>8.4f}  {ok:>5}")

    # 6. Gap statistics
    print()
    estadisticas_salto(n_max=100_000)
