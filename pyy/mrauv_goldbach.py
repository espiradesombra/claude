"""
mrauv_goldbach.py
=================
MRAUV (Modelo de Recorrido Acumulado por Velocidad)
Applied to Goldbach's conjecture via asymmetric fault measurement.

Author: Víctor Manzanares Alberola (VMA)
Theory: The number of Goldbach decompositions G(2n) is bounded below
        by D(n) - F_eff(n)/(2n), where D(n) is the MRAUV density
        prediction and F_eff(n) is the effective asymmetric fault
        caused by multiples of small primes interfering with the
        symmetry algorithm.

Conjecture: For all n > N0, D(n) > F_eff(n)/(2n) + epsilon > 0
            => Goldbach holds for 2n > 2*N0
"""

import math
from typing import List, Tuple


# ─────────────────────────────────────────────
# Core MRAUV functions
# ─────────────────────────────────────────────

def L(n: int) -> float:
    """Search corridor length around n."""
    return math.isqrt(n + 3) + 7


def m(n: int, K: int = None) -> float:
    """
    Overcount of composites in the search corridor.
    K defaults to floor(floor(sqrt(n)) * 9/24).
    Uses factorial decay: sum_{i=2}^{K} sqrt(n+3) / i!
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


def D(n: int) -> float:
    """
    MRAUV predicted prime density in [n, 2n].
    D(n) = (L(n) - m(n)) / (2n)
    """
    return (L(n) - m(n)) / (2 * n)


def sieve_primes(limit: int) -> List[int]:
    """Sieve of Eratosthenes up to limit."""
    if limit < 2:
        return []
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i, v in enumerate(is_prime) if v]


def F_eff(n: int, small_primes: List[int]) -> float:
    """
    Effective asymmetric fault at n.
    Counts multiples of small primes in [1, 2n] weighted by
    local prime density — these are the positions the symmetry
    algorithm may misclassify.

    F_eff(n) ≈ Σ_{p ≤ sqrt(2n)} floor(2n/p) · π(2n)/(2n)
    """
    bound = int(math.sqrt(2 * n)) + 1
    pi_2n = sum(1 for p in small_primes if p <= 2 * n)
    density_2n = pi_2n / (2 * n) if n > 0 else 0

    fault = 0.0
    for p in small_primes:
        if p > bound:
            break
        fault += (2 * n // p) * density_2n
    return fault


def goldbach_margin(n: int, small_primes: List[int]) -> float:
    """
    Margin = D(n) - F_eff(n)/(2n)
    Positive margin => Goldbach criterion satisfied at n.
    """
    d = D(n)
    f = F_eff(n, small_primes) / (2 * n)
    return d - f


# ─────────────────────────────────────────────
# MRAUV piecewise kinematic model
# ─────────────────────────────────────────────

def mrauv_fit(n0: int, n1: int, n2: int) -> Tuple[float, float, float]:
    """
    Fit MRAUV parameters (D0, V0, a0) from three sample points.
    Uses finite differences:
        V0 ≈ (D(n1) - D(n0)) / (n1 - n0)
        a0 ≈ (D(n2) - 2*D(n1) + D(n0)) / ((n1-n0)^2)
    Returns (D0, V0, a0)
    """
    D0 = D(n0)
    D1 = D(n1)
    D2 = D(n2)
    h = n1 - n0
    V0 = (D1 - D0) / h
    a0 = (D2 - 2 * D1 + D0) / (h * h)
    return D0, V0, a0


def mrauv_predict(n0: int, D0: float, V0: float, a0: float, delta_n: int) -> float:
    """
    Predict density at n0 + delta_n using kinematic model:
    D_pred = D0 + V0*delta_n + 0.5*a0*delta_n^2
    """
    return D0 + V0 * delta_n + 0.5 * a0 * delta_n * delta_n


# ─────────────────────────────────────────────
# Verification
# ─────────────────────────────────────────────

def verify_goldbach_MRAUV(
    n_max: int = 100_000,
    delta: int = 5_000,
    n_start: int = 1_000,
    verbose: bool = True
) -> bool:
    """
    Verify the MRAUV-Goldbach criterion over [n_start, n_max].

    For each segment, checks:
        D(n) > F_eff(n)/(2n)   =>   margin > 0

    Returns True if criterion holds throughout, False on first failure.
    """
    small_primes = sieve_primes(int(math.sqrt(2 * n_max)) + 2)

    if verbose:
        print(f"{'n':>8}  {'D(n)':>10}  {'F/2n':>10}  {'margin':>10}  {'status':>8}")
        print("-" * 55)

    n = n_start
    all_ok = True
    while n <= n_max:
        d = D(n)
        f_norm = F_eff(n, small_primes) / (2 * n)
        margin = d - f_norm
        status = "✅" if margin > 0 else "⚠️  ALERT"

        if verbose:
            print(f"{n:>8}  {d:>10.6f}  {f_norm:>10.6f}  {margin:>+10.6f}  {status}")

        if margin <= 0:
            all_ok = False

        n += delta

    if verbose:
        if all_ok:
            print("\n✅ Goldbach MRAUV criterion satisfied in entire range.")
        else:
            print("\n⚠️  Criterion failed at one or more points.")

    return all_ok


def count_goldbach_decompositions(two_n: int, primes_set: set) -> int:
    """
    Count exact Goldbach decompositions of 2n (brute force, for validation).
    Returns number of pairs (p, q) with p <= q, p+q = 2n, both prime.
    """
    n = two_n // 2
    count = 0
    primes_list = sorted(p for p in primes_set if p <= n)
    for p in primes_list:
        q = two_n - p
        if q in primes_set:
            count += 1
    return count


def compare_mrauv_vs_exact(
    two_n_values: List[int],
    verbose: bool = True
) -> None:
    """
    Compare MRAUV prediction with exact Goldbach decomposition count.
    """
    max_val = max(two_n_values)
    all_primes = sieve_primes(max_val)
    primes_set = set(all_primes)
    small_primes = sieve_primes(int(math.sqrt(max_val)) + 2)

    if verbose:
        print(f"{'2n':>10}  {'G_exact':>9}  {'D(n)·2n':>10}  {'F_eff':>9}  {'margin':>9}")
        print("-" * 55)

    for two_n in two_n_values:
        n = two_n // 2
        g_exact = count_goldbach_decompositions(two_n, primes_set)
        d_pred = D(n) * two_n
        f = F_eff(n, small_primes)
        margin = g_exact - f

        if verbose:
            print(f"{two_n:>10}  {g_exact:>9}  {d_pred:>10.1f}  {f:>9.1f}  {margin:>+9.1f}")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("MRAUV-Goldbach: Density criterion verification")
    print("=" * 60)
    print()

    # 1. Verify criterion over range
    verify_goldbach_MRAUV(n_max=50_000, delta=5_000)

    print()
    print("=" * 60)
    print("MRAUV vs exact Goldbach decompositions")
    print("=" * 60)
    print()

    # 2. Compare with exact counts
    sample = [100, 500, 1000, 5000, 10000, 50000]
    compare_mrauv_vs_exact([2 * n for n in sample])

    print()
    print("=" * 60)
    print("MRAUV kinematic fit (3-point)")
    print("=" * 60)
    print()

    # 3. Fit kinematic parameters at n=10000
    n0, n1, n2 = 10_000, 20_000, 30_000
    D0, V0, a0 = mrauv_fit(n0, n1, n2)
    print(f"Fit at n0={n0}: D0={D0:.6f}, V0={V0:.2e}, a0={a0:.2e}")
    for dn in [10_000, 20_000, 30_000, 50_000]:
        pred = mrauv_predict(n0, D0, V0, a0, dn)
        real = D(n0 + dn)
        err = abs(pred - real) / real * 100 if real else 0
        print(f"  n={n0+dn:6d}: predicted={pred:.6f}, real={real:.6f}, err={err:.3f}%")
