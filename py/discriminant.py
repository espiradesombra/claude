"""
discriminant.py
===============
Discriminant Method: Deterministic Factorization Filter

Author: Víctor Manzanares Alberola (VMA)

Given N to test for primality, the method reduces factorization
to finding perfect squares of a deformed discriminant Δ(S) in a
narrow algebraic window around √N.

Key formulas:
    N = (2v+3)(2b+3)   =>   reparametrize:
    S  = 2v
    M  = N - 9
    Δ(S) = S² + 6S - M

    Factorization exists  ⟺  Δ(S) = k²  for integer k ≥ 0

Deterministic stop condition (geometric justification):
    4S + 16 > 2√Δ + 1   =>   STOP (gap > square spacing, no further roots)

This condition is NOT probabilistic — it follows from the fact that
consecutive perfect squares differ by 2k+1, and once the "stride"
of Δ(S) exceeds the local gap, no new perfect squares can appear.

Comparison with classical methods:
    - Fermat factorization:  seeks x² - N = y²  (additive)
    - This method:           seeks S²+6S-M = k²  (multiplicative window)
    - Advantage:             deterministic stop, no trial division needed
"""

import math
from typing import Optional, Tuple, List


# ─────────────────────────────────────────────
# Core discriminant functions
# ─────────────────────────────────────────────

def discriminant(S: int, M: int) -> int:
    """
    Compute Δ(S) = S² + 6S - M
    where M = N - 9 and S = 2v for a candidate factor (2v+3).
    """
    return S * S + 6 * S - M


def is_perfect_square(n: int) -> Tuple[bool, int]:
    """Return (True, sqrt) if n is a perfect square, else (False, isqrt(n))."""
    if n < 0:
        return False, 0
    k = math.isqrt(n)
    return (k * k == n), k


def stop_condition(S: int, delta: int) -> bool:
    """
    Deterministic stop: 4S + 16 > 2*sqrt(|Δ|) + 1
    When true, no further perfect-square solutions exist.
    """
    if delta < 0:
        return True  # Δ negative: no real solution possible
    sqrt_delta = math.isqrt(delta)
    return (4 * S + 16) > (2 * sqrt_delta + 1)


# ─────────────────────────────────────────────
# Main factorization method
# ─────────────────────────────────────────────

def discriminant_factor(N: int) -> Tuple[Optional[int], Optional[int], dict]:
    """
    Attempt to factor N using the discriminant method.

    Parameters
    ----------
    N : integer to factor (should be odd, not divisible by 3)

    Returns
    -------
    (f1, f2) if N = f1 * f2 with f1,f2 > 1
    (None, None) if N is likely prime
    info dict with step count, final S, and other diagnostics
    """
    info = {"steps": 0, "final_S": None, "reason": None}

    # Quick small factor checks
    for p in [2, 3, 5, 7, 11, 13]:
        if N % p == 0 and N > p:
            info["reason"] = f"small factor {p}"
            return p, N // p, info
        if N == p:
            info["reason"] = "small prime"
            return None, None, info

    # N must be ≡ 1 (mod 6) or ≡ 5 (mod 6) for (2v+3)(2b+3) factorization
    # (both factors ≡ ±1 mod 6 => product ≡ ±1 mod 6 ≡ 1 or 5 mod 6)
    if N % 6 not in (1, 5):
        info["reason"] = "not in 6k±1 form"
        return None, None, info

    M = N - 9
    sqrt_N = math.isqrt(N)

    # Start S from largest possible (near √N)
    v0 = (sqrt_N - 3) // 2
    if v0 < 0:
        v0 = 0
    S = 2 * v0

    while S >= 0:
        delta = discriminant(S, M)
        info["steps"] += 1

        if delta >= 0:
            is_sq, k = is_perfect_square(delta)
            if is_sq:
                # Recover factors
                # S = 2v, k = S + 3 - 2b - 3 = S - 2b => b = (S - k) / 2
                # But Δ = S² + 6S - M = (S+3)² - 9 - M = (S+3)² - N
                # So (S+3)² - k² = N => (S+3-k)(S+3+k) = N
                a = S + 3 - k
                b_val = S + 3 + k
                if a > 1 and b_val > 1 and a * b_val == N:
                    info["final_S"] = S
                    info["reason"] = "discriminant square"
                    return a, b_val, info
                # Try alternate recovery
                v = S // 2
                b = (k - S - 3) // 2
                f1 = 2 * v + 3
                f2 = 2 * b + 3
                if f1 > 1 and f2 > 1 and f1 * f2 == N:
                    info["final_S"] = S
                    info["reason"] = "discriminant square (alt)"
                    return f1, f2, info

        # Deterministic stop check
        if stop_condition(S, delta):
            info["final_S"] = S
            info["reason"] = "deterministic stop"
            break

        S -= 2

    return None, None, info


# ─────────────────────────────────────────────
# Primality filter (wrapper)
# ─────────────────────────────────────────────

def discriminant_primality(N: int) -> dict:
    """
    Use discriminant method as primality pre-filter.
    Returns classification and, if composite, the factors found.
    """
    f1, f2, info = discriminant_factor(N)

    result = {
        "N": N,
        "composite": f1 is not None,
        "factors": (f1, f2) if f1 else None,
        "steps": info["steps"],
        "reason": info["reason"],
    }

    return result


# ─────────────────────────────────────────────
# Δ trajectory analysis
# ─────────────────────────────────────────────

def delta_trajectory(N: int, show_all: bool = False) -> None:
    """
    Show the trajectory of Δ(S) as S decreases from √N.
    Highlights where perfect squares occur and where stop fires.
    """
    M = N - 9
    sqrt_N = math.isqrt(N)
    v0 = (sqrt_N - 3) // 2
    S = 2 * v0

    print(f"\nΔ trajectory for N = {N}  (M = {M})")
    print(f"{'S':>6}  {'Δ(S)':>12}  {'√Δ':>8}  {'perfect sq':>11}  {'stop':>5}")
    print("-" * 52)

    while S >= 0:
        delta = discriminant(S, M)
        is_sq = False
        k_str = ""
        if delta >= 0:
            is_sq, k = is_perfect_square(delta)
            k_str = str(k) if is_sq else ""
        stop = stop_condition(S, delta)

        marker = " ← FACTOR" if is_sq and k_str else (" ← STOP" if stop else "")

        if show_all or is_sq or stop:
            print(f"{S:>6}  {delta:>12}  {k_str:>8}  {'YES' if is_sq else 'no':>11}  "
                  f"{'YES' if stop else 'no':>5}{marker}")

        if stop:
            break
        S -= 2


# ─────────────────────────────────────────────
# Benchmark
# ─────────────────────────────────────────────

def benchmark(N_list: List[int], verbose: bool = True) -> None:
    """
    Run discriminant method on a list of integers and show results.
    """
    if verbose:
        print(f"\n{'N':>12}  {'result':>22}  {'steps':>6}  {'reason'}")
        print("-" * 65)

    for N in N_list:
        r = discriminant_primality(N)
        if r["composite"]:
            f1, f2 = r["factors"]
            result_str = f"{f1} × {f2}"
        else:
            result_str = "prime candidate"

        if verbose:
            print(f"{N:>12}  {result_str:>22}  {r['steps']:>6}  {r['reason']}")


def stress_test(limit: int = 1000) -> None:
    """
    Compare discriminant method vs trial division up to limit.
    Reports any discrepancies.
    """
    def is_prime_trial(n):
        if n < 2: return False
        if n == 2: return True
        if n % 2 == 0: return False
        for i in range(3, int(n**0.5)+1, 2):
            if n % i == 0: return False
        return True

    mismatches = 0
    for N in range(5, limit + 1, 2):  # odd numbers only
        is_prime_disc = not discriminant_primality(N)["composite"]
        is_prime_ref = is_prime_trial(N)
        if is_prime_disc != is_prime_ref:
            print(f"  MISMATCH at N={N}: discriminant={is_prime_disc}, trial={is_prime_ref}")
            mismatches += 1

    if mismatches == 0:
        print(f"✅ Stress test passed: all {limit//2} odd numbers up to {limit} correct.")
    else:
        print(f"⚠️  {mismatches} mismatches found.")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("Discriminant Method — Deterministic Factorization Filter")
    print("=" * 65)

    # 1. Benchmark on classic composites and primes
    test_cases = [
        15, 21, 35, 77, 91, 143, 221, 323, 391, 437,
        1001, 1763, 3599, 10403, 15251,
        # Known primes
        97, 101, 997, 1009, 9973, 10007
    ]
    benchmark(test_cases)

    print()

    # 2. Trajectory visualization for a composite
    delta_trajectory(143, show_all=True)
    delta_trajectory(10403, show_all=False)

    print()

    # 3. Stress test
    print("Stress test vs trial division (odd numbers up to 2000):")
    stress_test(2000)

    print()

    # 4. Step efficiency analysis
    print("\nStep efficiency (discriminant vs trial division):")
    print(f"{'N':>10}  {'disc_steps':>11}  {'trial_steps':>12}  {'speedup':>8}")
    print("-" * 48)

    def trial_steps(N):
        steps = 0
        for i in range(3, int(N**0.5)+1, 2):
            steps += 1
            if N % i == 0:
                return steps
        return steps

    for N in [143, 1001, 10403, 100621, 999983]:
        r = discriminant_primality(N)
        ds = r["steps"]
        ts = trial_steps(N)
        speedup = ts / ds if ds > 0 else float('inf')
        print(f"{N:>10}  {ds:>11}  {ts:>12}  {speedup:>7.1f}x")
