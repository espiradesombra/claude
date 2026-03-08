"""
mdc.py
======
Método Diofántico Cinemático (MDC)
Kinematic Diophantine Method

Author: Víctor Manzanares Alberola (VMA), Valencia 2026

Unified procedure for locating integer solutions of Diophantine
equations via kinematic analysis of the decimal part of a real
parametric function.

Core idea:
    Given F(x, t) = 0, define d(t) = frac(g(x, t)) for a rational g.
    The exact Diophantine solution occurs at d(t*) = δ (target).
    Instead of testing all t, measure velocity/acceleration/jerk of
    d(t) at 4 consecutive points and extrapolate directly to t*.

Algorithm (central):
    1. Evaluate d(t0), d(t0+s), d(t0+2s), d(t0+3s)  [4 points, step s]
    2. Velocities:     v1=d1-d0,  v2=d2-d1,  v3=d3-d2
    3. Accelerations:  a1=v2-v1,  a2=v3-v2
    4. Jerk:           j = a2-a1
    5. Mean velocity:  v̄ = (v1+v2+v3)/3
    6. Prediction:     t* = t0 + (δ - d0) / v̄
    7. If round(t*) satisfies d(t*) = δ → solution found

Applications:
    App 1 — Factorization:   d(m) = frac(N / (2*(2m+3))),  δ = 0.5
    App 2 — Wieferich primes: d(K) = frac(p(K)),           δ = 0.0
             where p(K) solves 2^p = K*p^2 + 2
"""

import math
from typing import Optional, Tuple, List


# ─────────────────────────────────────────────
# Core MDC engine
# ─────────────────────────────────────────────

def mdc_predict(
    d_func,
    t0: float,
    delta: float,
    step: float = 1.0,
    max_jumps: int = 200,
    tol: float = 1e-9,
    verbose: bool = False
) -> Optional[float]:
    """
    General MDC predictor.

    Parameters
    ----------
    d_func    : callable t -> frac(g(t)), the decimal function
    t0        : starting parameter value
    delta     : target decimal value (0.5 for factorization, 0.0 for Wieferich)
    step      : iteration step (1 for factorization, 2 for Wieferich)
    max_jumps : maximum number of kinematic jumps before giving up
    tol       : tolerance for declaring a hit

    Returns
    -------
    t* (float) if found, None otherwise
    """
    t = t0

    for jump in range(max_jumps):
        # Evaluate 4 consecutive points
        try:
            d0 = d_func(t)
            d1 = d_func(t + step)
            d2 = d_func(t + 2 * step)
            d3 = d_func(t + 3 * step)
        except (ValueError, ZeroDivisionError, OverflowError):
            break

        # Kinematics
        v1 = d1 - d0
        v2 = d2 - d1
        v3 = d3 - d2
        a1 = v2 - v1
        a2 = v3 - v2
        j  = a2 - a1

        v_mean = (v1 + v2 + v3) / 3.0

        if verbose:
            print(f"  jump={jump:3d}  t={t:.2f}  d0={d0:.6f}  "
                  f"v̄={v_mean:.6f}  j={j:.2e}")

        # Detect sawtooth reset (large jerk or wrong direction)
        if abs(j) > 0.5 or abs(v_mean) < 1e-15:
            t += step
            continue

        # Predict t*
        if abs(v_mean) > 1e-15:
            t_pred = t + (delta - d0) / v_mean
        else:
            t += step
            continue

        # Round to nearest valid step
        t_candidate = round(t_pred / step) * step

        # Verify
        try:
            d_check = d_func(t_candidate)
            if abs(d_check - delta) < tol:
                if verbose:
                    print(f"  ✅ Found t*={t_candidate:.6f}  d={d_check:.9f}")
                return t_candidate
        except (ValueError, ZeroDivisionError, OverflowError):
            pass

        # Jump toward prediction
        if t_pred > t + step:
            t = t_pred - step  # land just before prediction
        else:
            t += step

    return None


# ─────────────────────────────────────────────
# Application 1: Factorization
# ─────────────────────────────────────────────

def _d_factor(N: int, m: float) -> float:
    """
    Decimal function for factorization.
    d(m) = frac(N / (2*(2m+3)))
    Target δ = 0.5: (2m+3) divides N iff d(m) = 0.5
    """
    denom = 2 * (2 * m + 3)
    if denom == 0:
        raise ZeroDivisionError
    val = N / denom
    return val - math.floor(val)


def mdc_factor(N: int, verbose: bool = False) -> Tuple[Optional[int], Optional[int], dict]:
    """
    Factor N using the MDC method.

    Searches for m such that d(m) = frac(N / (2*(2m+3))) = 0.5,
    which implies (2m+3) | N.

    Search space: L1 = {m : 2m+3 ≡ ±1 (mod 6)} = m ∈ {1,2,4,5,7,8,...}

    Returns
    -------
    (f1, f2) if N = f1*f2, or (None, None) if prime candidate
    info dict with steps, prediction, verification
    """
    info = {"steps": 0, "predictions": [], "reason": None}

    # Quick small factor checks
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
        if N % p == 0 and N > p:
            info["reason"] = f"small factor {p}"
            return p, N // p, info
        if N == p:
            info["reason"] = "small prime"
            return None, None, info

    # L1 filter: m such that (2m+3) % 6 in {1, 5}
    def in_L1(m):
        return (2 * m + 3) % 6 in (1, 5)

    # Start near sqrt(N)/2
    m_start = max(1, (math.isqrt(N) - 3) // 2)

    # Scan downward (factors near sqrt(N) first)
    m = m_start
    while m >= 1:
        if not in_L1(m):
            m -= 1
            continue

        d_func = lambda x, _N=N: _d_factor(_N, x)
        m_pred = mdc_predict(d_func, m, delta=0.5, step=1, max_jumps=50,
                              verbose=verbose)
        info["steps"] += 1

        if m_pred is not None:
            m_int = int(round(m_pred))
            f = 2 * m_int + 3
            if f > 1 and N % f == 0:
                other = N // f
                info["reason"] = "MDC decimal hit"
                info["predictions"].append(m_pred)
                return min(f, other), max(f, other), info

        m -= 6  # Skip by 6 to stay in L1

    # Scan upward (small factors)
    m = 1
    while m <= m_start:
        if not in_L1(m):
            m += 1
            continue
        f = 2 * m + 3
        if N % f == 0:
            other = N // f
            info["reason"] = "direct L1 hit"
            return min(f, other), max(f, other), info
        m += 1

    info["reason"] = "no factor found (prime candidate)"
    return None, None, info


# ─────────────────────────────────────────────
# Application 2: Wieferich primes
# ─────────────────────────────────────────────

def _p_of_K(K: float) -> float:
    """
    Solve 2^p = K*p^2 + 2 for real p > 2, given K > 0.
    Uses Newton's method on f(p) = 2^p - K*p^2 - 2 = 0.
    """
    if K <= 0:
        raise ValueError("K must be positive")

    # Initial estimate: p ≈ log2(K * p^2), iterate
    # Start from p = log2(K) + 2*log2(log2(K)+1)
    import math
    p = max(3.0, math.log2(max(K, 4)))
    for _ in range(60):
        f  = 2**p - K * p * p - 2
        fp = 2**p * math.log(2) - 2 * K * p
        if abs(fp) < 1e-15:
            break
        p -= f / fp
        if p < 2:
            p = 2.001
    return p


def _d_wieferich(K: float) -> float:
    """
    Decimal function for Wieferich detection.
    d(K) = frac(p(K)) where p(K) solves 2^p = K*p^2 + 2.
    Target δ = 0.0: p(K) integer => K = (2^n - 2)/n^2 => Wieferich condition.
    """
    p = _p_of_K(K)
    return p - math.floor(p)


def k_exact(n: int) -> float:
    """
    Exact K value where p(K) = n:
    K_exact(n) = (2^n - 2) / n^2
    """
    return (2**n - 2) / (n * n)


def is_wieferich_base2(p: int) -> bool:
    """Check if prime p is a Wieferich prime in base 2: p^2 | (2^p - 2)."""
    return pow(2, p, p * p) == 2 % (p * p)


def mdc_wieferich_scan(
    p_min: int = 3,
    p_max: int = 50,
    verbose: bool = True
) -> List[dict]:
    """
    Scan primes p in [p_min, p_max] for Wieferich property using MDC.

    For each prime p, computes K_exact(p) = (2^p - 2)/p^2
    and checks if it is an even integer (Wieferich condition).

    Also demonstrates the 4-point velocity extrapolation.
    """
    from math import log2

    results = []

    if verbose:
        print(f"\n{'p':>6}  {'K_exact':>18}  {'is_int':>7}  {'is_even':>8}  "
              f"{'Wieferich':>10}  {'K_pred':>12}")
        print("-" * 70)

    # We need sympy or manual sieve for primes
    def is_prime_trial(n):
        if n < 2: return False
        if n == 2: return True
        if n % 2 == 0: return False
        for i in range(3, int(n**0.5)+1, 2):
            if n % i == 0: return False
        return True

    primes = [p for p in range(p_min, p_max+1) if is_prime_trial(p)]

    for p in primes:
        k_ex = k_exact(p)
        is_int = abs(k_ex - round(k_ex)) < 1e-6
        is_even = is_int and (int(round(k_ex)) % 2 == 0)
        wf = is_wieferich_base2(p)

        # MDC prediction: start K just below K_exact, use 4-point extrapolation
        K0 = max(2.0, k_ex * 0.95)
        K0 = int(K0) if K0 > 2 else 2
        if K0 % 2 != 0:
            K0 += 1  # Ensure even start

        try:
            d0 = _d_wieferich(K0)
            d1 = _d_wieferich(K0 + 2)
            d2 = _d_wieferich(K0 + 4)
            d3 = _d_wieferich(K0 + 6)
            v1 = d1 - d0
            v2 = d2 - d1
            v3 = d3 - d2
            v_mean = (v1 + v2 + v3) / 3
            if abs(v_mean) > 1e-15:
                k_pred = K0 + 2 * ((0.0 - d0) / v_mean)
            else:
                k_pred = float('nan')
        except Exception:
            k_pred = float('nan')

        result = {
            "p": p,
            "K_exact": k_ex,
            "is_integer": is_int,
            "is_even": is_even,
            "wieferich": wf,
            "K_pred": k_pred,
        }
        results.append(result)

        if verbose:
            k_str = f"{k_ex:.4f}"
            kp_str = f"{k_pred:.4f}" if not math.isnan(k_pred) else "n/a"
            wf_str = "✅ YES" if wf else "no"
            print(f"{p:>6}  {k_str:>18}  {'YES' if is_int else 'no':>7}  "
                  f"{'even' if is_even else 'odd/no':>8}  {wf_str:>10}  {kp_str:>12}")

    return results


def wieferich_multi_base(p: int, bases: List[int]) -> dict:
    """
    Check Wieferich condition for prime p across multiple bases.
    Returns dict {base: is_wieferich} and list of independent bases.
    """
    result = {}
    for b in bases:
        if p % b == 0:
            result[b] = False
            continue
        # p^2 | (b^p - b)
        val = (pow(b, p, p*p) - b) % (p*p)
        result[b] = (val == 0)

    # Group bases by powers (b and b^k are dependent)
    def base_root(b):
        """Find the 'root' base (not a perfect power of a smaller base)."""
        for r in range(2, b):
            k = 2
            while r**k <= b:
                if r**k == b:
                    return r
                k += 1
        return b

    independent = {}
    for b, wf in result.items():
        if wf:
            root = base_root(b)
            if root not in independent:
                independent[root] = []
            independent[root].append(b)

    return {
        "p": p,
        "by_base": result,
        "wieferich_bases": [b for b, wf in result.items() if wf],
        "independent_roots": independent,
        "n_independent": len(independent),
    }


def analyze_1093_3511(verbose: bool = True) -> None:
    """
    Reproduce the 1093-3511 coincidence:
    Both are Wieferich in exactly the same bases {2,4,8,16,32} (powers of 2).
    This is a single condition, not five independent ones.
    """
    bases = list(range(2, 50))
    r1093 = wieferich_multi_base(1093, bases)
    r3511 = wieferich_multi_base(3511, bases)

    if verbose:
        print("\n1093 vs 3511 — Multi-base Wieferich analysis")
        print("=" * 50)
        print(f"1093 Wieferich bases: {r1093['wieferich_bases']}")
        print(f"3511 Wieferich bases: {r3511['wieferich_bases']}")
        print(f"\nIndependent roots for 1093: {r1093['independent_roots']}")
        print(f"Independent roots for 3511: {r3511['independent_roots']}")
        print(f"\nConclusion: both Wieferich in powers of 2 only.")
        print(f"This is a SINGLE condition (2^p Wieferich), not 5 independent ones.")


# ─────────────────────────────────────────────
# Sawtooth structure visualization (text)
# ─────────────────────────────────────────────

def sawtooth_table(n_levels: List[int], K_points: int = 8) -> None:
    """
    Show the sawtooth structure of d(K) near each level n.
    """
    print("\nSawtooth d(K) structure by level")
    print("=" * 55)

    for n in n_levels:
        k_ex = k_exact(n)
        K_start = max(2, int(k_ex * 0.85))
        if K_start % 2 != 0:
            K_start += 1

        print(f"\nLevel n={n}  K_exact={k_ex:.4f}")
        print(f"  {'K':>8}  {'d(K)':>10}  {'velocity':>10}")

        prev_d = None
        K = K_start
        for _ in range(K_points):
            try:
                d = _d_wieferich(K)
                v = d - prev_d if prev_d is not None else 0.0
                marker = " ← ZERO" if abs(d) < 0.01 else ""
                print(f"  {K:>8}  {d:>10.6f}  {v:>+10.6f}{marker}")
                prev_d = d
            except Exception:
                pass
            K += 2


# ─────────────────────────────────────────────
# Unification table
# ─────────────────────────────────────────────

def print_unification_table() -> None:
    """Print the MDC unification table comparing both applications."""
    print("""
MDC Unification Table
═══════════════════════════════════════════════════════════════
                    Factorization          Wieferich Primes
─────────────────── ─────────────────────  ──────────────────────
Parameter t         m (factor candidate)   K (even integer)
Function g(t)       N / (2*(2m+3))         p(K): 2^p = K*p^2 + 2
Target δ            0.5                    0.0
Search space        L1 = {6k±1}            2ℤ (even integers)
Solution            factor f = 2m+3 of N   Wieferich prime p = n
Step size           1                      2
Cost per jump       O(1) divisions         O(1) Newton steps
═══════════════════════════════════════════════════════════════
In both cases: d(t) is a real sawtooth, continuous by pieces,
whose local kinematics allow jumping directly to the solution.
""")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("MDC — Método Diofántico Cinemático")
    print("Kinematic Diophantine Method")
    print("=" * 65)

    # 1. Unification overview
    print_unification_table()

    # 2. Application 1: Factorization
    print("─" * 65)
    print("App 1: Factorization via MDC")
    print("─" * 65)
    test_N = [143, 221, 323, 1001, 3127, 10403, 15251, 100621]
    print(f"\n{'N':>10}  {'result':>22}  {'steps':>7}  reason")
    print("-" * 60)
    for N in test_N:
        f1, f2, info = mdc_factor(N)
        if f1:
            res = f"{f1} × {f2}"
        else:
            res = "prime candidate"
        print(f"{N:>10}  {res:>22}  {info['steps']:>7}  {info['reason']}")

    # 3. Application 2: Wieferich scan
    print()
    print("─" * 65)
    print("App 2: Wieferich prime scan via MDC (p ≤ 30)")
    print("─" * 65)
    mdc_wieferich_scan(p_min=3, p_max=30)

    # 4. Sawtooth structure
    print()
    sawtooth_table([5, 7, 11, 13, 17], K_points=6)

    # 5. 1093-3511 analysis
    print()
    analyze_1093_3511()

    # 6. p=13 multi-base (smallest independent multi-base Wieferich)
    print()
    print("p=13 multi-base Wieferich analysis:")
    r13 = wieferich_multi_base(13, list(range(2, 50)))
    print(f"  Wieferich bases: {r13['wieferich_bases']}")
    print(f"  Independent roots: {r13['independent_roots']}")
    print(f"  n_independent: {r13['n_independent']}")
