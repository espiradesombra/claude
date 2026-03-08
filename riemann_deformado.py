"""
riemann_deformado.py
====================
Two deformed Riemann prime estimators by Víctor Manzanares Alberola (VMA)

From: "Cribas, cotas y estructuras modulares de los primos" (2025)

─────────────────────────────────────────────────────────
BACKGROUND
─────────────────────────────────────────────────────────

The classical Riemann prime-counting approximation is:

    R(n) = Σ_{k=1}^{K}  μ(k)/k · Li(n^(1/k))

where μ is the Möbius function and Li is the logarithmic integral.

VMA proposed two deformations that move the Möbius sign and
the 1/k factor to different positions:

    R̂(n) = Σ_{k=1}^{K}  Li( μ(k) · n^(1/k) )      [VMA-hat]
    R̃(n) = Σ_{k=1}^{K}  Li( μ(k) · n^(1/(k+1)) )  [VMA-tilde]

Key observations (validated numerically by VMA):
  - R(n) classic:  converges in 10–20 iterations
  - R̂(n):         needs 30–50 iterations, more oscillatory
  - R̃(n):         smoother than R̂, converges a bit faster
  - At K=100, all three are within ±3 of π(n) at n=100,000

VMA's insight: "The /k is to avoid iterating so much.
                If I iterate the same, my models compete equally."

─────────────────────────────────────────────────────────
DENSITY MODEL  (Sigo en mis trece, §final)
─────────────────────────────────────────────────────────

VMA also derived a density-layer model:

    π(x) ≈ x · Σ_{n=0}^{N}  D₀/2ⁿ · σ(x, 10ⁿ, s)

where σ is a sigmoid transition and D₀ = 4/9 ≈ 0.444...

The constant D₀ = 4/9 emerges from:
    2.197802... × 0.445 ≈ 1
    0.445 ≈ 4/9
    2.197802... = 1/(4/9) = 9/4

Stop criterion (modular, most efficient):
    N such that 10^N ≤ x

─────────────────────────────────────────────────────────
"""

import math
from typing import List, Optional


# ─────────────────────────────────────────────
# Möbius function
# ─────────────────────────────────────────────

def mobius(n: int) -> int:
    """
    Möbius function μ(n):
      μ(1) = 1
      μ(n) = 0  if n has a squared prime factor
      μ(n) = (-1)^k  if n is a product of k distinct primes
    """
    if n == 1: return 1
    factors = 0
    d = 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            factors += 1
            if n % d == 0:
                return 0  # squared factor
        d += 1
    if n > 1:
        factors += 1
    return (-1) ** factors


# ─────────────────────────────────────────────
# Logarithmic integral Li(x)
# ─────────────────────────────────────────────

def Li(x: float, terms: int = 50) -> float:
    """
    Logarithmic integral Li(x) = ∫₂^x dt/ln(t)
    Computed via series expansion around x=e:
        Li(x) ≈ Ei(ln x)
    Uses the Ramanujan series for Ei.
    """
    if x <= 0:
        return float('nan')
    if x == 1:
        return float('-inf')
    if abs(x - 1.0) < 1e-10:
        return 0.0

    # Use the offset logarithmic integral: Li(x) = Ei(ln x)
    # Ei(t) = γ + ln|t| + Σ_{k=1}^∞ t^k/(k·k!)
    # where γ = 0.5772156649...
    t = math.log(x) if x > 0 else float('nan')

    if math.isnan(t):
        return float('nan')

    EULER_GAMMA = 0.5772156649015328606

    # Ei(t) series
    ei = EULER_GAMMA + math.log(abs(t)) if t != 0 else float('-inf')
    term = 1.0
    factorial = 1.0
    for k in range(1, terms + 1):
        factorial *= k
        term = t ** k / (k * factorial)
        ei += term
        if abs(term) < 1e-15 * abs(ei):
            break

    return ei


def Li_safe(x: float) -> float:
    """
    Li with safety for negative or complex-inducing arguments.
    When μ(k) = -1 and x > 0, argument becomes negative →
    handled by returning Li(-|x|) with sign convention.
    """
    if x > 1:
        return Li(x)
    elif x > 0:
        return Li(x)
    elif x < 0:
        # Branch for negative argument (VMA's "philosophical" complex branch)
        return -Li(-x) if -x > 1 else 0.0
    else:
        return 0.0


# ─────────────────────────────────────────────
# Classical Riemann R(n)
# ─────────────────────────────────────────────

def R_clasico(n: float, K: int = 50) -> float:
    """
    Classical Riemann formula:
        R(n) = Σ_{k=1}^{K}  μ(k)/k · Li(n^(1/k))
    """
    total = 0.0
    for k in range(1, K + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        arg = n ** (1.0 / k)
        li_val = Li_safe(arg)
        total += (mu / k) * li_val
    return total


# ─────────────────────────────────────────────
# VMA deformations
# ─────────────────────────────────────────────

def R_hat(n: float, K: int = 50) -> float:
    """
    VMA deformation R̂(n):
        R̂(n) = Σ_{k=1}^{K}  Li( μ(k) · n^(1/k) )

    Differs from classical: the Möbius sign enters the
    argument of Li instead of multiplying the whole term.
    More oscillatory, needs more iterations to converge.
    """
    total = 0.0
    for k in range(1, K + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        arg = mu * (n ** (1.0 / k))
        total += Li_safe(arg)
    return total


def R_tilde(n: float, K: int = 50) -> float:
    """
    VMA deformation R̃(n):
        R̃(n) = Σ_{k=1}^{K}  Li( μ(k) · n^(1/(k+1)) )

    Uses exponent 1/(k+1) instead of 1/k.
    Smoother convergence than R̂, the exponent grows slower
    so terms decrease faster with k.
    """
    total = 0.0
    for k in range(1, K + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        arg = mu * (n ** (1.0 / (k + 1)))
        total += Li_safe(arg)
    return total


# ─────────────────────────────────────────────
# VMA density-layer model
# ─────────────────────────────────────────────

D0 = 4 / 9  # ≈ 0.4444...  — VMA density constant
             # Origin: 2.1978... × 0.445 ≈ 1, 0.445 ≈ 4/9

def sigma_transition(x: float, center: float, s: float = 1e4) -> float:
    """
    Sigmoid transition: σ(x, c, s) = 1 / (1 + exp((x - c) / s))
    """
    z = (x - center) / s
    z = max(-500, min(500, z))
    return 1.0 / (1.0 + math.exp(z))


def pi_density_layers(x: float, s: float = None) -> float:
    """
    VMA density-layer estimator:
        π(x) ≈ x · Σ_{n=0}^{N}  D₀/2ⁿ · σ(x, 10ⁿ, s)

    Stop criterion (modular): N such that 10^N ≤ x
    s defaults to x / 10 for smooth transitions.

    The formula can be read as: at each order of magnitude
    (10^n), a new density layer D₀/2ⁿ activates via a
    sigmoid gate, contributing to the total prime density.
    """
    if x < 2:
        return 0.0
    if s is None:
        s = max(x / 10, 1.0)

    total = 0.0
    n = 0
    while 10 ** n <= x:
        layer = D0 / (2 ** n)
        center = 10.0 ** n
        total += layer * sigma_transition(x, center, s)
        n += 1

    return x * total


# ─────────────────────────────────────────────
# Comparison table
# ─────────────────────────────────────────────

def tabla_comparativa(n_values: List[int], K: int = 50,
                      include_density: bool = True) -> None:
    """
    Compare R(n), R̂(n), R̃(n), and density model vs π(n) real.
    """
    from math import floor

    # π(n) via sieve (for validation)
    def pi_sieve(n: int) -> int:
        if n < 2: return 0
        sieve = bytearray([1]) * (n + 1)
        sieve[0] = sieve[1] = 0
        p = 2
        while p * p <= n:
            if sieve[p]:
                for j in range(p * p, n + 1, p):
                    sieve[j] = 0
            p += 1
        return sum(sieve)

    print(f"\nComparative table  (K={K} iterations)")
    print(f"{'n':>10}  {'π(n)':>8}  {'R(n)cl':>10}  "
          f"{'R̂(n)':>10}  {'R̃(n)':>10}  {'layers':>10}")
    print("-" * 64)

    for n in n_values:
        pi_n = pi_sieve(n)
        r_cl = R_clasico(n, K)
        r_hat = R_hat(n, K)
        r_til = R_tilde(n, K)
        lay = pi_density_layers(n) if include_density else float('nan')

        def fmt(v, ref):
            if math.isnan(v): return f"{'n/a':>10}"
            err = v - ref
            return f"{v:>8.1f}({err:+.1f})"

        print(f"{n:>10}  {pi_n:>8}  "
              f"{fmt(r_cl, pi_n):>14}  "
              f"{fmt(r_hat, pi_n):>14}  "
              f"{fmt(r_til, pi_n):>14}  "
              f"{fmt(lay, pi_n):>14}")


def convergencia_K(n: int = 100_000,
                   K_values: List[int] = None) -> None:
    """
    Show how all three estimators converge as K increases.
    """
    if K_values is None:
        K_values = [5, 10, 20, 50, 100]

    from math import isqrt

    # π(n) via sieve
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    p = 2
    while p * p <= n:
        if sieve[p]:
            for j in range(p * p, n + 1, p):
                sieve[j] = 0
        p += 1
    pi_n = sum(sieve)

    print(f"\nConvergence with K  (n={n:,}, π(n)={pi_n})")
    print(f"{'K':>5}  {'R(n) cl':>10}  {'err':>6}  "
          f"{'R̂(n)':>10}  {'err':>6}  "
          f"{'R̃(n)':>10}  {'err':>6}")
    print("-" * 58)

    for K in K_values:
        rc = R_clasico(n, K)
        rh = R_hat(n, K)
        rt = R_tilde(n, K)
        print(f"{K:>5}  {rc:>10.1f}  {rc-pi_n:>+6.1f}  "
              f"{rh:>10.1f}  {rh-pi_n:>+6.1f}  "
              f"{rt:>10.1f}  {rt-pi_n:>+6.1f}")


def carrera_estimadores(n: int = 100_000, K: int = 50) -> None:
    """
    Race: which estimator is closest to π(n)? (VMA §carrera)
    """
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    p = 2
    while p * p <= n:
        if sieve[p]:
            for j in range(p * p, n + 1, p):
                sieve[j] = 0
        p += 1
    pi_n = sum(sieve)

    competitors = [
        ("π(n) real",    pi_n,              "🎯 Meta"),
        ("R(n) clásico", R_clasico(n, K),   "🥇"),
        ("R̃(n) tilde",  R_tilde(n, K),     "🥈"),
        ("R̂(n) hat",    R_hat(n, K),       "🥉"),
        ("Density D₀",   pi_density_layers(n), ""),
    ]

    print(f"\n🏁 Carrera de estimadores  (n={n:,}, K={K})")
    print(f"{'Modelo':20s}  {'resultado':>12}  {'error abs':>10}  {'pos'}")
    print("-" * 55)

    ranked = sorted(competitors[1:], key=lambda x: abs(x[1] - pi_n))
    print(f"{'π(n) real':20s}  {pi_n:>12}  {'—':>10}")
    for i, (name, val, medal) in enumerate(ranked):
        err = val - pi_n
        print(f"{name:20s}  {val:>12.1f}  {err:>+10.1f}  {medal}")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("Riemann Deformado — R̂(n) y R̃(n)")
    print("Author: Víctor Manzanares Alberola")
    print("=" * 65)

    print(f"\nD₀ = 4/9 = {D0:.6f}  (VMA density constant)")
    print(f"Origin: 9/4 = 2.25,  1/2.25 = {1/2.25:.6f}")

    # Comparative table
    n_vals = [1_000, 10_000, 100_000]
    tabla_comparativa(n_vals, K=50)

    # Convergence with K
    convergencia_K(n=100_000)

    # Race at n=100,000
    carrera_estimadores(n=100_000, K=50)

    # Quick single values
    print(f"\nQuick check Li(x):")
    for x in [10, 100, 1000]:
        from math import log
        li = Li(x)
        approx = x / log(x)
        print(f"  Li({x:6d}) = {li:10.4f}   x/ln(x) = {approx:10.4f}")
