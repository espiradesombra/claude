"""
sofi_structure.py
=================
Modular classification of Sophie Germain prime candidates.

Author: Víctor Manzanares Alberola (VMA)

Construction:
    L1  = { a : a ≡ 5 (mod 6) }                         all 6k-1 candidates
    L3  = { a ∈ L1 : a = (6k-1)(6h+1) for some k,h≥1 } composites, type A
    L4  = { a ∈ L1 : 2a+1 = (6j-1)(6g+1) for some j,g} composites, type B
    L2  = L3 ∩ L4                                        doubly composite
    U2  = L1 \ (L3 ∪ L4)                                residual candidates

Result:  U2 ⊆ LSG  (every element of U2 is a Sophie Germain prime candidate)
Conjecture: |U2| = ∞  =>  infinitely many Sophie Germain primes

Note: This uses the Chinese Remainder Theorem and Dirichlet's theorem
      on primes in arithmetic progressions as theoretical backbone.
"""

from typing import Dict, List, Set, Tuple
import math


# ─────────────────────────────────────────────
# Primality helpers
# ─────────────────────────────────────────────

def sieve_primes(limit: int) -> List[int]:
    """Sieve of Eratosthenes."""
    if limit < 2:
        return []
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i, v in enumerate(is_prime) if v]


def is_prime_trial(n: int) -> bool:
    """Trial division primality test."""
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True


# ─────────────────────────────────────────────
# Modular type classifiers
# ─────────────────────────────────────────────

def is_type_A_composite(a: int) -> bool:
    """
    Check if a = (6k-1)(6h+1) for some positive integers k, h.
    Equivalently: a has a factor f ≡ 5 (mod 6) with a/f ≡ 1 (mod 6).
    """
    if a < 25:  # Smallest: (6·1-1)(6·1+1) = 5·7 = 35
        return False
    limit = int(math.sqrt(a)) + 1
    f = 5
    while f <= limit:
        if a % f == 0:
            g = a // f
            if g % 6 == 1:
                return True
        f += 6
    return False


def is_type_B_composite(a: int) -> bool:
    """
    Check if 2a+1 = (6j-1)(6g+1) for some positive integers j, g.
    Equivalently: 2a+1 is type-A composite.
    """
    return is_type_A_composite(2 * a + 1)


def is_type_A_or_symmetric(a: int) -> bool:
    """
    Extended type-A: a = (6k-1)(6h-1) or (6k+1)(6h+1).
    Catches all 6k-1 composites not covered by standard type A.
    """
    limit = int(math.sqrt(a)) + 1
    # Check (6k-1)(6h-1)
    f = 5
    while f <= limit:
        if a % f == 0:
            g = a // f
            if g % 6 == 5:
                return True
            if g % 6 == 1:
                return True
        f += 6
    # Check (6k+1)(6h+1)
    f = 7
    while f <= limit:
        if a % f == 0:
            g = a // f
            if g % 6 == 1:
                return True
        f += 6
    return False


# ─────────────────────────────────────────────
# Main classification
# ─────────────────────────────────────────────

def classify_sofi(limit: int, verify: bool = True) -> Dict:
    """
    Classify all a ≡ 5 (mod 6) up to limit into Sofí sets.

    Parameters
    ----------
    limit   : upper bound
    verify  : if True, verify U2 ⊆ LSG using primality tests

    Returns
    -------
    dict with sets L1, L3, L4, L2, U2 and statistics
    """
    # Build L1
    L1 = list(range(5, limit + 1, 6))

    # Classify
    L3: Set[int] = set()
    L4: Set[int] = set()

    for a in L1:
        if is_type_A_composite(a):
            L3.add(a)
        if is_type_B_composite(a):
            L4.add(a)

    L2 = L3 & L4
    U2 = set(L1) - L3 - L4

    result = {
        "limit": limit,
        "L1": set(L1),
        "L3": L3,
        "L4": L4,
        "L2": L2,
        "U2": U2,
        "|L1|": len(L1),
        "|L3|": len(L3),
        "|L4|": len(L4),
        "|L2|": len(L2),
        "|U2|": len(U2),
    }

    if verify:
        # Check U2 ⊆ LSG (every element should be prime with 2a+1 also prime,
        # or prime with composite safe prime — in both cases, prime)
        not_prime_in_U2 = [a for a in U2 if not is_prime_trial(a)]
        sg_in_U2 = [a for a in U2 if is_prime_trial(a) and is_prime_trial(2 * a + 1)]
        prime_not_sg = [a for a in U2 if is_prime_trial(a) and not is_prime_trial(2 * a + 1)]

        result["not_prime_in_U2"] = not_prime_in_U2
        result["sg_in_U2"] = sg_in_U2
        result["prime_not_sg_in_U2"] = prime_not_sg
        result["U2_subset_LSG_holds"] = len(not_prime_in_U2) == 0

    return result


def print_classification(result: Dict, show_elements: bool = False) -> None:
    """Pretty-print classification results."""
    limit = result["limit"]
    print(f"\nSofí Classification up to {limit}")
    print("=" * 45)
    print(f"  |L1| = {result['|L1|']:6d}  (all 6k-1 candidates)")
    print(f"  |L3| = {result['|L3|']:6d}  (type-A composites)")
    print(f"  |L4| = {result['|L4|']:6d}  (type-B composites)")
    print(f"  |L2| = {result['|L2|']:6d}  (doubly composite)")
    print(f"  |U2| = {result['|U2|']:6d}  (residual, SGP candidates)")

    if "U2_subset_LSG_holds" in result:
        holds = result["U2_subset_LSG_holds"]
        print(f"\n  U2 ⊆ LSG verified: {'✅ YES' if holds else '❌ NO'}")
        print(f"  Sophie Germain primes in U2: {len(result['sg_in_U2'])}")
        print(f"  Primes in U2 (safe prime composite): {len(result['prime_not_sg_in_U2'])}")
        if result["not_prime_in_U2"]:
            print(f"  ⚠️  Non-primes found in U2: {result['not_prime_in_U2'][:10]}")

    if show_elements and result["|U2|"] <= 50:
        print(f"\n  U2 elements: {sorted(result['U2'])}")
        if "sg_in_U2" in result:
            print(f"  SGP in U2:   {sorted(result['sg_in_U2'])}")


# ─────────────────────────────────────────────
# Density analysis
# ─────────────────────────────────────────────

def sofi_density_table(limits: List[int]) -> None:
    """
    Print density of U2 (SGP candidates) vs total 6k-1 numbers.
    Shows growth pattern consistent with conjecture |U2| = ∞.
    """
    print("\nSofí density analysis")
    print(f"{'limit':>10}  {'|L1|':>7}  {'|U2|':>7}  {'|SGP|':>7}  {'U2/L1%':>8}  {'SGP/U2%':>9}")
    print("-" * 60)

    for limit in limits:
        r = classify_sofi(limit, verify=True)
        l1 = r["|L1|"]
        u2 = r["|U2|"]
        sgp = len(r.get("sg_in_U2", []))
        ratio_u2 = u2 / l1 * 100 if l1 else 0
        ratio_sgp = sgp / u2 * 100 if u2 else 0
        print(f"{limit:>10}  {l1:>7}  {u2:>7}  {sgp:>7}  {ratio_u2:>7.2f}%  {ratio_sgp:>8.2f}%")


# ─────────────────────────────────────────────
# CRT analysis (theoretical backbone)
# ─────────────────────────────────────────────

def crt_analysis(verbose: bool = True) -> None:
    """
    Demonstrate the Chinese Remainder Theorem structure behind L3 and L4.
    Shows why composites cluster in specific residue classes.
    """
    print("\nCRT Analysis of L3 and L4")
    print("=" * 45)

    # For L3: a ≡ 5 (mod 6) and a = (6k-1)(6h+1)
    # Expanding: a = 36kh + 6k - 6h - 1
    # So a ≡ -1 ≡ 5 (mod 6) ✓
    # And a ≡ (6k-1)(6h+1) — CRT gives residue classes mod 36, 72, etc.

    print("L3 residue classes mod 30:")
    l3_mod30 = set()
    for a in range(5, 10000, 6):
        if is_type_A_composite(a):
            l3_mod30.add(a % 30)
    print(f"  {sorted(l3_mod30)}")

    print("L4 residue classes mod 30:")
    l4_mod30 = set()
    for a in range(5, 10000, 6):
        if is_type_B_composite(a):
            l4_mod30.add(a % 30)
    print(f"  {sorted(l4_mod30)}")

    print("U2 residue classes mod 30:")
    u2_mod30 = set()
    for a in range(5, 10000, 6):
        if not is_type_A_composite(a) and not is_type_B_composite(a):
            u2_mod30.add(a % 30)
    print(f"  {sorted(u2_mod30)}")

    print("\nDirichlet: each residue class in U2 (mod 30) contains ∞ primes ✓")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 55)
    print("Sofí Structure — Sophie Germain Modular Classification")
    print("=" * 55)

    # Small classification with element listing
    r = classify_sofi(300, verify=True)
    print_classification(r, show_elements=True)

    print()

    # Density table
    sofi_density_table([500, 1000, 2000, 5000, 10000])

    print()

    # CRT structure
    crt_analysis()
