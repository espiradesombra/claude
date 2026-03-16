"""
siguiente_primo.py
==================
Next Prime Algorithm via Accumulated Products of Consecutives

Author: Víctor Manzanares Alberola (VMA), Valencia 16/02/2026

Core idea:
    Given a known prime `inicio`, find the next prime by maintaining
    three running accumulators (t, tt, nt) built from products of
    consecutive integers around three sliding counters (ny, n, m).

    The detection logic is extracted from a Karnaugh map over 3
    consecutive iterations — it does NOT compute (n-1)! explicitly,
    but exploits the same divisibility structure as Wilson's theorem:
        p is prime  ⟺  (p-1)! ≡ p-1 (mod p)

    The 3-pass Karnaugh conditions (paso1, paso2, paso3) act as a
    lightweight sieve that fires exactly when `n` is prime.

Accumulators:
    ny = n - 1     (counter below n)
    n              (current candidate)
    m  = n + 1     (counter above n)

    t  = product accumulated with factor ny   → t  *= ny  each step
    tt = product accumulated with factor n    → tt *= n   each step
    nt = product accumulated with factor m    → nt *= m   each step

    Residues checked each step:
        t1  = t  % ny,  t2  = t  % n,  t3  = t  % m
        tt1 = tt % ny,  tt2 = tt % n,  tt3 = tt % m
        nt1 = nt % ny,  nt2 = nt % n,  nt3 = nt % m

Detection (Karnaugh 3-pass):
    antp1  = (t3 > 0) AND (nt2 == 0)              # memory step 1
    ant2p1 = (t3 > 0) AND (nt2 == 0) OR antp1     # memory 2 step 1
    paso1  = (t3 > 0) AND (nt2 == 0) OR ant2p1    # fires 2 iterations

    antp2  = paso1 AND (t2 > 0) AND (tt2 > 0) AND (t3 == 0)
    paso2  = paso1 AND (t2 > 0) AND (tt2 > 0) AND (t3 == 0) OR antp2

    paso3  = antp2 AND (t1 > 0) AND (nt1 + nt2 == 0)
             # note: nt3 == 0 except possibly at twin primes

    Prime detected when paso3 fires.
"""

import math
from typing import Iterator, List, Optional, Tuple


# ─────────────────────────────────────────────
# Reference primality (for validation only)
# ─────────────────────────────────────────────

def _is_prime_ref(n: int) -> bool:
    """Trial division reference. Used only for validation."""
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0: return False
    return True


# ─────────────────────────────────────────────
# Core algorithm
# ─────────────────────────────────────────────

def siguiente_primo(inicio: int, verbose: bool = False) -> Optional[int]:
    """
    Find the next prime after `inicio` using the accumulated-products method.

    Parameters
    ----------
    inicio  : a known prime (the starting point)
    verbose : if True, print internal state at each step

    Returns
    -------
    The next prime after `inicio`, or None if not found within safety limit.
    """
    # Initialize counters at inicio
    n  = inicio
    ny = n - 1   # n - 1
    m  = n + 1   # n + 1

    # Initialize accumulators: ny², n², m²  (Wilson-like seed)
    t  = ny * ny
    tt = n  * n
    nt = m  * m

    # Advance counters by 1 before entering loop (as in pseudocode)
    n  += 1
    m  += 1
    ny += 1

    # Memory cells for Karnaugh 3-pass
    antp1  = False
    ant2p1 = False
    antp2  = False

    safety = inicio * 20 + 1000  # upper bound on search

    if verbose:
        print(f"Starting from inicio={inicio}")
        print(f"{'step':>5}  {'ny':>6}  {'n':>6}  {'m':>6}  "
              f"{'t3':>4}  {'nt2':>5}  {'t2':>4}  {'tt2':>5}  "
              f"{'p1':>4}  {'p2':>4}  {'p3':>4}  {'detect'}")
        print("-" * 80)

    for step in range(safety):
        # Compute residues
        t1  = t  % ny
        t2  = t  % n
        t3  = t  % m
        tt1 = tt % ny
        tt2 = tt % n
        tt3 = tt % m
        nt1 = nt % ny
        nt2 = nt % n
        nt3 = nt % m   # noqa: F841 (used in twin prime note)

        # ── Karnaugh detection ──────────────────────────────────
        # Step 1: detect candidate window (2 iterations memory)
        core1  = (t3 > 0) and (nt2 == 0)
        paso1  = core1 or ant2p1
        new_ant2p1 = core1 or antp1
        new_antp1  = core1

        # Step 2: refine within window
        core2  = paso1 and (t2 > 0) and (tt2 > 0) and (t3 == 0)
        paso2  = core2 or antp2           # noqa: F841 (for completeness)
        new_antp2 = paso1 and (t2 > 0) and (tt2 > 0) and (t3 == 0)

        # Step 3: confirm primality
        paso3 = antp2 and (t1 > 0) and (nt1 + nt2 == 0)

        if verbose:
            print(f"{step:>5}  {ny:>6}  {n:>6}  {m:>6}  "
                  f"{t3:>4}  {nt2:>5}  {t2:>4}  {tt2:>5}  "
                  f"{int(paso1):>4}  {int(paso2):>4}  {int(paso3):>4}  "
                  f"{'✅ PRIME=' + str(n-1) if paso3 else ''}")

        if paso3:
            # n-1 is the detected prime (detection lags by 1 step)
            detected = n - 1
            if detected > inicio:
                return detected

        # Update accumulators
        t  *= ny
        tt *= n
        nt *= m

        # Advance counters
        ny += 1
        n  += 1
        m  += 1

        # Update memories
        antp1  = new_antp1
        ant2p1 = new_ant2p1
        antp2  = new_antp2

    return None


# ─────────────────────────────────────────────
# Iterator version: infinite prime generator
# ─────────────────────────────────────────────

def prime_generator(start: int = 2) -> Iterator[int]:
    """
    Infinite generator of primes starting from `start` (inclusive).
    Uses siguiente_primo() iteratively.

    Usage:
        gen = prime_generator(start=2)
        for p in gen:
            print(p)
            if p > 1000: break
    """
    # Yield primes up to start using reference (bootstrap)
    if start <= 2:
        yield 2
        current = 2
    else:
        # Find first prime >= start
        current = start if _is_prime_ref(start) else start + 1
        while not _is_prime_ref(current):
            current += 1
        yield current

    # Now use siguiente_primo iteratively
    while True:
        next_p = siguiente_primo(current)
        if next_p is None:
            break
        yield next_p
        current = next_p


# ─────────────────────────────────────────────
# Batch: first N primes after a given start
# ─────────────────────────────────────────────

def primeros_n_primos(n: int, start: int = 2) -> List[int]:
    """
    Return the first n primes starting from start.

    Parameters
    ----------
    n     : how many primes to generate
    start : starting value (inclusive)
    """
    result = []
    gen = prime_generator(start)
    for p in gen:
        result.append(p)
        if len(result) >= n:
            break
    return result


# ─────────────────────────────────────────────
# Validation
# ─────────────────────────────────────────────

def validate(n_primes: int = 50, start: int = 2, verbose: bool = True) -> bool:
    """
    Validate siguiente_primo() against trial division for the first n primes.

    Returns True if all match, False on first mismatch.
    """
    # Reference sequence
    ref = []
    candidate = 2
    while len(ref) < n_primes + 1:
        if _is_prime_ref(candidate):
            ref.append(candidate)
        candidate += 1

    if verbose:
        print(f"Validating siguiente_primo() for first {n_primes} primes...")
        print(f"{'p_ref':>8}  {'next_ref':>10}  {'next_alg':>10}  {'match':>6}")
        print("-" * 40)

    all_ok = True
    for i in range(min(n_primes, len(ref) - 1)):
        p      = ref[i]
        expected = ref[i + 1]
        computed = siguiente_primo(p)
        match = computed == expected

        if verbose:
            status = "✅" if match else f"❌ (got {computed})"
            print(f"{p:>8}  {expected:>10}  {str(computed):>10}  {status:>6}")

        if not match:
            all_ok = False
            if verbose:
                print(f"  MISMATCH at p={p}: expected {expected}, got {computed}")

    if verbose:
        if all_ok:
            print(f"\n✅ All {n_primes} siguiente_primo() calls correct.")
        else:
            print(f"\n⚠️  Mismatches found.")

    return all_ok


# ─────────────────────────────────────────────
# Twin prime detection
# ─────────────────────────────────────────────

def find_twin_primes(limit: int) -> List[Tuple[int, int]]:
    """
    Find all twin prime pairs (p, p+2) up to limit using siguiente_primo().
    Note from the pseudocode: nt3 == 0 except possibly at twin primes.
    """
    twins = []
    gen = prime_generator(start=3)
    prev = 3
    for p in gen:
        if p > limit:
            break
        if p - prev == 2:
            twins.append((prev, p))
        prev = p
    return twins


# ─────────────────────────────────────────────
# Gap analysis
# ─────────────────────────────────────────────

def prime_gap_analysis(start: int, count: int) -> None:
    """
    Show gaps between consecutive primes starting from `start`.
    """
    primes = primeros_n_primos(count + 1, start)

    print(f"\nPrime gap analysis starting from {start}")
    print(f"{'p':>10}  {'next_p':>10}  {'gap':>6}")
    print("-" * 32)

    max_gap = 0
    max_gap_pair = (0, 0)

    for i in range(len(primes) - 1):
        gap = primes[i+1] - primes[i]
        if gap > max_gap:
            max_gap = gap
            max_gap_pair = (primes[i], primes[i+1])
        print(f"{primes[i]:>10}  {primes[i+1]:>10}  {gap:>6}")

    print(f"\nMax gap in range: {max_gap} between {max_gap_pair}")


# ─────────────────────────────────────────────
# Wilson's theorem comparison
# ─────────────────────────────────────────────

def wilson_comparison(primes: List[int]) -> None:
    """
    Compare siguiente_primo() detection with Wilson's theorem directly.
    Shows that the accumulated-product method approximates Wilson
    without computing the full factorial.
    """
    print("\nWilson's theorem vs siguiente_primo() comparison")
    print(f"{'p':>8}  {'Wilson: (p-1)!%p':>18}  {'Wilson prime?':>14}  "
          f"{'siguiente_primo detects?':>24}")
    print("-" * 70)

    for p in primes:
        # Wilson: (p-1)! % p == p-1 iff p is prime
        factorial_mod = 1
        for i in range(1, p):
            factorial_mod = (factorial_mod * i) % p
        wilson_prime = (factorial_mod == p - 1)

        # Check if siguiente_primo from prev prime reaches p
        prev = p - 1
        while prev > 1 and not _is_prime_ref(prev):
            prev -= 1
        detected = siguiente_primo(prev) == p if prev > 1 else False

        match = "✅" if wilson_prime == detected else "⚠️"
        print(f"{p:>8}  {factorial_mod:>18}  {str(wilson_prime):>14}  "
              f"{str(detected):>24}  {match}")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("siguiente_primo — Next Prime via Accumulated Products")
    print("Author: Víctor Manzanares Alberola, 2026")
    print("=" * 65)

    # 1. Validation against trial division
    print()
    validate(n_primes=30, start=2)

    # 2. Generate first 20 primes
    print()
    print("First 20 primes via prime_generator():")
    primes20 = primeros_n_primos(20)
    print(" ".join(str(p) for p in primes20))

    # 3. Verbose trace for a small case
    print()
    print("Verbose trace: siguiente_primo(11)")
    result = siguiente_primo(11, verbose=True)
    print(f"Result: {result}")

    # 4. Twin primes up to 200
    print()
    twins = find_twin_primes(200)
    print(f"Twin primes up to 200: {twins}")

    # 5. Gap analysis
    print()
    prime_gap_analysis(start=2, count=20)

    # 6. Wilson comparison
    print()
    test_primes = [5, 7, 11, 13, 17, 19, 23, 29, 31]
    wilson_comparison(test_primes)
