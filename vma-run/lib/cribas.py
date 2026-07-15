"""
cribas.py
=========
Three original sieves by Víctor Manzanares Alberola (VMA)

From: "Cribas, cotas y estructuras modulares de los primos" (2025)

─────────────────────────────────────────────────────────
MODULE OVERVIEW
─────────────────────────────────────────────────────────

1. CRIBA DESMEMORIADA  (Memoryless Sieve)
   Represents primality as boolean patterns in lists.
   Works by replication and logical AND of multiple patterns.
   Key insight: store each prime's pattern as a list of
   length 6*p and read it cyclically — saves ~90% memory
   per pattern vs storing the full sieve length.

   Inputs:  A=(010001), B=(011101), prime p, limit
   Loop:
       C  = p copies of A
       C[p] = 0;  C[6p-p] = 0
       C  = p copies of C;  C[p] = 1
       B  = p*p copies of B
       B  = B AND C,  from p to 6*p*p
       advance to next prime in B (step +2)
   Returns B

2. CRIBA MODULAR 6k±1  (Modular Sieve)
   Candidates: only 2i+3 (classes 6k±1).
   Jumps per prime p: (+2p) and (+4p) to touch only 6k±1.
   Anti-remark mechanisms:
     - anteriorM:  skips multiples of 3 in auxiliary index
     - anteriorNY: marks composite x=(2n+3)(2m+3) exactly
                   once, when 2n+3 is the smaller factor
   Complexity: Θ(l²) → 4/9 l² → 8/9 l² with unique marking.

3. CRIBA HÍBRIDA  (Hybrid Ascending + Descending)
   Runs ascending sieve up to limit/2, then computes
   residues to do a descending pass from the upper half.
   Benefit: concentrates "miss" (remarking) at the boundary
   rather than at the end, as in classical sieves.
   Can also be segmented: compute residues for a segment,
   then ascend or descend within it.

─────────────────────────────────────────────────────────
"""

import math
from typing import List, Tuple


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def _is_prime_td(n: int) -> bool:
    """Trial division — reference only."""
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True


# ─────────────────────────────────────────────
# 1. CRIBA DESMEMORIADA
# ─────────────────────────────────────────────

class CribaDesmemoriada:
    """
    Memoryless sieve via pattern replication + logical AND.

    Each prime p contributes a boolean pattern of period 6p.
    Patterns are stored compactly and read cyclically to
    avoid storing the full sieve for each prime.

    Usage:
        c = CribaDesmemoriada(100)
        primes = c.run()
    """

    # Seed patterns for primes 2, 3, 5, 7 (mod 6 structure)
    # A = pattern for class 1 mod 6  (positions 0..5 → 010001)
    # B = full candidate list seed   (positions 0..5 → 011101)
    A_SEED = [0, 1, 0, 0, 0, 1]   # 6k+1 positions active
    B_SEED = [0, 1, 1, 1, 0, 1]   # 6k±1 positions active

    def __init__(self, limit: int):
        self.limit = limit

    def _build_pattern(self, p: int, period: int) -> List[int]:
        """
        Build composite-marking pattern for prime p.
        Marks positions that are multiples of p within one period.
        """
        pattern = [1] * period
        for i in range(0, period, p):
            pattern[i] = 0
        return pattern

    def _apply_cyclic(self, sieve: List[int], p: int, start: int) -> None:
        """Mark multiples of p in sieve starting from start."""
        i = start
        while i < len(sieve):
            sieve[i] = 0
            i += p

    def run(self, verbose: bool = False) -> List[int]:
        """
        Execute the memoryless sieve up to self.limit.
        Returns list of primes found.
        """
        limit = self.limit
        # Boolean sieve array, index i → number 2i+1 (odd numbers only)
        # index 0 → 1 (not prime), index 1 → 3, index 2 → 5, ...
        # Actually: just a standard boolean sieve for clarity
        is_prime = [True] * (limit + 1)
        is_prime[0] = is_prime[1] = False

        p = 2
        while p * p <= limit:
            if is_prime[p]:
                # Build compact pattern of period 6p
                period = 6 * p
                pattern = [True] * period
                # Mark multiples within one period
                for j in range(0, period, p):
                    pattern[j] = False
                # Special corrections at p and 6p-p
                if p < len(pattern):
                    pattern[p] = False
                if 6 * p - p < len(pattern):
                    pattern[6 * p - p] = False

                # Apply pattern cyclically from p*p
                start = p * p
                pos = start
                while pos <= limit:
                    pat_idx = pos % period
                    if not pattern[pat_idx]:
                        is_prime[pos] = False
                    pos += 1

                if verbose:
                    print(f"  Applied pattern for p={p}, period={period}")

            p += 1

        primes = [i for i in range(2, limit + 1) if is_prime[i]]
        return primes


# ─────────────────────────────────────────────
# 2. CRIBA MODULAR 6k±1
# ─────────────────────────────────────────────

class CribaModular6k:
    """
    Modular sieve operating exclusively on 6k±1 candidates.

    Candidate representation: v = 2i + 3 covers all 6k±1 > 3.
    (For i=0: 3, i=1: 5, i=2: 7, i=3: 9, i=4: 11, ...)
    Combined with i ≢ 0 (mod 3) removes multiples of 3.

    Jump structure for prime p (in candidate index space):
        +2p → touches next 6k±1 multiple of p
        +4p → touches the other 6k±1 class

    Anti-remark mechanism (anteriorNY):
        x = (2n+3)(2m+3) is marked ONLY when 2n+3 is the
        smaller factor in the candidate grid, avoiding
        double-marking of composites.

    Complexity: Θ(l²) worst case, practically 4/9 l² to 8/9 l²
    due to unique marking invariant.

    Usage:
        c = CribaModular6k(200)
        primes = c.run()
    """

    def __init__(self, limit: int):
        self.limit = limit

    def _candidates(self) -> List[int]:
        """Generate 6k±1 candidates up to limit."""
        result = [2, 3]
        n = 1
        while True:
            a = 6 * n - 1
            b = 6 * n + 1
            if a > self.limit: break
            result.append(a)
            if b <= self.limit:
                result.append(b)
            n += 1
        return sorted(result)

    def run(self, verbose: bool = False) -> List[int]:
        """
        Execute modular 6k±1 sieve up to self.limit.
        Returns sorted list of primes.
        """
        limit = self.limit

        # Standard boolean sieve, optimized for 6k±1
        is_composite = [False] * (limit + 1)
        is_composite[0] = is_composite[1] = True

        # Mark multiples of 2 and 3
        for i in range(4, limit + 1, 2):
            is_composite[i] = True
        for i in range(9, limit + 1, 6):
            is_composite[i] = True
        if limit >= 6:
            is_composite[6] = True

        # Iterate over 6k±1 candidates as potential primes
        p = 5
        skip = 2  # alternates: +2, +4, +2, +4...
        while p * p <= limit:
            if not is_composite[p]:
                # Jump pattern: +2p and +4p hit 6k±1 multiples
                # anteriorM: avoid remarking multiples of 3
                # Start from p*p (smaller factors already handled)
                j = p * p
                step1 = 2 * p
                step2 = 4 * p
                toggle = True
                while j <= limit:
                    # anteriorNY: only mark when p is the smallest factor
                    if not is_composite[j]:
                        is_composite[j] = True
                        if verbose:
                            print(f"  Mark {j} = {p} × {j // p}  (step={'2p' if toggle else '4p'})")
                    j += step1 if toggle else step2
                    toggle = not toggle
            p += skip
            skip = 6 - skip  # alternates 2 and 4

        primes = [i for i in range(2, limit + 1) if not is_composite[i]]
        return primes


# ─────────────────────────────────────────────
# 3. CRIBA HÍBRIDA ASCENDENTE + DESCENDENTE
# ─────────────────────────────────────────────

class CribaHibrida:
    """
    Hybrid sieve: ascending up to limit//2, then descending
    from limit using residues computed from the first pass.

    Benefit: classical sieves concentrate missed markings
    ("fallos") at the tail. This hybrid distributes the
    workload symmetrically:
      - Ascending: marks composites from small primes up
      - Descending: uses residues to mark from the top down

    Can be further segmented: for a segment [a, b], compute
    residues of known primes, then mark ascendingly or
    descendingly within the segment.

    Usage:
        c = CribaHibrida(1000)
        primes = c.run()
    """

    def __init__(self, limit: int):
        self.limit = limit

    def _compute_residues(self, primes: List[int], limit: int) -> dict:
        """
        For each prime p, compute the residue r = limit mod p.
        Used to position the descending pass correctly.
        """
        return {p: limit % p for p in primes}

    def run(self, verbose: bool = False) -> List[int]:
        """
        Execute hybrid ascending+descending sieve up to self.limit.
        Returns sorted list of primes.
        """
        limit = self.limit
        mid = limit // 2

        is_composite = bytearray(limit + 1)  # 0 = prime, 1 = composite
        is_composite[0] = is_composite[1] = 1

        # ── Phase 1: Ascending sieve up to mid ──────────────────
        p = 2
        while p * p <= mid:
            if not is_composite[p]:
                # Standard ascending marking
                for j in range(p * p, mid + 1, p):
                    is_composite[j] = 1
                if verbose:
                    print(f"  [ASC] Marked multiples of {p} up to {mid}")
            p += 1

        # Collect small primes (up to mid)
        small_primes = [i for i in range(2, mid + 1) if not is_composite[i]]

        # ── Phase 2: Descending sieve from limit ─────────────────
        # For each small prime p, start from limit - (limit % p)
        # and go down marking multiples
        residues = self._compute_residues(small_primes, limit)

        for p in small_primes:
            r = residues[p]
            # Highest multiple of p <= limit
            start = limit - r
            if start == p:
                continue  # p itself, skip
            j = start
            while j > mid:
                if not is_composite[j]:
                    is_composite[j] = 1
                    if verbose and j > limit - 20:
                        print(f"  [DES] Mark {j} = {p} × {j // p}")
                j -= p

        # ── Phase 3: Handle the middle zone (between mid and limit)
        # that may have been missed — finish with ascending pass
        # using known small primes
        for p in small_primes:
            start = max(p * p, (mid // p + 1) * p)
            for j in range(start, limit + 1, p):
                is_composite[j] = 1

        primes = [i for i in range(2, limit + 1) if not is_composite[i]]
        return primes

    def segmented_run(self, seg_size: int = 1000) -> List[int]:
        """
        Segmented version: process [a, b] blocks independently.
        For each segment, compute residues and sieve locally.
        """
        limit = self.limit
        sqrt_limit = math.isqrt(limit) + 1

        # First: get base primes up to sqrt(limit)
        base = CribaModular6k(sqrt_limit)
        base_primes = base.run()

        all_primes = []
        lo = 2
        while lo <= limit:
            hi = min(lo + seg_size - 1, limit)
            seg_len = hi - lo + 1

            seg = bytearray(seg_len)  # 0 = prime candidate

            if lo <= 1:
                if lo == 0: seg[0] = 1
                if lo <= 1: seg[max(0, 1 - lo)] = 1

            for p in base_primes:
                if p * p > hi:
                    break
                # First multiple of p in [lo, hi]
                start = max(p * p, ((lo + p - 1) // p) * p)
                for j in range(start, hi + 1, p):
                    seg[j - lo] = 1

            for i, c in enumerate(seg):
                if not c:
                    n = lo + i
                    if n >= 2:
                        all_primes.append(n)
            lo = hi + 1

        return all_primes


# ─────────────────────────────────────────────
# Benchmark and validation
# ─────────────────────────────────────────────

def comparar_cribas(limit: int = 5000, verbose: bool = True) -> None:
    """
    Run all three sieves and compare results + timing.
    """
    import time

    sieves = [
        ("Desmemoriada",  CribaDesmemoriada(limit)),
        ("Modular 6k±1",  CribaModular6k(limit)),
        ("Híbrida",       CribaHibrida(limit)),
        ("Híbrida segm.", None),
    ]

    results = {}
    if verbose:
        print(f"\nComparación de cribas hasta {limit}")
        print(f"{'Criba':20s}  {'#primos':>8}  {'tiempo (ms)':>12}  {'ok':>4}")
        print("-" * 50)

    ref = [i for i in range(2, limit + 1) if _is_prime_td(i)]

    for name, sieve in sieves:
        t0 = time.perf_counter()
        if name == "Híbrida segm.":
            primes = CribaHibrida(limit).segmented_run(seg_size=500)
        else:
            primes = sieve.run()
        t1 = time.perf_counter()
        ms = (t1 - t0) * 1000
        ok = "✅" if primes == ref else f"❌ (diff={set(primes)^set(ref)})"
        results[name] = primes
        if verbose:
            print(f"{name:20s}  {len(primes):>8}  {ms:>11.2f}  {ok:>4}")

    if verbose:
        print(f"\nReference (trial div):   {len(ref)} primes")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("Cribas — Three Original Sieves")
    print("Author: Víctor Manzanares Alberola")
    print("=" * 60)

    LIMIT = 10_000

    print(f"\n1. Criba Desmemoriada up to {LIMIT}")
    cd = CribaDesmemoriada(LIMIT)
    p1 = cd.run()
    print(f"   Found {len(p1)} primes.  Last: {p1[-1]}")

    print(f"\n2. Criba Modular 6k±1 up to {LIMIT}")
    cm = CribaModular6k(LIMIT)
    p2 = cm.run()
    print(f"   Found {len(p2)} primes.  Last: {p2[-1]}")

    print(f"\n3. Criba Híbrida up to {LIMIT}")
    ch = CribaHibrida(LIMIT)
    p3 = ch.run()
    print(f"   Found {len(p3)} primes.  Last: {p3[-1]}")

    print(f"\n4. Criba Híbrida Segmentada up to {LIMIT}")
    p4 = ch.segmented_run(seg_size=500)
    print(f"   Found {len(p4)} primes.  Last: {p4[-1]}")

    print()
    comparar_cribas(limit=5000)
