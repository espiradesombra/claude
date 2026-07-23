"""
siguiente_primo.py
==================
Next prime after a known integer — VMA / rescat 2026-07-23

Motor fiable (per defecte): roda 6k±1 + prova de divisors.
  - Correcte per a qualsevol inicio ≥ 1
  - Sense dependències

Mode experimental (method="karnaugh"):
  Algorisme original VMA amb acumuladors i mapa de Karnaugh.
  Es valida sempre amb _is_prime_ref; si falla, es fa fallback a la roda.

Autor material: Víctor Manzanares Alberola
Fix rescat: Grok (xAI) — l'heurística Karnaugh sola retornava falsos
(ex. siguiente_primo(10)=14). Ara el camí per defecte és correcte.
"""

from __future__ import annotations

import math
from typing import Iterator, List, Optional, Tuple


# ─────────────────────────────────────────────
# Reference primality
# ─────────────────────────────────────────────

def _is_prime_ref(n: int) -> bool:
    """Trial division 6k±1 — referència i motor fiable."""
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


# ─────────────────────────────────────────────
# Motor fiable: roda 6k±1
# ─────────────────────────────────────────────

def siguiente_primo_wheel(inicio: int) -> int:
    """
    Següent prim estrictament major que `inicio`.
    Recorre només candidats 6k±1 (després de 2,3).
    """
    if inicio < 2:
        return 2
    n = inicio + 1
    if n <= 2:
        return 2
    if n == 3:
        return 3
    # alinear a candidat 6k±1
    r = n % 6
    if r == 0:
        n += 1  # 6k → 6k+1
    elif r == 2:
        n += 3  # 6k+2 → 6k+5
    elif r == 3:
        n += 2  # 6k+3 → 6k+5
    elif r == 4:
        n += 1  # 6k+4 → 6k+5
    # r in {1,5}: ja és candidat

    # alternar +2 / +4 per recórrer 6k±1
    # si n ≡ 5 (mod 6): seqüència n, n+2, n+6, n+8, ... → passos 2,4,2,4...
    # si n ≡ 1 (mod 6): seqüència n, n+4, n+6, n+10, ... → passos 4,2,4,2...
    step = 2 if (n % 6 == 5) else 4
    while True:
        if _is_prime_ref(n):
            return n
        n += step
        step = 6 - step


# ─────────────────────────────────────────────
# Experimental Karnaugh (VMA original, validat)
# ─────────────────────────────────────────────

def siguiente_primo_karnaugh(inicio: int, verbose: bool = False) -> Optional[int]:
    """
    Intent de detecció per acumuladors + Karnaugh (pseudocodi VMA).
    Pot fallar; el caller ha de validar.
    """
    if inicio < 2:
        inicio = 2

    n = inicio
    ny = n - 1
    m = n + 1

    t = ny * ny
    tt = n * n
    nt = m * m

    n += 1
    m += 1
    ny += 1

    antp1 = False
    ant2p1 = False
    antp2 = False

    safety = max(inicio * 50 + 2000, 10_000)

    for step in range(safety):
        if ny <= 0 or n <= 0 or m <= 0:
            break
        t1 = t % ny
        t2 = t % n
        t3 = t % m
        tt2 = tt % n
        nt1 = nt % ny
        nt2 = nt % n

        core1 = (t3 > 0) and (nt2 == 0)
        paso1 = core1 or ant2p1
        new_ant2p1 = core1 or antp1
        new_antp1 = core1

        core2 = paso1 and (t2 > 0) and (tt2 > 0) and (t3 == 0)
        new_antp2 = core2
        paso3 = antp2 and (t1 > 0) and (nt1 + nt2 == 0)

        if verbose:
            print(f"step={step} n={n} paso3={paso3} detect_candidate={n-1}")

        if paso3:
            detected = n - 1
            if detected > inicio and _is_prime_ref(detected):
                return detected

        # creixement controlat (evita overflow boig en Python és ok, però
        # per velocitat reiniciem mòdul si cal — aquí productes complets)
        t *= ny
        tt *= n
        nt *= m

        ny += 1
        n += 1
        m += 1

        antp1 = new_antp1
        ant2p1 = new_ant2p1
        antp2 = new_antp2

    return None


# ─────────────────────────────────────────────
# API pública
# ─────────────────────────────────────────────

def siguiente_primo(
    inicio: int,
    verbose: bool = False,
    method: str = "wheel",
) -> Optional[int]:
    """
    Troba el següent prim després de `inicio`.

    method:
      "wheel"    — roda 6k±1 (defecte, sempre correcte)
      "karnaugh" — experimental VMA; si falla, fallback a wheel
      "auto"     — prova karnaugh, valida, fallback
    """
    if method == "wheel":
        return siguiente_primo_wheel(inicio)

    if method in ("karnaugh", "auto"):
        exp = siguiente_primo_karnaugh(inicio, verbose=verbose)
        if exp is not None and exp > inicio and _is_prime_ref(exp):
            # comprovar que no n'hem saltat cap
            expected = siguiente_primo_wheel(inicio)
            if exp == expected:
                return exp
            if verbose:
                print(f"karnaugh got {exp}, expected {expected} — fallback wheel")
        elif verbose:
            print("karnaugh miss — fallback wheel")
        return siguiente_primo_wheel(inicio)

    raise ValueError(f"method desconegut: {method}")


# Compat: nom antic heurística (mini annex)
def siguiente_primo_heuristica(p0: int) -> bool:
    """
    Heurística booleana de l'annex curt (una sola iteració).
    Conservada per compatibilitat; NO és un generador de prims.
    """
    n = p0
    m = n + 1
    ny = n - 1
    if ny <= 0:
        return False
    t = ny ** 2
    tt = n ** 2
    nt = m ** 2
    restos = [
        t % ny, t % n, t % m,
        tt % ny, tt % n, tt % m,
        nt % ny, nt % n, nt % m,
    ]
    bools = [r == 0 for r in restos]
    paso1 = bools[0] or bools[3]
    paso2 = bools[1] and not bools[4]
    paso3 = bools[2] ^ bools[5]
    return paso1 and paso2 and paso3


# ─────────────────────────────────────────────
# Generator / batch
# ─────────────────────────────────────────────

def prime_generator(start: int = 2) -> Iterator[int]:
    """Generador infinit de prims des de `start` (inclusiu si és prim)."""
    if start <= 2:
        yield 2
        current = 2
    else:
        current = start if _is_prime_ref(start) else siguiente_primo_wheel(start - 1)
        if current < start:
            current = siguiente_primo_wheel(current)
        # si start no és prim, current ja és el primer >= start
        if current < start:
            current = siguiente_primo_wheel(start - 1)
        yield current

    while True:
        current = siguiente_primo_wheel(current)
        yield current


def primeros_n_primos(n: int, start: int = 2) -> List[int]:
    """Primers n prims a partir de start (inclusiu si és prim)."""
    result: List[int] = []
    for p in prime_generator(start):
        result.append(p)
        if len(result) >= n:
            break
    return result


def validate(n_primes: int = 50, start: int = 2, verbose: bool = True) -> bool:
    """Valida siguiente_primo (wheel) vs trial division."""
    ref: List[int] = []
    c = 2
    while len(ref) < n_primes + 1:
        if _is_prime_ref(c):
            ref.append(c)
        c += 1

    if verbose:
        print(f"Validating siguiente_primo() for first {n_primes} primes...")
        print(f"{'p_ref':>8}  {'next_ref':>10}  {'next_alg':>10}  {'match':>6}")
        print("-" * 40)

    all_ok = True
    for i in range(min(n_primes, len(ref) - 1)):
        p = ref[i]
        expected = ref[i + 1]
        computed = siguiente_primo(p)
        match = computed == expected
        if verbose:
            status = "OK" if match else f"FAIL (got {computed})"
            print(f"{p:>8}  {expected:>10}  {str(computed):>10}  {status}")
        if not match:
            all_ok = False
            if verbose:
                print(f"  MISMATCH at p={p}: expected {expected}, got {computed}")
                break

    if verbose:
        print("\nAll OK." if all_ok else "\nMismatches found.")
    return all_ok


def find_twin_primes(limit: int) -> List[Tuple[int, int]]:
    twins = []
    prev = None
    for p in prime_generator(3):
        if p > limit:
            break
        if prev is not None and p - prev == 2:
            twins.append((prev, p))
        prev = p
    return twins


def prime_gap_analysis(start: int, count: int) -> None:
    primes = primeros_n_primos(count + 1, start)
    print(f"\nPrime gap analysis starting from {start}")
    print(f"{'p':>10}  {'next_p':>10}  {'gap':>6}")
    print("-" * 32)
    max_gap = 0
    max_pair = (0, 0)
    for i in range(len(primes) - 1):
        gap = primes[i + 1] - primes[i]
        if gap > max_gap:
            max_gap = gap
            max_pair = (primes[i], primes[i + 1])
        print(f"{primes[i]:>10}  {primes[i+1]:>10}  {gap:>6}")
    print(f"\nMax gap: {max_gap} between {max_pair}")


def wilson_comparison(primes: List[int]) -> None:
    print("\nWilson's theorem vs siguiente_primo()")
    print(f"{'p':>8}  {'(p-1)!%p':>12}  {'Wilson?':>8}  {'detected?':>10}")
    print("-" * 50)
    for p in primes:
        fac = 1
        for i in range(1, p):
            fac = (fac * i) % p
        wilson = fac == p - 1
        prev = p - 1
        while prev > 1 and not _is_prime_ref(prev):
            prev -= 1
        detected = siguiente_primo(prev) == p if prev >= 2 else (p == 2)
        print(f"{p:>8}  {fac:>12}  {str(wilson):>8}  {str(detected):>10}")


if __name__ == "__main__":
    print("=" * 65)
    print("siguiente_primo — wheel 6k±1 (fiable) + karnaugh experimental")
    print("=" * 65)
    print()
    validate(n_primes=40, start=2)
    print()
    print("First 20 primes:", primeros_n_primos(20))
    print("siguiente_primo(10) =", siguiente_primo(10))
    print("siguiente_primo(2)  =", siguiente_primo(2))
    print("twins ≤ 100:", find_twin_primes(100))
