# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC PARABOLA HYBRID v2
  Bisecció Adaptativa amb Densitat Logarítmica
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)
  VERSIÓ: Hybrid-2.0
  COL·LABORACIÓ D'ESCRIPTURA: Claude (Anthropic)

  ─────────────────────────────────────────────────────────────────────────
  NOVETAT CENTRAL v2: BISECCIÓ ADAPTATIVA (proposta de Víctor)
  ─────────────────────────────────────────────────────────────────────────

  v1 usava segments harmònics fixos: [m_max/1], [m_max/2], [m_max/3]...
  La densitat de mostres era uniforme en tots els segments.

  v2 usa BISECCIÓ ADAPTATIVA amb densitat variable:

  PER QUÈ TÉ SENTIT LA DISTRIBUCIÓ LOGARÍTMICA?
  ─────────────────────────────────────────────
  Si N = p × q amb p ≤ q, la probabilitat que p caiga en [m, 2m] és
  proporcional a la densitat de factors en eixe interval:

      densitat(m) ∝ 1/k = 1/(N//(2m+3)) ≈ (2m+3)/N

  Per tant hi ha MÉS factors en la zona baixa (m petit, factor menut)
  que en la zona alta. En particular:
      · m ∈ [1, m_max/2]       → aproximadament el 50% dels factors
      · m ∈ [m_max/2, m_max]   → aproximadament l'altre 50%
  (no és uniforme: la meitat inferior té molts més valors de k possibles)

  REGLA D'ADAPTACIÓ DE MOSTRES:
  ─────────────────────────────
  El nombre de mostres s'adapta a l'amplada del rang:
      amplada > LLINDAR_GRAN  → N_MOSTRA_BASE    (poques, molt espaiades)
      amplada ∈ [LLINDAR_MITJ, LLINDAR_GRAN] → N_MOSTRA_BASE × 3
      amplada < LLINDAR_MITJ  → k-sweep directe (exhaustiu i ràpid)

  PRIORITAT DE CERCA (clau per a la distribució logarítmica):
  ────────────────────────────────────────────────────────────
  La bisecció genera subrangs en aquest ordre de prioritat:
      1. La meitat INFERIOR del rang actual [m_ini, m_mig]
         (zona on la densitat logarítmica és MÀXIMA)
      2. La meitat SUPERIOR [m_mig, m_fi]
         (menys densa, però cal cobrir-la)

  Açò garanteix que els factors menuts (que estadísticament son més
  probables en semiprimers RSA-like) s'explorin primer.

  CONNEXIÓ MDC → DISCRIMINANT S±k (igual que v1, però en cada subrang):
  ──────────────────────────────────────────────────────────────────────
  Per cada "zona calenta" detectada per d(m) ≈ 0:
      S = m + (N//(2m+3) − 3)//2
      Δ(S) = S² + 6S − (N−9)   ?= k²
      Si sí: p = 2v+3, q = 2b+3  (en O(1), exacte)
==============================================================================
"""

from fractions import Fraction
import math
import time
from collections import deque


# ============================================================================
#  CONSTANTS
# ============================================================================
PRIMOS    = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)

LIM_KSWEEP    = 2_000_000   # valors de k en el k-sweep inicial
LIM_CRIBA     = 2_000_000   # criba de factors petits (<4M)

# Llindars d'adaptació de mostres (amplada del rang en unitats de m)
LLINDAR_GRAN  = 5_000_000   # rang gran → N_BASE mostres
LLINDAR_MITJ  = 50_000      # rang menut → k-sweep directe
N_MOSTRA_BASE = 30          # mostres per a rangs grans
N_MOSTRA_MIG  = 90          # mostres per a rangs mitjans

# Umbral de "zona calenta" (|d(m)| < UMBRAL → candidat fort)
UMBRAL_FASE   = Fraction(1, 5)   # = 0.20 (més permissiu que v1)
RADI_DISC     = 30               # radi de S en la cerca del discriminant


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL
# ============================================================================

def pasa_rueda(m: int) -> bool:
    """D(m)=2m+3 coprimer amb primorial p23? Densitat activa ~8.07%."""
    return m >= 1 and math.gcd(2 * m + 3, PRIMORIAL) == 1


# ============================================================================
#  BLOC 2 — FUNCIÓ D'ONA (Mètode B)
# ============================================================================

def d_frac(m: int, N: int) -> Fraction:
    """
    d(m) = Fraction(N mod D, D) − ½,  D = 2*(2m+3).
    Fraction exacta → zero errors d'arrodoniment.
    d(m) → 0 quan 2m+3 s'acosta a un factor de N.
    """
    D = 2 * (2 * m + 3)
    return Fraction(N % D, D) - Fraction(1, 2)


def check_factor(m: int, N: int) -> bool:
    """Verificació directa exacta: (2m+3) divideix N?"""
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


# ============================================================================
#  BLOC 3 — K-SWEEP (Mètode A)
# ============================================================================

def k_sweep(N: int, m_ini: int, m_fi: int) -> int | None:
    """
    Iteració sobre pendents k = N//(2m+3) en lloc de m.
    Exhaustiu i determinista per a [m_ini, m_fi].
    O(n_k) on n_k = k_hi − k_lo + 1.
    """
    if m_ini < 1: m_ini = 1
    if m_fi < m_ini: return None
    k_lo = max(1, N // (2 * m_fi  + 3))
    k_hi =        N // (2 * m_ini + 3)
    for k in range(k_lo, k_hi + 1):
        if k <= 0: continue
        c = N // k
        if c < 3 or c % 2 == 0: continue
        if N % c == 0 and 1 < c < N:
            return c
    return None


# ============================================================================
#  BLOC 4 — DISCRIMINANT S±k (Mètode C)
# ============================================================================

def discriminant_exacte(S: int, N: int) -> tuple:
    """
    Comprova si Δ(S) = S² + 6S − (N−9) és k² (quadrat perfecte).

    MATEMÀTICA (Libro 5, Teorema 1):
      N = (2v+3)(2b+3),  S = v+b,  k = v−b
      Δ(S) = k²  →  v=(S+k)//2, b=(S−k)//2
      p=2v+3, q=2b+3  →  p×q=N exactament

    RETORNA: (p, q) o (None, None).
    """
    if S < 0: return None, None
    delta = S * S + 6 * S - (N - 9)
    if delta < 0: return None, None
    k = math.isqrt(delta)
    if k * k != delta: return None, None
    Sp_k, Sm_k = S + k, S - k
    if Sp_k % 2 != 0 or Sm_k % 2 != 0: return None, None
    v, b = Sp_k // 2, Sm_k // 2
    if b < 0: return None, None
    p, q = 2*v+3, 2*b+3
    if p > 1 and q > 1 and p != N and q != N and p*q == N:
        return min(p, q), max(p, q)
    return None, None


def discriminant_des_de_m(m_zona: int, N: int) -> tuple:
    """
    Usa m com a pista per calcular S ≈ m + (N//(2m+3)−3)//2,
    aleshores explora S±RADI_DISC per trobar el quadrat perfecte.

    Connexió MDC → Discriminant: el MDC detecta un m calent,
    el discriminant el confirma algebraicament en O(RADI_DISC) passos.
    """
    p_cand = 2 * m_zona + 3
    if p_cand <= 1 or p_cand >= N: return None, None
    q_cand = N // p_cand
    if q_cand < 3: return None, None
    b = (q_cand - 3) // 2
    S_base = m_zona + b   # S estimat des de la zona calenta

    for ds in range(RADI_DISC + 1):
        for signe in ([0] if ds == 0 else [1, -1]):
            S = S_base + signe * ds
            if S < 0: continue
            p, q = discriminant_exacte(S, N)
            if p is not None:
                return p, q
    return None, None


# ============================================================================
#  BLOC 5 — DETECTOR DE ZONA CALENTA (Mètode B, amb mostres adaptatives)
# ============================================================================

def escanejar_rang_adaptatiu(N: int, m_ini: int, m_fi: int,
                              n_mostres: int) -> list[int]:
    """
    Escandeix [m_ini, m_fi] amb n_mostres punts equidistants.

    Per cada punt m vàlid (que passa la roda):
      · Verificació directa (check_factor) → factor immediat
      · Mesura d(m) exacta (Fraction)
      · Si |d(m)| < UMBRAL_FASE → zona calenta → afegir a llista
      · Salt cinemàtic: V=d(m+1)−d(m), A=d(m+2)−2d(m+1)+d(m)
        Si A≠0: m_proj = m − round(V/A) → afegir si cau en el rang

    RETORNA: llista de m candidats (zones calentes), ordenada per |d|.
    Retorna [−factor] si el factor és trobat directament (senyal especial).
    """
    amplada = m_fi - m_ini
    if amplada <= 0: return []
    pas = max(1, amplada // n_mostres)

    zonas = []   # (|d0|, m)

    m = m_ini
    while m <= m_fi:
        if pasa_rueda(m):

            # Verificació directa
            if check_factor(m, N):
                return [-(2 * m + 3)]   # senyal: factor trobat!

            # Mesura d(m) (Fraction exacta)
            D  = 2 * (2 * m + 3)
            R  = N % D
            d0 = Fraction(R, D) - Fraction(1, 2)

            abs_d0 = abs(d0)
            if abs_d0 < UMBRAL_FASE:
                zonas.append((abs_d0, m))

            # Salt cinemàtic (3 punts)
            if m + 2 <= m_fi:
                D1, D2 = 2*(2*(m+1)+3), 2*(2*(m+2)+3)
                d1 = Fraction(N % D1, D1) - Fraction(1, 2)
                d2 = Fraction(N % D2, D2) - Fraction(1, 2)
                V  = d1 - d0
                A  = d2 - 2*d1 + d0
                if A != Fraction(0):
                    despl = V / A
                    di = int(despl + Fraction(1,2)) if despl >= 0 \
                         else -int(-despl + Fraction(1,2))
                    m_proj = m - di
                    if m_ini <= m_proj <= m_fi and pasa_rueda(m_proj):
                        dp = abs(d_frac(m_proj, N))
                        if dp < UMBRAL_FASE * 2:
                            zonas.append((dp, m_proj))

        m += pas

    # Ordenar per |d| ascendent, eliminar duplicats
    vistos = set()
    resultat = []
    for _, mc in sorted(zonas, key=lambda x: x[0]):
        if mc not in vistos:
            vistos.add(mc)
            resultat.append(mc)

    return resultat


# ============================================================================
#  BLOC 6 — BISECCIÓ ADAPTATIVA (el cor de v2)
# ============================================================================

def biseccio_adaptativa(N: int, m_ini: int, m_fi: int,
                         verbose: bool = False) -> int | None:
    """
    Bisecció adaptativa del rang [m_ini, m_fi] amb densitat logarítmica.

    ALGORISME (cua de prioritat per amplada decreixent):
    ────────────────────────────────────────────────────
    Mantenim una cua de subrangs ordenats per prioritat.

    Per a cada subrang [a, b]:

      1. ADAPTACIÓ DE MOSTRES:
           amplada > LLINDAR_GRAN  → N_MOSTRA_BASE mostres (escàs)
           amplada ∈ [LLINDAR_MITJ, LLINDAR_GRAN] → N_MOSTRA_MIG mostres
           amplada < LLINDAR_MITJ  → k-sweep directe (exhaustiu)

      2. K-SWEEP DIRECTE si amplada petita:
           n_k del rang és manejable → exhaustiu, determinista.

      3. DETECTOR MDC si amplada gran:
           Escaneja amb n_mostres adaptatives.
           Per cada zona calenta: aplica discriminant S±k.

      4. BISECCIÓ si no trobat:
           Partim el rang en dos subrangs:
               [a, mig] i [mig, b]
           ORDRE DE PRIORITAT (distribució logarítmica):
               a. Subrang INFERIOR [a, mig] → primer (densitat > superior)
               b. Subrang SUPERIOR [mig, b] → segon

    PER QUÈ PRIMERO L'INFERIOR?
    ────────────────────────────
    La densitat de valors de k per unitat de m és:
        dk/dm = d/dm [N/(2m+3)] ≈ −2N/(2m+3)²
    Com més petit és m (factors menuts), la funció N/(2m+3) canvia
    més ràpidament → hi ha MÉS valors de k per unitat de m →
    MÉS oportunitats de trobar un factor en rangs baixos de m.

    Per tant, explorar primer el rang baix és estadísticament òptim.

    RETORNA: el factor (enter) o None.
    """
    # Cua de subrangs: (prioritat, m_ini, m_fi)
    # Prioritat MENOR = explorar PRIMER
    # Usem amplada invertida: rangs petits (zona baixa explorada) primer
    # Però la regla principal és: meitat inferior sempre és prioritat 0,
    # la meitat superior és prioritat 1 dins del mateix nivell.

    # Usem deque amb parells (nivell_biseccio, m_ini, m_fi)
    # Nivell 0 = rang original, nivell 1 = primers fills, etc.
    # Dins del mateix nivell: inferior avant, superior darrere.

    cua = deque()
    cua.appendleft((0, m_ini, m_fi))   # rang complet, prioritat màxima

    while cua:
        nivell, a, b = cua.popleft()
        amplada = b - a

        if a >= b or b < 1 or a > (math.isqrt(N)-3)//2:
            continue

        if verbose:
            print(f"    Bisec niv={nivell}  [{a:>12,}, {b:>12,}]"
                  f"  amp={amplada:>10,}", end="  ")

        # ── CAS 1: rang prou menut → k-sweep directe (exhaustiu) ──────────
        if amplada <= LLINDAR_MITJ:
            factor = k_sweep(N, a, b)
            if verbose:
                print(f"k-sweep n_k≈{min(amplada*2,999999):,}")
            if factor:
                return factor
            continue

        # ── CAS 2: rang gran → detector MDC adaptatiu ──────────────────────
        n_mostres = N_MOSTRA_BASE if amplada > LLINDAR_GRAN else N_MOSTRA_MIG

        zones = escanejar_rang_adaptatiu(N, a, b, n_mostres)

        if verbose:
            if zones and zones[0] < 0:
                print(f"MDC directe! Factor={-zones[0]:,}")
            else:
                print(f"MDC n={n_mostres} → {len(zones)} zones calentes")

        for m_zona in zones:

            # Senyal especial: factor trobat directament
            if m_zona < 0:
                f = -m_zona
                if N % f == 0 and f > 1 and f < N:
                    return f
                continue

            # Verificació directa al voltant de la zona
            for dm in range(-2, 3):
                m_t = m_zona + dm
                if m_t >= 1 and pasa_rueda(m_t) and check_factor(m_t, N):
                    return 2 * m_t + 3

            # Discriminant S±k des de la zona calenta
            p_d, q_d = discriminant_des_de_m(m_zona, N)
            if p_d is not None:
                return p_d

        # ── CAS 3: bisecció → encuem els dos subrangs ──────────────────────
        mig = (a + b) // 2

        # ORDRE LOGARÍTMIC: meitat inferior PRIMER (densitat major)
        cua.appendleft((nivell + 1, a,   mig))   # inferior → davant
        cua.append(    (nivell + 1, mig, b  ))   # superior → darrere

    return None


# ============================================================================
#  BLOC 7 — ORQUESTRADOR PRINCIPAL
# ============================================================================

def mdc_hybrid_v2(N: int, verbose: bool = True) -> tuple:
    """
    MDC PARABOLA HYBRID v2 — Bisecció Adaptativa.

    ┌─────────────────────────────────────────────────────────────────────┐
    │  F0: Precondicions (trivials + criba roda fins 2M)                 │
    │  F1: K-sweep [m_lim, m_max]  → factors equilibrats prop de √N     │
    │  F2: Bisecció adaptativa [1, m_lim]  → distribució logarítmica    │
    │       · Rangs grans: detector MDC (Fraction) + discriminant S±k   │
    │       · Rangs petits: k-sweep directe (exhaustiu)                  │
    │       · Prioritat: meitat inferior → meitat superior (per nivell)  │
    └─────────────────────────────────────────────────────────────────────┘
    """
    def log(msg):
        if verbose: print(msg)

    def retornar(p, q, etapa):
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose:
            bar = '═' * 68
            print(f"\n{bar}")
            print(f"  🎯 FACTOR TROBAT! [{etapa}]")
            print(f"     p = {p:,}")
            print(f"     q = {q:,}")
            print(f"     t = {t_ms:.3f} ms")
            print(f"{bar}")
        return (min(p, q), max(p, q)), t_ms

    if verbose:
        bar = '═' * 68
        print(f"\n{bar}")
        print(f"  MDC HYBRID v2 — Bisecció Adaptativa + Distribució 1/k")
        print(f"  N = {N:,}  ({len(str(N))} dígits  |  {N.bit_length()} bits)")
        print(f"{bar}")

    t0 = time.perf_counter()

    # ── F0: PRECONDICIONS ──────────────────────────────────────────────────
    log("\n  [F0] Precondicions + criba roda...")

    for p in PRIMOS:
        if N % p == 0 and p < N:
            log(f"    → Factor primorial: {p}")
            return retornar(p, N // p, "F0 trivial")

    r = math.isqrt(N)
    if r * r == N:
        return retornar(r, r, "F0 quadrat perfecte")

    m_max = (math.isqrt(N) - 3) // 2

    # Criba roda per a factors < 4M (cas asimètric: p << q)
    lim_c = min(m_max, LIM_CRIBA)
    log(f"    → Criba roda fins D={2*lim_c+3:,}  (m_max={m_max:,})")
    for m_c in range(1, lim_c + 1):
        if pasa_rueda(m_c):
            D = 2 * m_c + 3
            if N % D == 0 and D < N:
                log(f"    → Factor criba: {D:,}  (m={m_c:,})")
                return retornar(D, N // D, "F0 criba")

    if lim_c >= m_max:
        t_ms = (time.perf_counter() - t0) * 1000
        log(f"    → Criba cobreix tot el rang. N és primer.")
        return (None, None), t_ms

    t_criba = (time.perf_counter() - t0) * 1000
    log(f"    → {t_criba:.1f} ms. Factor no trobat per criba.")

    # ── F1: K-SWEEP PROP DE √N (factors molt equilibrats) ─────────────────
    log(f"\n  [F1] K-sweep {LIM_KSWEEP:,} k prop de √N...")

    k_top = max(1, N // (2 * m_max + 3))
    k_lim = k_top + LIM_KSWEEP
    d_lim = N // k_lim if k_lim > 0 else 2 * m_max + 3
    m_lim = max(lim_c + 1, (d_lim - 3) // 2)

    log(f"    → Rang k-sweep: m ∈ [{m_lim:,}, {m_max:,}]")

    factor_f1 = k_sweep(N, m_lim, m_max)
    if factor_f1:
        return retornar(factor_f1, N // factor_f1, "F1 k-sweep")

    t_f1 = (time.perf_counter() - t0) * 1000
    log(f"    → {t_f1:.1f} ms. Res trobat prop de √N.")

    # ── F2: BISECCIÓ ADAPTATIVA [1, m_lim] ────────────────────────────────
    log(f"\n  [F2] Bisecció adaptativa [1, {m_lim:,}]")
    log(f"    → Llindars: gran>{LLINDAR_GRAN:,} | menut<{LLINDAR_MITJ:,}")
    log(f"    → Mostres: gran={N_MOSTRA_BASE} | mig={N_MOSTRA_MIG} | menut=k-sweep")
    log(f"    → Prioritat: meitat inferior PRIMER (densitat logarítmica)")

    factor_f2 = biseccio_adaptativa(N, 1, m_lim, verbose=verbose)

    if factor_f2:
        return retornar(factor_f2, N // factor_f2, "F2 bisecció adaptativa")

    t_ms = (time.perf_counter() - t0) * 1000
    log(f"\n  ✗ No trobat. Temps total: {t_ms:.2f} ms")
    log(f"    (N podria ser primer, o cal augmentar LIM_KSWEEP/LLINDAR_MITJ)")
    return (None, None), t_ms


# ============================================================================
#  BLOC 8 — BATERIA DE PROVES
# ============================================================================

def executar_bateria(verbose_ind: bool = False):
    proves = [
        (3 * 5,                                      "trivial 3×5"),
        (101 * 103,                                  "clàssic 101×103"),
        (100003 * 100019,                            "target 100003×100019 (11d)"),
        (1_000_003 * 1_000_033,                      "equilibrat 13d"),
        (9_999_991 * 9_999_973,                      "equilibrat 14d"),
        (1_548_586_332_452_843,                      "factor 59 (16d, asimètric)"),
        (999_983 * 1_000_000_000_000_003,            "asimètric 6d×16d"),
        (100_000_000_003 * 100_000_000_019,          "equilibrat 22d"),
        (10_000_000_000_037 * 10_000_000_000_099,    "equilibrat 26d"),
    ]

    print("\n" + "█" * 68)
    print("  🚀 BATERIA — MDC HYBRID v2 (Bisecció Adaptativa)")
    print("  Distribució logarítmica 1/k · MDC + Discriminant S±k")
    print("█" * 68)

    ok = fail = 0
    for N, desc in proves:
        print(f"\n{'─'*68}")
        print(f"  📋 {desc}  ({len(str(N))} dígits)")
        (p, q), t_ms = mdc_hybrid_v2(N, verbose=verbose_ind)
        if p is not None:
            assert p * q == N, f"ERROR! {p}×{q}≠{N}"
            print(f"  ✅  {p:,} × {q:,}   ({t_ms:.2f} ms)")
            ok += 1
        else:
            print(f"  ⚠️   No trobat  ({t_ms:.2f} ms)")
            fail += 1

    print(f"\n{'█'*68}")
    print(f"  Encerts: {ok}/{ok+fail}")
    print("█" * 68)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

if __name__ == "__main__":

    # Prova verbose individual (el cas del prompt original)
    N = 100003 * 100019
    print(f"Prova verbose: N = {N:,} = 100.003 × 100.019")
    (p, q), t_ms = mdc_hybrid_v2(N, verbose=True)
    if p:
        print(f"\nResultat: {p:,} × {q:,}  ({t_ms:.3f} ms)")

    print()
    executar_bateria(verbose_ind=False)
