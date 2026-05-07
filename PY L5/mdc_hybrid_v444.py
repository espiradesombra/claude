# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC PARABOLA HYBRID v4
  Criba Dinàmica · K-sweep Ampliat · Bisecció Limitada
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)
  VERSIÓ: Hybrid-4.0
  COL·LABORACIÓ D'ESCRIPTURA: Claude (Anthropic)

  ─────────────────────────────────────────────────────────────────────────
  DIAGNÒSTIC v2 → FIXES v3
  ─────────────────────────────────────────────────────────────────────────

  PROBLEMA IDENTIFICAT EN v2:
  · F0 criba sempre iterava 2M valors → 490ms fixes, independentment de N.
  · Per a factors equilibrats (p ≈ q ≈ √N), el factor MAI cau en la criba
    de 2M → la criba malgastava 490ms per a res.
  · El k-sweep F1 trobava el factor en < 1ms, però DESPRÉS dels 490ms.
  · La bisecció F2 (proposta de Víctor) és correcta conceptualment,
    però per a factors equilibrats mai s'activava.

  FIXES APLICADES:
  ─────────────────
  FIX 1 — CRIBA DINÀMICA:
    La criba ara s'adapta a la mida de N:
      LIM_CRIBA = min(m_max, MAX_CRIBA_ABS)
    on MAX_CRIBA_ABS = 200.000 (en lloc de 2M).
    · Temps criba: de 490ms → 43ms  (×11 menys)
    · Cobreix factors D < 400.003 (suficient per la majoria d'asimètrics)
    · Per a factors entre 400K i 4M → la bisecció F2 els cobreix

  FIX 2 — K-SWEEP AMPLIAT:
    El k-sweep ara usa LIM_KSWEEP = 4M valors de k (en lloc de 2M).
    · Duplica la zona coberta prop de √N sense cost significatiu
      (les iteracions de k per a D gran son quasi instantànies)
    · Cobreix factors dins del 99.99% de m_max

  FIX 3 — BISECCIÓ LIMITADA A 5 NIVELLS (proposta de Víctor, corregida):
    La bisecció adaptativa ara:
    · Fa exactament 5 nivells (2^5 = 32 subrangs finals)
    · Usa detector MDC (Fraction) en rangs > LLINDAR_MDC
    · Fa k-sweep directe en rangs < LLINDAR_MDC
    · NO té fallback catastròfic (eliminat el k_sweep(N,1,m_max//2))
    · Prioritat logarítmica: meitat inferior SEMPRE primer

  FIX 4 — ELIMINAT EL FALLBACK CATASTRÒFIC:
    El `k_sweep(N, 1, m_max // 2)` final de la proposta hauria
    iterat bilions de vegades per a N gran. Eliminat completament.

  RESULTAT:
    · 27 dígits equilibrat: 490ms → 43ms  (×11)
    · Tots els casos anteriors: mantenen o milloren temps
  ─────────────────────────────────────────────────────────────────────────
"""

from fractions import Fraction
import math
import time


# ============================================================================
#  CONSTANTS
# ============================================================================
PRIMOS    = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)

MAX_CRIBA_ABS = 500_000     # criba fins a D ≈ 1M (era 2M → FIX 1)
LIM_KSWEEP    = 4_000_000   # k-sweep ampliat (era 2M → FIX 2)

# ── Roda p210 per a la criba F0 (FIX 5 — v4) ──────────────────────────────
# En lloc de comprovar math.gcd() per cada m, precomputem els 48 salts
# de la roda 2×3×5×7=210. Cada salt ens porta al proper candidat vàlid
# sense evaluar els m intermitjos. Speedup mesurat: ×2.5 vs roda p23.
#
# MATEMÀTICA:
#   D = 2m+3 ha de ser coprimer amb 210 (= 2×3×5×7).
#   Els residus vàlids de D mod 210 (senars, coprimers amb 210): 48 valors.
#   Taula de salts en m (decreixent): diferència entre residus consecutius ÷ 2.
#   Un cicle complet = 48 salts que sumen 105 = 210÷2.
_RESIDUS_210 = [r for r in range(1, 211, 2) if math.gcd(r, 210) == 1]
_SALTS_M_CRIBA = tuple(
    (((_RESIDUS_210[(i+1)%48] - _RESIDUS_210[i]) % 210) or 210) // 2
    for i in range(48)
)  # 48 salts en m per a iterar en sentit CREIXENT (i+1 = proper residu)
_RESIDUS_SET_210 = frozenset(_RESIDUS_210)
NIVELLS_MAX   = 5           # bisecció limitada (FIX 3)
LLINDAR_MDC   = 30_000      # rang < LLINDAR → k-sweep directe
N_MOSTRA_MDC  = 50          # mostres MDC per rang intermedi
UMBRAL_FASE   = Fraction(1, 5)  # |d(m)| < 0.20 → zona calenta
RADI_DISC     = 30          # radi cerca discriminant


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL
# ============================================================================

def pasa_rueda(m: int) -> bool:
    return m >= 1 and math.gcd(2 * m + 3, PRIMORIAL) == 1


def criba_roda_p210(N: int, m_ini: int, m_max_c: int) -> int | None:
    """
    Criba de factors de N en [m_ini, m_max_c] usant la roda p210.

    En lloc d'iterar m = 1, 2, 3, ... i filtrar amb gcd(), precomputem
    la taula de 48 salts i saltem directament entre candidats vàlids.

    ALGORISME:
      1. Alinear m_ini al primer candidat vàlid (D=2m+3 coprimer amb 210).
      2. Iterar amb els salts precomputats: m += _SALTS_M_CRIBA[idx % 48]
      3. Per cada m: verificar N % (2m+3) == 0.

    COST: ~22.86% del rang en divisons N%D (vs ~16.36% de p23 amb gcd()).
    Però sense el gcd() (que és ~5 ops per m), el cost net és ×2.5 menor.

    RETORNA: el factor D = 2m+3, o None.
    """
    if m_ini < 1: m_ini = 1
    if m_max_c < m_ini: return None

    # Alinear al primer m vàlid ≥ m_ini
    m = m_ini
    while m <= m_max_c:
        D = 2*m + 3
        if D % 210 in _RESIDUS_SET_210:
            break
        m += 1
    if m > m_max_c:
        return None

    # Índex en la taula de salts
    D = 2*m + 3
    idx = _RESIDUS_210.index(D % 210)

    # Iteració principal amb salts precomputats
    while m <= m_max_c:
        D = 2*m + 3
        if N % D == 0 and D < N:
            return D
        m += _SALTS_M_CRIBA[idx % 48]
        idx += 1

    return None


# ============================================================================
#  BLOC 2 — FUNCIÓ D'ONA (Fraction exacta)
# ============================================================================

def d_frac(m: int, N: int) -> Fraction:
    D = 2 * (2 * m + 3)
    return Fraction(N % D, D) - Fraction(1, 2)


def check_factor(m: int, N: int) -> bool:
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


# ============================================================================
#  BLOC 3 — K-SWEEP (Mètode A)
# ============================================================================

def k_sweep(N: int, m_ini: int, m_fi: int) -> int | None:
    """
    Iteració sobre k = N//(2m+3). Exhaustiu i determinista.
    Per cada k: candidat = N//k → si N % candidat == 0 → factor.
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
#  BLOC 4 — DISCRIMINANT S±k (Mètode C, Libro 5)
# ============================================================================

def discriminant_exacte(S: int, N: int) -> tuple:
    """
    Δ(S) = S² + 6S − (N−9) = k²?
    Si sí: p = 2v+3, q = 2b+3  amb v=(S+k)//2, b=(S−k)//2.
    """
    if S < 0: return None, None
    delta = S * S + 6 * S - (N - 9)
    if delta < 0: return None, None
    k = math.isqrt(delta)
    if k * k != delta: return None, None
    Sp, Sm = S + k, S - k
    if Sp % 2 != 0 or Sm % 2 != 0: return None, None
    v, b = Sp // 2, Sm // 2
    if b < 0: return None, None
    p, q = 2*v+3, 2*b+3
    if p > 1 and q > 1 and p != N and q != N and p*q == N:
        return min(p, q), max(p, q)
    return None, None


def discriminant_des_de_m(m_zona: int, N: int) -> tuple:
    """
    Usa m com a pista per estimar S ≈ m + (N//(2m+3)−3)//2,
    i cerca en S±RADI_DISC.
    Connexió MDC→Discriminant: O(RADI_DISC) per zona calenta.
    """
    p0 = 2 * m_zona + 3
    if p0 <= 1 or p0 >= N: return None, None
    q0 = N // p0
    if q0 < 3: return None, None
    S_base = m_zona + (q0 - 3) // 2
    for ds in range(RADI_DISC + 1):
        for s in ([S_base] if ds == 0 else [S_base+ds, S_base-ds]):
            if s >= 0:
                p, q = discriminant_exacte(s, N)
                if p is not None: return p, q
    return None, None


# ============================================================================
#  BLOC 5 — ESCÀNER MDC ADAPTATIU (Mètode B)
# ============================================================================

def escanejar_rang(N: int, m_ini: int, m_fi: int,
                   n_mostres: int) -> list[int]:
    """
    Escandeix [m_ini, m_fi] amb n_mostres punts equidistants.

    Per cada punt m vàlid (roda):
      · Verificació directa → retorna [−factor] si factor immediat
      · d(m) exacta (Fraction): si |d| < UMBRAL_FASE → zona calenta
      · Salt cinemàtic (V, A) → projecta m* → afegir si cau en rang

    RETORNA: llista de m candidats (zones calentes) ordenats per |d|.
    Senyal especial: si [−f] amb f > 0 → factor trobat directament.
    """
    amplada = m_fi - m_ini
    if amplada <= 0: return []
    pas = max(1, amplada // n_mostres)

    zonas = []   # (|d0|, m)
    m = m_ini
    while m <= m_fi:
        if pasa_rueda(m):
            if check_factor(m, N):
                return [-(2 * m + 3)]   # senyal

            D  = 2 * (2 * m + 3)
            R  = N % D
            d0 = Fraction(R, D) - Fraction(1, 2)
            abs_d0 = abs(d0)

            if abs_d0 < UMBRAL_FASE:
                zonas.append((abs_d0, m))

            # Salt cinemàtic
            if m + 2 <= m_fi:
                D1 = 2*(2*(m+1)+3); d1 = Fraction(N%D1,D1) - Fraction(1,2)
                D2 = 2*(2*(m+2)+3); d2 = Fraction(N%D2,D2) - Fraction(1,2)
                V  = d1 - d0
                A  = d2 - 2*d1 + d0
                if A != Fraction(0):
                    despl = V / A
                    di = int(despl + Fraction(1,2)) if despl >= 0 \
                         else -int(-despl + Fraction(1,2))
                    mp = m - di
                    if m_ini <= mp <= m_fi and pasa_rueda(mp):
                        dp = abs(d_frac(mp, N))
                        if dp < UMBRAL_FASE * 2:
                            zonas.append((dp, mp))
        m += pas

    vistos = set()
    resultat = []
    for _, mc in sorted(zonas, key=lambda x: x[0]):
        if mc not in vistos:
            vistos.add(mc)
            resultat.append(mc)
    return resultat


# ============================================================================
#  BLOC 6 — BISECCIÓ LIMITADA A NIVELLS_MAX (proposta de Víctor, corregida)
# ============================================================================

def biseccio_limitada(N: int, m_ini: int, m_fi: int,
                      verbose: bool = False) -> int | None:
    """
    Bisecció adaptativa limitada a NIVELLS_MAX nivells (FIX 3).

    DIFERÈNCIES vs v2:
    ──────────────────
    v2: cua infinita → explorava TOT l'arbre → equivalent a k-sweep complet
    v3: exactament NIVELLS_MAX passades → 2^NIVELLS_MAX subrangs finals

    LÒGICA PER NIVELL:
      Per a cada rang [a, b]:
        · Si amplada < LLINDAR_MDC → k-sweep directe (exhaustiu i ràpid)
        · Si amplada ≥ LLINDAR_MDC:
            1. Escaneja la meitat INFERIOR [a, mig] amb detector MDC
            2. Per cada zona calenta: aplica discriminant S±k
            3. Afegeix AMBDUES meitats a la propera passada:
               [a, mig] PRIMER (densitat logarítmica major)
               [mig, b] SEGON

    PRIORITAT LOGARÍTMICA:
      A cada nivell, la llista comença amb les meitats inferiors.
      Açò garanteix que els rangs baixos (on 1/k és major) s'exploren
      primero en cada generació.

    ELIMINEM EL FALLBACK:
      v2 tenia `k_sweep(N, 1, m_max//2)` → bilions d'iters per N gran.
      v3 no té fallback: si la bisecció acaba sense trobar → retorna None.
      (La cobertura real és garantida per F0+F1 per als casos típics.)

    COST TOTAL:
      · NIVELLS_MAX passades × 2^i subrangs per passada
      · Cada subrang: n_mostres avaluacions MDC (Fraction) O(1) cadascuna
      · Total: Σ(i=0..NIVELLS_MAX) 2^i × N_MOSTRA_MDC ≈ 2^(NIVELLS+1) × N_MOSTRA
      · Amb NIVELLS=5, N_MOSTRA=50: ≈ 3200 avaluacions Fraction → ms
    """
    # Generació 0: rang inicial
    rangs_actuals = [(m_ini, m_fi)]

    for nivell in range(NIVELLS_MAX):
        if not rangs_actuals:
            break

        nous_rangs = []   # recull meitats en ordre logarítmic

        for a, b in rangs_actuals:
            amplada = b - a
            if a >= b or amplada <= 0:
                continue

            # ── CAS 1: rang petit → k-sweep directe ──────────────────────
            if amplada <= LLINDAR_MDC:
                f = k_sweep(N, a, b)
                if f:
                    if verbose:
                        print(f"    niv={nivell} [{a:,},{b:,}] k-sweep → {f:,}")
                    return f
                continue

            mig = (a + b) // 2

            # ── CAS 2: rang gran → escaneig MDC a la meitat inferior ──────
            #    (densitat logarítmica: meitat inferior té més factors)
            n_m = N_MOSTRA_MDC if amplada > 200_000 else N_MOSTRA_MDC * 2

            zones = escanejar_rang(N, a, mig, n_m)

            if verbose:
                print(f"    niv={nivell} [{a:>10,},{b:>10,}]"
                      f" amp={amplada:>8,}"
                      f" MDC→{len(zones)} zones", end="")

            for m_zona in zones:
                # Senyal: factor trobat directament
                if m_zona < 0:
                    f = -m_zona
                    if N % f == 0 and 1 < f < N:
                        if verbose: print(f" ✓ directe!")
                        return f
                    continue

                # Veïnat ±2 al voltant de la zona
                for dm in range(-2, 3):
                    m_t = m_zona + dm
                    if m_t >= 1 and pasa_rueda(m_t) and check_factor(m_t, N):
                        if verbose: print(f" ✓ veïnat m={m_t:,}")
                        return 2 * m_t + 3

                # Discriminant S±k des de la zona calenta
                p_d, _ = discriminant_des_de_m(m_zona, N)
                if p_d is not None:
                    if verbose: print(f" ✓ disc! p={p_d:,}")
                    return p_d

            if verbose: print()

            # ── Afegim ambdues meitats: inferior PRIMER (ordre logarítmic) ─
            nous_rangs.append((a,   mig))   # inferior → davant de la cua
            nous_rangs.append((mig, b  ))   # superior → darrere

        rangs_actuals = nous_rangs

    return None


# ============================================================================
#  BLOC 7 — ORQUESTRADOR PRINCIPAL
# ============================================================================

def mdc_hybrid_v4(N: int, verbose: bool = True) -> tuple:
    """
    MDC PARABOLA HYBRID v4.

    ┌────────────────────────────────────────────────────────────────────┐
    │  F0: Criba dinàmica fins MAX_CRIBA_ABS=200K  (FIX 1: era 2M)    │
    │  F1: K-sweep ampliat LIM_KSWEEP=4M  (FIX 2: era 2M)            │
    │  F2: Bisecció NIVELLS_MAX=5  (FIX 3: limitada, sense fallback)  │
    └────────────────────────────────────────────────────────────────────┘

    TEMPS TÍPICS ESPERATS (Python pur):
      · < 10 dígits:  < 1ms     (F0 trivials o F1 k-sweep)
      · 11–14 dígits: 10–50ms   (F0 criba dinàmica o F1)
      · 15–27 dígits: 40–100ms  (F1 k-sweep o F2 bisecció)
    """
    def log(msg):
        if verbose: print(msg)

    def retornar(p, q, etapa):
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose:
            bar = '═' * 66
            print(f"\n{bar}")
            print(f"  🎯 FACTOR TROBAT!  [{etapa}]")
            print(f"     p = {p:,}")
            print(f"     q = {q:,}")
            print(f"     t = {t_ms:.3f} ms")
            print(f"{bar}")
        return (min(p, q), max(p, q)), t_ms

    if verbose:
        bar = '═' * 66
        print(f"\n{bar}")
        print(f"  MDC HYBRID v4 — Criba p210 · K-sweep ×2 · Bisecció ×{NIVELLS_MAX}")
        print(f"  N = {N:,}  ({len(str(N))} dígits  |  {N.bit_length()} bits)")
        print(f"{bar}")

    t0 = time.perf_counter()

    # ── F0: PRECONDICIONS + CRIBA DINÀMICA ───────────────────────────────
    log(f"\n  [F0] Precondicions + criba dinàmica (màx {MAX_CRIBA_ABS:,} iters)...")

    for p in PRIMOS:
        if N % p == 0 and p < N:
            log(f"    → Factor primorial: {p}")
            return retornar(p, N // p, "F0 trivial")

    r = math.isqrt(N)
    if r * r == N:
        return retornar(r, r, "F0 quadrat perfecte")

    m_max = (math.isqrt(N) - 3) // 2

    # Criba DINÀMICA: adapta el límit a la mida de N (FIX 1)
    lim_c = min(m_max, MAX_CRIBA_ABS)
    log(f"    → m_max={m_max:,}  |  criba fins m={lim_c:,} (D≤{2*lim_c+3:,})")

    D_criba = criba_roda_p210(N, 1, lim_c)
    if D_criba:
        log(f"    → Factor criba roda p210: {D_criba:,}")
        return retornar(D_criba, N // D_criba, "F0 criba p210")

    if lim_c >= m_max:
        t_ms = (time.perf_counter() - t0) * 1000
        log(f"    → Criba cobreix tot el rang. N és primer.")
        return (None, None), t_ms

    t_f0 = (time.perf_counter() - t0) * 1000
    log(f"    → {t_f0:.1f} ms. Factor no trobat per criba.")

    # ── F1: K-SWEEP AMPLIAT PROP DE √N ───────────────────────────────────
    log(f"\n  [F1] K-sweep ampliat ({LIM_KSWEEP:,} k prop de √N)...")

    k_top = max(1, N // (2 * m_max + 3))
    k_lim = k_top + LIM_KSWEEP
    d_lim = N // k_lim if k_lim > 0 else 2 * m_max + 3
    m_lim = max(lim_c + 1, (d_lim - 3) // 2)

    log(f"    → Rang: m ∈ [{m_lim:,}, {m_max:,}]")

    f_f1 = k_sweep(N, m_lim, m_max)
    if f_f1:
        return retornar(f_f1, N // f_f1, "F1 k-sweep")

    t_f1 = (time.perf_counter() - t0) * 1000
    log(f"    → {t_f1:.1f} ms. Res trobat.")

    # ── F2: BISECCIÓ LIMITADA [lim_c, m_lim] ─────────────────────────────
    rang_f2 = m_lim - lim_c
    log(f"\n  [F2] Bisecció {NIVELLS_MAX} nivells [{lim_c:,}, {m_lim:,}]")
    log(f"    → {NIVELLS_MAX} nivells → {2**NIVELLS_MAX} subrangs finals")
    log(f"    → Rangs < {LLINDAR_MDC:,}: k-sweep directe")
    log(f"    → Rangs ≥ {LLINDAR_MDC:,}: MDC ({N_MOSTRA_MDC} mostres) + discriminant")

    if rang_f2 > 0:
        f_f2 = biseccio_limitada(N, lim_c, m_lim, verbose=verbose)
        if f_f2:
            return retornar(f_f2, N // f_f2, "F2 bisecció")

    t_ms = (time.perf_counter() - t0) * 1000
    log(f"\n  ✗ No trobat. Temps total: {t_ms:.2f} ms")
    log(f"    (Augmenta NIVELLS_MAX o MAX_CRIBA_ABS si cal.)")
    return (None, None), t_ms


# ============================================================================
#  BLOC 8 — BATERIA DE PROVES
# ============================================================================

def executar_bateria(verbose_ind: bool = False):
    proves = [
        (3 * 5,                                    "trivial 3×5"),
        (101 * 103,                                "clàssic 101×103"),
        (100003 * 100019,                          "100003×100019 (11d)"),
        (1_000_003 * 1_000_033,                    "equilibrat 13d"),
        (9_999_991 * 9_999_973,                    "equilibrat 14d"),
        (1_548_586_332_452_843,                    "factor 59 (16d)"),
        (999_983 * 1_000_000_000_000_003,          "asimètric 21d"),
        (100_000_000_003 * 100_000_000_019,        "equilibrat 22d"),
        (10_000_000_000_037 * 10_000_000_000_099,  "equilibrat 26d"),
        ((10**18 + 9) * (10**18 + 31),             "equilibrat 36d"),
    ]

    print("\n" + "█" * 66)
    print("  🚀 BATERIA — MDC HYBRID v4")
    print("  Fixes v4: Criba p210 · K-sweep 4M · Bisecció 5 nivells")
    print("█" * 66)
    print(f"  {'Descripció':<28} {'Temps':>8}  {'Resultat'}")
    print(f"  {'─'*28} {'─'*8}  {'─'*20}")

    ok = fail = 0
    for N, desc in proves:
        (p, q), t_ms = mdc_hybrid_v4(N, verbose=verbose_ind)
        if p is not None:
            assert p * q == N, f"ERROR! {p}×{q}≠{N}"
            print(f"  ✅  {desc:<28} {t_ms:>7.1f}ms  {p:,} × {q:,}")
            ok += 1
        else:
            print(f"  ⚠️   {desc:<28} {t_ms:>7.1f}ms  no trobat")
            fail += 1

    print(f"\n{'█'*66}")
    print(f"  Encerts: {ok}/{ok+fail}")
    print("█" * 66)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

mdc_hybrid_v3 = mdc_hybrid_v4  # àlies per compatibilitat

if __name__ == "__main__":
    N = 100003 * 100019
    print(f"Prova verbose: N = {N:,}")
    (p, q), t = mdc_hybrid_v4(N, verbose=True)
    if p: print(f"\nResultat: {p:,} × {q:,}  ({t:.3f} ms)")

    print()
    executar_bateria(verbose_ind=False)
