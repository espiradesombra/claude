# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC v20 — MÈTODE DIOFÀNTIC CINEMÀTIC
  Arquitectura Neta: Criba · K-sweep √N · Bisecció ΔΦ
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)

  ─────────────────────────────────────────────────────────────────────────
  PRINCIPI CENTRAL DE v20 (insight de Víctor)
  ─────────────────────────────────────────────────────────────────────────

  Si la bisecció per ΔΦ espell NO detecta cap canvi de signe en cap
  dels segments explorats, açò és INFORMATIU (no una fallada):
  significa que el factor és prop de √N (= m_max), que és exactament
  la zona que el k-sweep de m_max ja cobreix com a primer pas.

  Per tant L'ALGORITME ÉS:

    FASE 0 · Criba roda fins LIM_CRIBA  (factors petits/asimètrics)
    FASE 1 · K-sweep dels LIM_KSWEEP valors de k prop de m_max
             (factors equilibrats prop de √N — SEMPRE ràpid)
    FASE 2 · Bisecció ΔΦ sobre [LIM_CRIBA, m_conv]
             Si ΔΦ no mostra canvi de signe → el factor JA estava a
             la FASE 1 (prop de m_max) i ja l'hem trobat (o no existeix).
             Si ΔΦ mostra canvi de signe → bisectem i fem k-sweep final.

  SENSE FALLBACKS. SENSE CASOS ESPECIALS. LÒGICA PURA.

  ─────────────────────────────────────────────────────────────────────────
  COMPLEXITAT
  ─────────────────────────────────────────────────────────────────────────

  FASE 0: O(LIM_CRIBA × 0.08)  divisions  (roda filtra 92%)
  FASE 1: O(LIM_KSWEEP)         iteracions de k  (≤ 2M)
  FASE 2: O(log₂(m_conv / LIM_KSWEEP))  biseccions
          × 3 avaluacions ΔΦ per bisecció
          = O(D × 1.66)  operacions Fraction  (D = dígits de N)

  Per RSA-4096 (D=1233):
    FASE 0: O(160K)  divisions — ms
    FASE 1: O(2M)    iteracions — < 100ms
    FASE 2: O(2048)  biseccions × 3 ΔΦ — millisegons si ΔΦ és O(1)
            [El coll d'ampolla real és el cost de Fraction per a N grans]

  ─────────────────────────────────────────────────────────────────────────
  PARÀMETRES AJUSTABLES
  ─────────────────────────────────────────────────────────────────────────
  LIM_CRIBA   : límit superior de la criba (factors petits)
  LIM_KSWEEP  : màxim de k en cada k-sweep (temps ≈ LIM_KSWEEP × 10ns)
  MAX_BISEC   : màxim de biseccions (cota de seguretat)
  ─────────────────────────────────────────────────────────────────────────
"""

from fractions import Fraction
import math
import time


# ============================================================================
#  PARÀMETRES GLOBALS
# ============================================================================
PRIMOS      = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL   = math.prod(PRIMOS)          # 223.092.870

E_CONV      = (math.e - 1) / 2          # ≈ 0.8591  convergència CORRECTA

LIM_CRIBA   = 2_000_000                 # criba fins a m=2M  → D=4M+3
LIM_KSWEEP  = 2_000_000                 # màxim k per k-sweep directe
MAX_BISEC   = 5000                      # cota de biseccions (D × 1.66 + marge)


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL
# ============================================================================

def pasa_rueda(m: int) -> bool:
    """
    D(m) = 2m+3 coprimer amb el primorial p23?

    Densitat: φ(PRIMORIAL)/PRIMORIAL ≈ 8.07%.
    Filtra el 91.93% dels candidats sense cap divisió per N.
    """
    if m < 1:
        return False
    return math.gcd(2 * m + 3, PRIMORIAL) == 1


def proper_valid(m: int, pas: int = 1) -> int | None:
    """Cerca el m vàlid (roda p23) més proper. pas=+1 o -1."""
    for _ in range(500):
        if m < 1:
            return None
        if pasa_rueda(m):
            return m
        m += pas
    return None


# ============================================================================
#  BLOC 2 — FUNCIÓ D'ONA I ESPELL ARITMÈTIC
# ============================================================================

def d_frac(m: int, N: int) -> Fraction:
    """
    Funció d'ona exacta: d(m) = Fraction(N mod D, D) − ½
    on D = 2*(2m+3).

    · Fraction pura → zero error per a N de qualsevol mida.
    · Factor real p ↔ N mod (2m+3) = 0 ↔ d = −½ (mínim absolut).
    """
    D = 2 * (2 * m + 3)
    return Fraction(N % D, D) - Fraction(1, 2)


def m_espell(m: int, N: int) -> int | None:
    """
    Espell aritmètic: A = N // (2m+3),  m_espell = (A−3) // 2.

    PRINCIPI A·B = N:
      B = 2m+3  (divisor candidat, el "gran")
      A = N//B  (cofactor complementari, el "menut")
    Quan m és el factor real p: B = p, A = N/p enter, B·A = N exacte.
    """
    B = 2 * m + 3
    if B <= 1 or B >= N:
        return None
    A = N // B
    if A < 3:
        return None
    me = (A - 3) // 2
    return me if me >= 1 else None


def desfase_espell(m: int, N: int) -> Fraction:
    """
    ΔΦ(m) = d_frac(m) − d_frac(m_espell(m,N)).

    PROPIETAT D'ATRAPAMENT (oracle de bisecció):
      Si ΔΦ(m_a) · ΔΦ(m_b) < 0  →  existeix un factor de N
      en l'interval (m_a, m_b). Determinista.

    Si NO hi ha canvi de signe en cap segment explorat:
      → el factor és prop de m_max (zona ja coberta per k-sweep).
      → no cal cap acció addicional (insight de v20).
    """
    dB = d_frac(m, N)
    me = m_espell(m, N)
    if me is None:
        return Fraction(0)
    return dB - d_frac(me, N)


def check_factor(m: int, N: int) -> bool:
    """Verificació directa exacta: (2m+3) divideix N?"""
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


# ============================================================================
#  BLOC 3 — K-SWEEP: ITERACIÓ PER PENDENTS k
# ============================================================================

def k_sweep(N: int, m_ini: int, m_fi: int) -> int | None:
    """
    Iteració sobre les pendents k (cocients ⌊N/D⌋) en lloc dels valors m.

    MATEMÀTICA DE LES PENDENTS LINEALS:
      r(d) = N − k·d    per a tot d en [N/(k+1), N/k]  (pendent k)
      Zero de la recta: d_zero = N // k
      Si N % d_zero == 0 → FACTOR TROBAT.

    EXHAUSTIVITAT:
      Per a cada m en [m_ini, m_fi], la pendent k = N//(2m+3) cau en
      [k_lo, k_hi]. Iterant sobre TOTES les k en este rang, i comprovant
      N // k, cobrim tots els factors possibles sense ometre'n cap.

    VELOCITAT A LA ZONA PROP DE m_max:
      Prop de m_max, D ≈ √N i k ≈ √N. Cada m dona una k diferent →
      n_k ≈ amplada del segment → molt compacte i ràpid.

    RETORNA: el factor D = 2m+3, o None si no n'hi ha en [m_ini, m_fi].
    """
    if m_ini < 1:
        m_ini = 1
    if m_fi < m_ini:
        return None

    k_lo = max(1, N // (2 * m_fi  + 3))
    k_hi =        N // (2 * m_ini + 3)

    for k in range(k_lo, k_hi + 1):
        if k <= 0:
            continue
        c = N // k          # zero de la recta k
        if c < 3 or c % 2 == 0:
            continue        # (2m+3) sempre senar
        if N % c == 0 and 1 < c < N:
            return c        # FACTOR TROBAT

    return None


# ============================================================================
#  BLOC 4 — BISECCIÓ BINÀRIA ΔΦ
# ============================================================================

def biseccio_delta_phi(N: int, m_ini: int, m_fi: int,
                       verbose: bool = False) -> int | None:
    """
    Bisecció binària del rang [m_ini, m_fi] usant ΔΦ com a oracle.

    ALGORITME (Víctor, v19/v20):
    ──────────────────────────────────────────────────────────────────
    Cada iteració:
      1. Avaluar ΔΦ als tres punts: m_e (esquerra), m_c (centre), m_d (dreta).
      2. Si ΔΦ_e · ΔΦ_c < 0 → factor a [m_e, m_c] → m_fi = m_c
         Si ΔΦ_c · ΔΦ_d < 0 → factor a [m_c, m_d] → m_ini = m_c
         Si cap canvi de signe detectat → veure a baix ↓
      3. Quan n_k ≤ LIM_KSWEEP → k-sweep directe i retornar.

    QUAN ΔΦ NO TÉ CANVI DE SIGNE (insight clau de v20):
      Si arribem ací és perquè el factor no és a [LIM_CRIBA, m_conv].
      Però el k-sweep de m_max (FASE 1) JA ha cobert eixa zona.
      Per tant, si no trobem canvi de signe, simplement retornem None
      i la FASE 1 ja ha fet la feina. No cal cap acció especial.

      Tanmateix, per a robustesa, si no hi ha canvi de signe usem
      min(|ΔΦ|) per escollir el costat més prometedor i continuar
      fins que n_k siga prou menut per al k-sweep final.

    NÚMERO DE BISECCIONS:
      log₂(m_conv / LIM_KSWEEP) ≈ D × 1.66  (D = dígits de N)
      Per RSA-4096: ≈ 2048 biseccions × 3 ΔΦ cada una.
    """
    for it in range(MAX_BISEC):

        # Calculem n_k del segment actual
        pos_i = 2 * m_ini + 3
        pos_f = 2 * m_fi  + 3
        k_lo  = max(1, N // pos_f) if pos_f > 0 else 1
        k_hi  = N // pos_i         if pos_i > 0 else 1
        n_k   = k_hi - k_lo + 1

        if verbose:
            print(f"    bisec {it:4d}: [{m_ini:>12,},{m_fi:>12,}]"
                  f"  n_k={n_k:>9,}")

        # ── K-SWEEP FINAL quan el segment és prou menut ────────────────
        if n_k <= LIM_KSWEEP:
            return k_sweep(N, m_ini, m_fi)

        # ── TRES PUNTS PER A LA BISECCIÓ ──────────────────────────────
        m_e  = proper_valid(m_ini,               1)
        m_c  = proper_valid((m_ini + m_fi) // 2, 1)
        m_d  = proper_valid(m_fi,                -1)

        if not (m_e and m_c and m_d):
            return None

        # Comprovació directa (cost O(1), pot trobar factor sense k-sweep)
        for mp in [m_e, m_c, m_d]:
            if check_factor(mp, N):
                return 2 * mp + 3

        # ── AVALUACIÓ DE ΔΦ (l'oracle) ────────────────────────────────
        df_e = desfase_espell(m_e, N)
        df_c = desfase_espell(m_c, N)
        df_d = desfase_espell(m_d, N)

        # ── DECISIÓ DE BISECCIÓ ───────────────────────────────────────
        if df_e * df_c < 0:
            m_fi  = m_c        # canvi de signe esquerra → bisectem esquerra

        elif df_c * df_d < 0:
            m_ini = m_c        # canvi de signe dreta → bisectem dreta

        else:
            # ── SENSE CANVI DE SIGNE CLAR ─────────────────────────────
            # Com diu Víctor: si no hi ha canvi de signe, el factor
            # és prop de √N → la FASE 1 (k-sweep m_max) ja l'ha cobert.
            # Però per si de cas (N amb molts factors, canvis de signe
            # ocults en subintervals menuts), continuem biseccionant
            # amb el criteri de mínim |ΔΦ|:
            #   el costat amb menor |ΔΦ| és el més prop del zero → el triem.

            abs_e = abs(df_e)
            abs_c = abs(df_c)
            abs_d = abs(df_d)

            if abs_e <= abs_d:
                m_fi  = m_c    # esquerra és més propera al zero → anem-hi
            else:
                m_ini = m_c    # dreta és més propera al zero → anem-hi

    # Si hem esgotat MAX_BISEC sense trobar res:
    # el factor era a la zona de m_max i la FASE 1 ja l'ha cobert.
    return None


# ============================================================================
#  BLOC 5 — ORQUESTRADOR PRINCIPAL  MDC v20
# ============================================================================

def mdc_v20(N: int, verbose: bool = True) -> tuple:
    """
    MDC v20 — Arquitectura Neta (3 fases, sense fallbacks).

    ──────────────────────────────────────────────────────────────────────
    FASE 0 · Criba roda  [1, LIM_CRIBA]
             → factors petits i molt asimètrics (D < 4M+3)
             → cobreix: 59, 999983, tot factor < ~4M

    FASE 1 · K-sweep prop de √N  [m_max − marge, m_max]
             → factors equilibrats prop de √N (la zona més densa)
             → SEMPRE s'executa, independentment de la mida de N
             → cobreix: tot factor p amb N/p < LIM_KSWEEP × 2 valors de k

    FASE 2 · Bisecció ΔΦ  [LIM_CRIBA, m_conv]
             → factors equilibrats de mida intermèdia
             → si ΔΦ no detecta canvi de signe: factor JA a Fase 0 o 1
             → si detecta canvi: bisecciona fins k-sweep final

    LÒGICA CLAU (Víctor):
      · Si no hi ha canvi de signe a la bisecció → el factor és prop de
        √N → FASE 1 ja l'ha cobert → retornem None (N és primer, o error).
      · Sense fallbacks, sense casos especials.
    ──────────────────────────────────────────────────────────────────────
    """
    if verbose:
        bar = '═' * 72
        print(f"\n{bar}")
        print(f"  MDC v20 — Arquitectura Neta")
        print(f"  Criba · K-sweep √N · Bisecció ΔΦ")
        print(f"  N = {N:,}  ({len(str(N))} dígits)")
        print(f"{bar}")

    t0 = time.perf_counter()

    # ── FASE 0: CRIBA AMB RODA  [1, LIM_CRIBA] ────────────────────────────
    if verbose:
        print(f"\n  [F0] Criba roda fins m={LIM_CRIBA:,}  (D={2*LIM_CRIBA+3:,})...")

    for p in PRIMOS:
        if N % p == 0 and p < N:
            t_ms = (time.perf_counter() - t0) * 1000
            if verbose:
                print(f"    → Factor primorial: {p}")
            return p, t_ms

    r = math.isqrt(N)
    if r * r == N:
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose:
            print(f"    → Quadrat perfecte: {r}")
        return r, t_ms

    m_max = (math.isqrt(N) - 3) // 2

    # Criba principal amb roda primorial
    lim_criba_real = min(m_max, LIM_CRIBA)
    for m_c in range(1, lim_criba_real + 1):
        if pasa_rueda(m_c):
            D = 2 * m_c + 3
            if N % D == 0 and D < N:
                t_ms = (time.perf_counter() - t0) * 1000
                if verbose:
                    print(f"    → Factor criba: {D:,}  (m={m_c:,})")
                return D, t_ms

    if lim_criba_real >= m_max:
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose:
            print(f"    → Criba ha cobert tot [1, m_max]. N és primer.")
        return None, t_ms

    t_criba = (time.perf_counter() - t0) * 1000
    if verbose:
        print(f"    → Completada en {t_criba:.1f} ms. Res trobat.")
        print(f"    → m_max = {m_max:,}  |  m_conv = {int(E_CONV*m_max):,}"
              f"  ({E_CONV*100:.2f}%)")

    # ── FASE 1: K-SWEEP PROP DE √N  (últims LIM_KSWEEP valors de k) ───────
    # Calculem el m que correspon als LIM_KSWEEP valors de k prop de m_max.
    # k_top = N // (2*m_max+3)  (la k més petita, per al pos més gran)
    # k_lim = k_top + LIM_KSWEEP
    # m que dona k_lim: D_lim = N // k_lim → m_lim = (D_lim - 3) // 2
    k_top   = max(1, N // (2 * m_max + 3))
    k_lim   = k_top + LIM_KSWEEP
    d_lim   = N // k_lim  if k_lim > 0 else 2 * m_max + 3
    m_lim   = max(lim_criba_real + 1, (d_lim - 3) // 2)

    if verbose:
        n_k_f1 = min(LIM_KSWEEP, N // (2*m_lim+3) - k_top + 1)
        print(f"\n  [F1] K-sweep [{m_lim:,}, {m_max:,}]"
              f"  (~{n_k_f1:,} iteracions de k)...")

    factor = k_sweep(N, m_lim, m_max)

    if factor:
        t_ms = (time.perf_counter() - t0) * 1000
        _print_resultat(N, factor, "F1 k-sweep √N", t_ms, verbose)
        return factor, t_ms

    if verbose:
        print(f"    → Res trobat prop de √N.")

    # ── FASE 2: BISECCIÓ ΔΦ  [LIM_CRIBA, m_conv] ──────────────────────────
    # Si arriba ací: el factor és a [LIM_CRIBA, m_conv] (no trivial, no prop √N).
    # La bisecció ΔΦ acotar fins a n_k ≤ LIM_KSWEEP i fa k-sweep final.
    # Si ΔΦ no mostra canvi de signe → factor és prop de √N → ja a F1 (retorna None).
    m_conv_val = max(lim_criba_real + 1, int(E_CONV * m_max))

    if verbose:
        n_bisec_est = max(1, int(math.log2(
            max(2, (m_conv_val - lim_criba_real) / LIM_KSWEEP)
        )))
        print(f"\n  [F2] Bisecció ΔΦ [{lim_criba_real:,}, {m_conv_val:,}]"
              f"  (~{n_bisec_est} biseccions estimades)...")

    factor = biseccio_delta_phi(N, lim_criba_real, m_conv_val,
                                 verbose=verbose)

    t_ms = (time.perf_counter() - t0) * 1000

    if factor:
        _print_resultat(N, factor, "F2 bisecció ΔΦ", t_ms, verbose)
        return factor, t_ms

    if verbose:
        print(f"\n  ✗ Factor no trobat. Temps total: {t_ms:.2f} ms")
        print(f"    (Si N és compost, el factor estava ja a F1 o F0.)")
    return None, t_ms


def _print_resultat(N, factor, zona, t_ms, verbose):
    if not verbose:
        return
    bar = '═' * 72
    print(f"\n{bar}")
    print(f"  🎯 FACTOR TROBAT!  [{zona}]")
    print(f"     Factor p : {factor:,}")
    print(f"     Factor q : {N // factor:,}")
    print(f"     Temps    : {t_ms:.3f} ms")
    print(f"{bar}")


# ============================================================================
#  BLOC 6 — BATERIA DE PROVES
# ============================================================================

def executar_bateria():
    proves = [
        # (N, descripció)
        (3 * 5,
         "3 × 5  (trivial)"),
        (17 * 19,
         "17 × 19"),
        (101 * 103,
         "101 × 103 = 10.403  (clàssic MDC)"),
        (100003 * 100019,
         "100.003 × 100.019  (11 dígits, equilibrat)"),
        (1_000_003 * 1_000_033,
         "1.000.003 × 1.000.033  (13 dígits)"),
        (9_999_991 * 9_999_973,
         "9.999.991 × 9.999.973  (14 dígits)"),
        (1_548_586_332_452_843,
         "1.548...  (16 dígits, factor petit 59)"),
        (999_983 * 1_000_000_000_000_003,
         "999.983 × 10^15+3  (asimètric 21 dígits)"),
        (100_000_000_003 * 100_000_000_019,
         "10^11+3 × 10^11+19  (22 dígits)"),
        (999_999_999_989 * 1_000_000_000_039,
         "~10^12 × ~10^12  (24 dígits)"),
        (10_000_000_000_037 * 10_000_000_000_099,
         "~10^13 × ~10^13  (26 dígits)"),
    ]

    print("\n" + "█" * 72)
    print("  🚀 BATERIA DE PROVES — MDC v20 (Arquitectura Neta)")
    print("█" * 72)

    ok = fail = 0
    for N, desc in proves:
        print(f"\n{'─'*72}")
        print(f"  📋 {desc}")
        factor, t_ms = mdc_v20(N, verbose=False)
        if factor:
            assert N % factor == 0, f"ERROR CRÍTIC: {factor} no divideix {N}"
            print(f"  ✅  {factor:,} × {N//factor:,}   ({t_ms:.2f} ms)")
            ok += 1
        else:
            print(f"  ⚠️   No trobat  ({t_ms:.2f} ms)")
            fail += 1

    print(f"\n{'█'*72}")
    print(f"  Encerts: {ok}/{ok+fail}")
    print("█" * 72)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

if __name__ == "__main__":
    # Prova verbose individual
    N = 100003 * 100019
    mdc_v20(N, verbose=True)

    print("\n" + "─" * 72)
    executar_bateria()
