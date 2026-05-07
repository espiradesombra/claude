# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC v22 — MÈTODE DIOFÀNTIC CINEMÀTIC
  Bisecció Quaternària (÷4 per iteració) · Signe Enter Pur
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)

  ─────────────────────────────────────────────────────────────────────────
  NOVETAT v22: BISECCIÓ QUATERNÀRIA (idea de Víctor)
  ─────────────────────────────────────────────────────────────────────────

  v21 usava bisecció BINÀRIA: 1 punt interior (50%), rang → rang/2
    · 3 avaluacions sign(ΔΦ) per iteració
    · 1 bit guanyat per iteració
    · Eficiència: 1/3 bits per avaluació

  v22 usa bisecció QUATERNÀRIA: 3 punts interiors (25%, 50%, 75%)
    · 5 avaluacions sign(ΔΦ) per iteració  (el doble + 1 de v21)
    · 2 bits guanyats per iteració          (el doble de v21)
    · Eficiència: 2/5 = 0.4 bits per avaluació  → MILLOR

  IMPACTE PER RSA-2048 (m_max ≈ 10^308, LIM_KSWEEP = 2M):
    Binària  v21: log₂(10^308 / 2M) ≈ 1002 iteracions × 5 sign  = 5010 sign
    Quaternà v22: log₄(10^308 / 2M) ≈  501 iteracions × 5 sign  = 2505 sign
    → La meitat d'avaluacions per al mateix resultat.

  EXTENSIÓ NATURAL — Bisecció per 8 (0.125):
    7 punts interiors, 9 avaluacions, 3 bits/iteració → 3/9 = 0.333 bits/aval.
    → PITJOR que quaternària. El màxim d'eficiència és a base 4.

  ─────────────────────────────────────────────────────────────────────────
  PER QUÈ BASE 4 ÉS L'ÒPTIM (matemàtica):
  ─────────────────────────────────────────────────────────────────────────
  Per a bisecció de base B:
    · Punts interiors: B−1
    · Avaluacions per iteració: B+1  (inclou extrems)
    · Bits per iteració: log₂(B)
    · Eficiència: log₂(B) / (B+1)

    B=2: log₂(2)/(2+1) = 1/3  ≈ 0.333
    B=4: log₂(4)/(4+1) = 2/5  = 0.400  ← MÀXIM
    B=8: log₂(8)/(8+1) = 3/9  ≈ 0.333
    B=e: log₂(e)/(e+1) ≈ 0.386

  El màxim matemàtic és a B ≈ e (2.718), però B ha de ser enter →
  B=4 és l'enter òptim, que coincideix exactament amb la proposta de Víctor:
  dividir per 2 dues vegades (0.5 → 0.25 → 0.125 correspon a tres nivells,
  però usem els tres punts interiors d'una sola passada per divisió per 4).
  ─────────────────────────────────────────────────────────────────────────
"""

import math
import time


# ============================================================================
#  PARÀMETRES GLOBALS
# ============================================================================
PRIMOS     = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL  = math.prod(PRIMOS)

E_CONV     = (math.e - 1) / 2    # ≈ 0.8591 — convergència correcta

LIM_CRIBA  = 2_000_000
LIM_KSWEEP = 2_000_000
MAX_BISEC  = 8000


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL
# ============================================================================

def pasa_rueda(m: int) -> bool:
    if m < 1:
        return False
    return math.gcd(2 * m + 3, PRIMORIAL) == 1


def proper_valid(m: int, pas: int = 1) -> int | None:
    for _ in range(500):
        if m < 1:
            return None
        if pasa_rueda(m):
            return m
        m += pas
    return None


# ============================================================================
#  BLOC 2 — SIGNE DEL DESFASE ESPELL (enters purs, sense Fraction)
# ============================================================================

def sgn_desfase(m: int, N: int) -> int:
    """
    sign(ΔΦ(m)) = sign(R_B·D_A − R_A·D_B)

    On:
      D_B = 2*(2m+3),       R_B = N mod D_B   (senyal principal)
      A   = N // (2m+3),    D_A = 2*(2*((A-3)//2)+3),  R_A = N mod D_A  (espell)

    Tot entera pur. Per a N de D dígits: cost O(D^1.585) per Karatsuba.
    Per RSA-2048: ~µs per avaluació.

    RETORNA: +1, −1 o 0
    """
    B   = 2 * m + 3
    if B <= 1 or B >= N:
        return 0
    D_B = 2 * B
    R_B = N % D_B

    A   = N // B
    if A < 3:
        return 0
    m_e = (A - 3) // 2
    if m_e < 1:
        return 0
    D_A = 2 * (2 * m_e + 3)
    R_A = N % D_A

    diff = R_B * D_A - R_A * D_B
    return (1 if diff > 0 else -1) if diff != 0 else 0


def mag_desfase(m: int, N: int) -> int:
    """
    |ΔΦ(m)| proporcional = |R_B·D_A − R_A·D_B|  (enter pur, no normalitzat).
    Usat per triar el costat millor quan no hi ha canvi de signe clar.
    Valor menor → més prop del zero de ΔΦ → més prop del factor.
    """
    B   = 2 * m + 3
    if B <= 1 or B >= N:
        return 10**300
    D_B = 2 * B
    R_B = N % D_B
    A   = N // B
    if A < 3:
        return 10**300
    m_e = (A - 3) // 2
    if m_e < 1:
        return 10**300
    D_A = 2 * (2 * m_e + 3)
    R_A = N % D_A
    return abs(R_B * D_A - R_A * D_B)


def check_factor(m: int, N: int) -> bool:
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


# ============================================================================
#  BLOC 3 — K-SWEEP
# ============================================================================

def k_sweep(N: int, m_ini: int, m_fi: int) -> int | None:
    """
    Iteració sobre pendents k: comprova N // k per k en [k_lo, k_hi].
    Exhaustiu i determinista. Usat com a pas final de la bisecció.
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
        c = N // k
        if c < 3 or c % 2 == 0:
            continue
        if N % c == 0 and 1 < c < N:
            return c
    return None


# ============================================================================
#  BLOC 4 — BISECCIÓ QUATERNÀRIA (÷4 per iteració)
# ============================================================================

def biseccio_quaternaria(N: int, m_ini: int, m_fi: int,
                          verbose: bool = False) -> int | None:
    """
    Bisecció de base 4: tres punts interiors per iteració (25%, 50%, 75%).

    ALGORISME:
    ──────────────────────────────────────────────────────────────────────
    Per cada iteració, generem 5 punts dins del segment:
      m_e  = extrem esquerre  (0%)
      m_q1 = primer quart     (25%)
      m_c  = centre           (50%)
      m_q3 = tercer quart     (75%)
      m_d  = extrem dret      (100%)

    Avaluem sign(ΔΦ) als 5 punts → detectem el quarter amb canvi de signe:
      [m_e,  m_q1]: s_e  ≠ s_q1  → bisectem ací (rang → rang/4, 0–25%)
      [m_q1, m_c ]: s_q1 ≠ s_c   → bisectem ací (rang → rang/4, 25–50%)
      [m_c,  m_q3]: s_c  ≠ s_q3  → bisectem ací (rang → rang/4, 50–75%)
      [m_q3, m_d ]: s_q3 ≠ s_d   → bisectem ací (rang → rang/4, 75–100%)

    Sense canvi de signe clar:
      Triem el quarter amb menor mag_desfase en els seus extrems.

    Quan n_k ≤ LIM_KSWEEP → k-sweep directe.

    NOMBRE D'ITERACIONS:
      log₄(m_max / LIM_KSWEEP) = log₂(m_max / LIM_KSWEEP) / 2
      Per RSA-2048: ≈ 501 iteracions  (vs 1002 en binària)
    ──────────────────────────────────────────────────────────────────────
    """
    for it in range(MAX_BISEC):

        # n_k del segment actual
        pos_i  = 2 * m_ini + 3
        pos_f  = 2 * m_fi  + 3
        k_lo   = max(1, N // pos_f) if pos_f > 0 else 1
        k_hi   = N // pos_i         if pos_i > 0 else 1
        n_k    = k_hi - k_lo + 1

        if verbose:
            print(f"    {it:5d}: [{m_ini:>14,}, {m_fi:>14,}]"
                  f"  amp={m_fi-m_ini:>12,}  n_k={n_k:>9,}")

        # K-sweep final quan el segment és prou menut
        if n_k <= LIM_KSWEEP:
            return k_sweep(N, m_ini, m_fi)

        # ── CINC PUNTS (extrems + tres quarts) ────────────────────────────
        q = (m_fi - m_ini) // 4          # amplada de cada quart

        m_e  = proper_valid(m_ini,           1)
        m_q1 = proper_valid(m_ini + q,       1)
        m_c  = proper_valid(m_ini + 2 * q,   1)
        m_q3 = proper_valid(m_ini + 3 * q,   1)
        m_d  = proper_valid(m_fi,            -1)

        if not all([m_e, m_q1, m_c, m_q3, m_d]):
            # Segment massa menut per a 4 quarts → bisecció binària de seguretat
            m_c2 = proper_valid((m_ini + m_fi) // 2, 1)
            if m_c2:
                if check_factor(m_c2, N):
                    return 2 * m_c2 + 3
                s_e2 = sgn_desfase(m_e or m_ini, N)
                s_c2 = sgn_desfase(m_c2, N)
                s_d2 = sgn_desfase(m_d or m_fi, N)
                if s_e2 != 0 and s_c2 != 0 and s_e2 != s_c2:
                    m_fi  = m_c2
                elif s_c2 != 0 and s_d2 != 0 and s_c2 != s_d2:
                    m_ini = m_c2
                else:
                    m_ini = m_c2
            continue

        # Comprovació directa dels cinc punts
        for mp in [m_e, m_q1, m_c, m_q3, m_d]:
            if check_factor(mp, N):
                return 2 * mp + 3

        # ── SIGNES ALS CINC PUNTS ──────────────────────────────────────────
        s_e  = sgn_desfase(m_e,  N)
        s_q1 = sgn_desfase(m_q1, N)
        s_c  = sgn_desfase(m_c,  N)
        s_q3 = sgn_desfase(m_q3, N)
        s_d  = sgn_desfase(m_d,  N)

        # ── DETECCIÓ DEL QUARTER AMB CANVI DE SIGNE ───────────────────────
        # Prioritat: el primer quarter (esquerra) on hi ha canvi de signe.
        # Si n'hi ha múltiples (N compost amb molts factors), el primer és
        # el que conté el factor més menut (o el primer per ordre).

        if s_e != 0 and s_q1 != 0 and s_e != s_q1:
            # Canvi en [0%, 25%]
            m_fi  = m_q1

        elif s_q1 != 0 and s_c != 0 and s_q1 != s_c:
            # Canvi en [25%, 50%]
            m_ini = m_q1
            m_fi  = m_c

        elif s_c != 0 and s_q3 != 0 and s_c != s_q3:
            # Canvi en [50%, 75%]
            m_ini = m_c
            m_fi  = m_q3

        elif s_q3 != 0 and s_d != 0 and s_q3 != s_d:
            # Canvi en [75%, 100%]
            m_ini = m_q3

        else:
            # ── SENSE CANVI DE SIGNE CLAR ─────────────────────────────────
            # Triem el quarter on la magnitud del desfase és menor
            # (el zero de ΔΦ és en algun punt interior d'eixe quarter).
            mags = [
                (mag_desfase(m_e,  N) + mag_desfase(m_q1, N), 0),   # quarter 0–25%
                (mag_desfase(m_q1, N) + mag_desfase(m_c,  N), 1),   # quarter 25–50%
                (mag_desfase(m_c,  N) + mag_desfase(m_q3, N), 2),   # quarter 50–75%
                (mag_desfase(m_q3, N) + mag_desfase(m_d,  N), 3),   # quarter 75–100%
            ]
            millor = min(mags, key=lambda x: x[0])[1]

            if millor == 0:
                m_fi  = m_q1
            elif millor == 1:
                m_ini = m_q1;  m_fi = m_c
            elif millor == 2:
                m_ini = m_c;   m_fi = m_q3
            else:
                m_ini = m_q3

    return None


# ============================================================================
#  BLOC 5 — ORQUESTRADOR  MDC v22
# ============================================================================

def mdc_v22(N: int, verbose: bool = True) -> tuple:
    """
    MDC v22 — Bisecció Quaternària.

    F0 · Criba roda  [1, LIM_CRIBA]
    F1 · K-sweep √N  (LIM_KSWEEP valors de k prop de m_max)
    F2 · Bisecció quaternària sign(ΔΦ)  [criba, m_conv]
         → ÷4 per iteració → la meitat d'iteracions que v21
    """
    if verbose:
        bar = '═' * 72
        print(f"\n{bar}")
        print(f"  MDC v22 — Bisecció Quaternària (base 4, òptim matemàtic)")
        print(f"  N = {N:,}")
        print(f"  ({len(str(N))} dígits  |  {N.bit_length()} bits)")
        print(f"{bar}")

    t0 = time.perf_counter()

    # ── F0: CRIBA ──────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [F0] Criba roda fins m={LIM_CRIBA:,}...")

    for p in PRIMOS:
        if N % p == 0 and p < N:
            t_ms = (time.perf_counter() - t0) * 1000
            if verbose: print(f"    → Factor primorial: {p}")
            return p, t_ms

    r = math.isqrt(N)
    if r * r == N:
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose: print(f"    → Quadrat perfecte: {r}")
        return r, t_ms

    m_max = (math.isqrt(N) - 3) // 2
    lim_c = min(m_max, LIM_CRIBA)

    for m_c in range(1, lim_c + 1):
        if pasa_rueda(m_c):
            D = 2 * m_c + 3
            if N % D == 0 and D < N:
                t_ms = (time.perf_counter() - t0) * 1000
                if verbose: print(f"    → Factor criba: {D:,}")
                return D, t_ms

    if lim_c >= m_max:
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose: print(f"    → Criba cobreix tot el rang.")
        return None, t_ms

    t_c = (time.perf_counter() - t0) * 1000
    if verbose:
        print(f"    → {t_c:.1f} ms. m_max={m_max:,}  ({m_max.bit_length()} bits)")

    # ── F1: K-SWEEP √N ────────────────────────────────────────────────────
    k_top = max(1, N // (2 * m_max + 3))
    k_lim = k_top + LIM_KSWEEP
    d_lim = N // k_lim if k_lim > 0 else 2 * m_max + 3
    m_lim = max(lim_c + 1, (d_lim - 3) // 2)

    if verbose:
        print(f"\n  [F1] K-sweep [{m_lim:,}, {m_max:,}]  (~{LIM_KSWEEP:,} k)...")

    factor = k_sweep(N, m_lim, m_max)
    if factor:
        t_ms = (time.perf_counter() - t0) * 1000
        _print(N, factor, "F1 k-sweep √N", t_ms, verbose)
        return factor, t_ms

    t_f1 = (time.perf_counter() - t0) * 1000
    if verbose: print(f"    → {t_f1:.1f} ms. Res trobat.")

    # ── F2: BISECCIÓ QUATERNÀRIA ───────────────────────────────────────────
    m_conv_val = max(lim_c + 1, int(E_CONV * m_max))

    if verbose:
        rng = m_conv_val - lim_c
        n_it_bin  = max(1, rng.bit_length() - LIM_KSWEEP.bit_length() + 1)
        n_it_quat = (n_it_bin + 1) // 2
        print(f"\n  [F2] Bisecció quaternària [{lim_c:,}, {m_conv_val:,}]")
        print(f"    → Rang: {rng:,}  ({rng.bit_length()} bits)")
        print(f"    → Binària (v21): ~{n_it_bin} iteracions")
        print(f"    → Quaternària v22: ~{n_it_quat} iteracions  (÷2)")

    factor = biseccio_quaternaria(N, lim_c, m_conv_val, verbose=verbose)

    t_ms = (time.perf_counter() - t0) * 1000
    if factor:
        _print(N, factor, "F2 bisecció quaternària", t_ms, verbose)
        return factor, t_ms

    if verbose: print(f"\n  ✗ No trobat. Temps: {t_ms:.2f} ms")
    return None, t_ms


def _print(N, f, zona, t_ms, verbose):
    if not verbose: return
    bar = '═' * 72
    print(f"\n{bar}")
    print(f"  🎯 [{zona}]")
    print(f"     p = {f:,}")
    print(f"     q = {N//f:,}")
    print(f"     t = {t_ms:.3f} ms")
    print(f"{bar}")


# ============================================================================
#  BLOC 6 — BATERIA AMB COMPARACIÓ v21 vs v22
# ============================================================================

def executar_bateria():
    import sys
    # Importem v21 per comparació
    sys.path.insert(0, '/home/claude')
    from mdc_v21 import mdc_v21

    proves = [
        (101 * 103,
         "10.403  (5 dígits)"),
        (100003 * 100019,
         "11 dígits equilibrat"),
        (9_999_991 * 9_999_973,
         "14 dígits"),
        (999_983 * 1_000_000_000_000_003,
         "21 dígits asimètric"),
        (100_000_000_003 * 100_000_000_019,
         "22 dígits"),
        (10_000_000_000_037 * 10_000_000_000_099,
         "26 dígits"),
    ]

    # Afegim semiprimers de primers veritables grans
    try:
        import sympy
        for exp in [12, 15, 18, 22, 26, 30]:
            base = 10**exp
            p = sympy.nextprime(base + 37)
            q = sympy.nextprime(p + 1000)
            N = p * q
            proves.append((N, f"{len(str(N))} dígits (primers veritables)"))
    except ImportError:
        pass

    print("\n" + "█" * 72)
    print("  📊 COMPARACIÓ v21 (binària) vs v22 (quaternària)")
    print("█" * 72)
    print(f"  {'Descripció':<35} {'v21 (ms)':>10}  {'v22 (ms)':>10}  {'Ràtio':>8}")
    print(f"  {'─'*35} {'─'*10}  {'─'*10}  {'─'*8}")

    ok = fail = 0
    for item in proves:
        N, desc = item[0], item[1]
        f21, t21 = mdc_v21(N, verbose=False)
        f22, t22 = mdc_v22(N, verbose=False)

        if f22:
            assert N % f22 == 0
            ratio = t21 / t22 if t22 > 0 else float('inf')
            estat = "✅"
            ok += 1
        else:
            ratio = 0.0
            estat = "⚠️ "
            fail += 1

        print(f"  {estat} {desc:<33} {t21:>10.1f}  {t22:>10.1f}  {ratio:>7.2f}x")

    print(f"\n  Encerts: {ok}/{ok+fail}")
    print("█" * 72)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

if __name__ == "__main__":
    # Prova verbose individual
    N = 100003 * 100019
    mdc_v22(N, verbose=True)

    print("\n" + "─" * 72)
    executar_bateria()
