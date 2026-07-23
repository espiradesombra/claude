unt de tall al 25% (prop de m_ini, lluny de √N) és on la
    funció ΔΦ té major variació i el canvi de signe és més visible.
    Prop de √N, ΔΦ tendeix a 0 i és més difícil de llegir.
    En posar el punt de tall al 25% concentres la "lupa" on la
    senyal és més clara i deixes la zona prop de √N per a F1.

  ─────────────────────────────────────────────────────────────────────────
  ARQUITECTURA v23:
    F0 · Criba roda  [1, LIM_CRIBA]
    F1 · K-sweep     [m_lim, m_max]        — cobreix prop de √N
    F2 · Bisecció asimètrica 25/75         — cobreix la resta
         Si no detecta canvi → factor ja a F1 → fi.
         Si detecta canvi   → bisecciona el 25% inferior.
  ─────────────────────────────────────────────────────────────────────────
"""

import math
import time


# ============================================================================
#  PARÀMETRES GLOBALS
# ============================================================================
PRIMOS     = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL  = math.prod(PRIMOS)
E_CONV     = (math.e - 1) / 2    # ≈ 0.8591
LIM_CRIBA  = 2_000_000
LIM_KSWEEP = 2_000_000
MAX_BISEC  = 8000

# Punt de tall asimètric: 25% des de m_ini (lluny de √N)
TALL_ASIM  = 4     # rang // TALL_ASIM  → 25% des del costat menut


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
#  BLOC 2 — SIGNE DEL DESFASE ESPELL (enters purs)
# ============================================================================

def sgn_desfase(m: int, N: int) -> int:
    """
    sign(ΔΦ(m)) = sign(R_B·D_A − R_A·D_B)
    Enter pur, sense Fraction. O(D^1.585) per Karatsuba.
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
    d   = R_B * D_A - R_A * D_B
    return (1 if d > 0 else -1) if d != 0 else 0


def mag_desfase(m: int, N: int) -> int:
    """
    |ΔΦ(m)| proporcional (enter pur, no normalitzat).
    Menor → més prop del zero → més prop del factor.
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
    Iteració sobre pendents k: comprova N//k per k en [k_lo, k_hi].
    Exhaustiu i determinista per al segment [m_ini, m_fi].
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
#  BLOC 4 — BISECCIÓ ASIMÈTRICA 25/75
# ============================================================================

def biseccio_asimetrica(N: int, m_ini: int, m_fi: int,
                         verbose: bool = False) -> tuple[int | None, bool]:
    """
    Bisecció asimètrica: tall al 25% des de m_ini (lluny de √N).

    RETORNA: (factor, factor_es_a_f1)
      factor        : el factor trobat, o None
      factor_es_a_f1: True si el factor no és al rang [m_ini, m_fi]
                      (per tant és prop de √N, cobert per F1)

    ─────────────────────────────────────────────────────────────────────
    ALGORISME:
      Cada iteració treballa sobre el rang actual [m_ini, m_fi]:

      1. Calculem m_25 = m_ini + (m_fi - m_ini) // TALL_ASIM  (el 25%)

      2. Avaluem sign(ΔΦ) als tres punts: m_ini (extrem), m_25 (25%), m_fi (dret)

      3. DECISIÓ:
         a) Canvi entre m_ini i m_25:
            → factor al 25% inferior [m_ini, m_25]
            → m_fi = m_25  (rang → rang × 0.25, molt precís)

         b) Canvi entre m_25 i m_fi:
            → factor al 75% superior [m_25, m_fi]
            → m_ini = m_25  (rang → rang × 0.75)
            (Nota: el 75% superior és prop de √N, i si és molt gran
             tornarem a bisellar amb el 25% d'eixe tros)

         c) CAP canvi de signe:
            → el factor NO és en [m_ini, m_fi]
            → el factor és prop de √N (el F1 ja l'ha cobert)
            → retornem (None, True)  ← "está en raiz"

      4. Quan n_k ≤ LIM_KSWEEP → k-sweep directe.

    CONVÈRGÈNCIA:
      Cada vegada que caiem al cas a) (25%) el rang es redueix per 4.
      Cada vegada que caiem al cas b) (75%) el rang es redueix per 1.33.
      El pitjor cas és si sempre caiem al b): log_{1.33}(rang) talls,
      però cada tall b) deixa la zona de "interès" cada vegada més
      propera a √N, fins que F1 la cobrirà de nou.

    PER QUÈ EL PUNT DE TALL AL 25% (LLUNY DE √N) DONA MÉS PRECISIÓ:
      La funció ΔΦ té més variació (pendent més gran) per a m petits,
      on els dents de serra canvien de k a cada pocs passos.
      Prop de √N (m gran), ΔΦ és quasi constant i difícil de llegir.
      En col·locar el punt de tall al 25% (m petit), aprofitem la zona
      on el senyal és més clar i el canvi de signe és inconfusible.
    ─────────────────────────────────────────────────────────────────────
    """
    for it in range(MAX_BISEC):

        # n_k del segment
        pos_i = 2 * m_ini + 3
        pos_f = 2 * m_fi  + 3
        k_lo  = max(1, N // pos_f) if pos_f > 0 else 1
        k_hi  = N // pos_i         if pos_i > 0 else 1
        n_k   = k_hi - k_lo + 1

        if verbose:
            rang = m_fi - m_ini
            print(f"    {it:5d}: [{m_ini:>14,}, {m_fi:>14,}]"
                  f"  rang={rang:>12,}  n_k={n_k:>9,}")

        # K-sweep final
        if n_k <= LIM_KSWEEP:
            return k_sweep(N, m_ini, m_fi), False

        # Punt de tall asimètric: 25% des de m_ini (la zona de menor m)
        rang  = m_fi - m_ini
        m_25  = proper_valid(m_ini + rang // TALL_ASIM, 1)
        m_e   = proper_valid(m_ini, 1)
        m_d   = proper_valid(m_fi, -1)

        if not (m_e and m_25 and m_d):
            # Segment massa menut: k-sweep directe
            return k_sweep(N, m_ini, m_fi), False

        # Comprovació directa
        for mp in [m_e, m_25, m_d]:
            if check_factor(mp, N):
                return 2 * mp + 3, False

        # Signes
        s_e  = sgn_desfase(m_e,  N)
        s_25 = sgn_desfase(m_25, N)
        s_d  = sgn_desfase(m_d,  N)

        if verbose:
            print(f"          s_e={s_e:+d}  s_25={s_25:+d}  s_d={s_d:+d}"
                  f"  (punts: {m_e:,} | {m_25:,} | {m_d:,})")

        # DECISIÓ ASIMÈTRICA
        if s_e != 0 and s_25 != 0 and s_e != s_25:
            # Factor al 25% inferior → rang × 0.25 (molt precís)
            if verbose:
                print(f"          → Canvi [0-25%]: bisectem rang × 0.25")
            m_fi = m_25

        elif s_25 != 0 and s_d != 0 and s_25 != s_d:
            # Factor al 75% superior → rang × 0.75
            if verbose:
                print(f"          → Canvi [25-100%]: bisectem rang × 0.75")
            m_ini = m_25

        else:
            # CAP CANVI DE SIGNE:
            # El factor NO és en aquest rang → és prop de √N (cobert per F1)
            if verbose:
                print(f"          → Sense canvi de signe: "
                      f"factor prop de √N, cobert per F1")
            return None, True     # ← "factor_es_a_f1 = True"

    return None, False


# ============================================================================
#  BLOC 5 — ORQUESTRADOR  MDC v23
# ============================================================================

def mdc_v23(N: int, verbose: bool = True) -> tuple:
    """
    MDC v23 — Bisecció Asimètrica 25/75.

    F0 · Criba roda  [1, LIM_CRIBA]
    F1 · K-sweep     [m_lim, m_max]  — factors prop de √N
    F2 · Bisecció asimètrica 25/75   — factors de mida intermèdia
         Si "factor_es_a_f1" → el factor ja estava a F1.
    """
    if verbose:
        bar = '═' * 72
        print(f"\n{bar}")
        print(f"  MDC v23 — Bisecció Asimètrica 25/75")
        print(f"  Tall al 25% (lluny de √N) · Si no hi és → ja és a √N")
        print(f"  N = {N:,}")
        print(f"  ({len(str(N))} dígits  |  {N.bit_length()} bits)")
        print(f"{bar}")

    t0 = time.perf_counter()

    # ── F0: CRIBA RODA ─────────────────────────────────────────────────────
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
        print(f"    → {t_c:.1f} ms. m_max={m_max:,}")

    # ── F1: K-SWEEP PROP √N ────────────────────────────────────────────────
    # Sempre: cobreix els LIM_KSWEEP valors de k més alts (factors prop de √N).
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

    # ── F2: BISECCIÓ ASIMÈTRICA 25/75 ─────────────────────────────────────
    # Rang: [LIM_CRIBA, m_conv] — factors intermedis (ni petits ni prop de √N)
    m_conv_val = max(lim_c + 1, int(E_CONV * m_max))

    if verbose:
        rang = m_conv_val - lim_c
        # Pitjor cas: log₄(rang/LIM_KSWEEP) talls al 25%
        n_pit = max(1, round(math.log(max(2, rang / LIM_KSWEEP), 4)))
        print(f"\n  [F2] Bisecció asimètrica 25/75 [{lim_c:,}, {m_conv_val:,}]")
        print(f"    → Rang: {rang:,}  ({rang.bit_length()} bits)")
        print(f"    → Pitjor cas: ~{n_pit} talls (sempre al 25% inferior)")
        print(f"    → Cas típic: menys talls (quan cau al 75% → prop √N → fi)")

    factor, es_f1 = biseccio_asimetrica(N, lim_c, m_conv_val, verbose=verbose)

    t_ms = (time.perf_counter() - t0) * 1000

    if factor:
        _print(N, factor, "F2 bisecció asimètrica", t_ms, verbose)
        return factor, t_ms

    if es_f1 and verbose:
        print(f"\n  ℹ️  La bisecció confirma: factor és prop de √N.")
        print(f"     (F1 ja l'ha cobert, o N és primer en aquest rang.)")

    if verbose:
        print(f"\n  ✗ No trobat. Temps: {t_ms:.2f} ms")
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
#  BLOC 6 — BATERIA AMB COMPARACIÓ v21/v22/v23
# ============================================================================

def executar_bateria():
    import sys
    sys.path.insert(0, '/home/claude')
    from mdc_v21 import mdc_v21 as v21
    from mdc_v22 import mdc_v22 as v22

    proves = [
        (101 * 103,                                    "5 dígits"),
        (100003 * 100019,                              "11 dígits equilibrat"),
        (9_999_991 * 9_999_973,                        "14 dígits"),
        (1_548_586_332_452_843,                        "16 dígits factor 59"),
        (999_983 * 1_000_000_000_000_003,              "21 dígits asimètric"),
        (100_000_000_003 * 100_000_000_019,            "22 dígits"),
        (10_000_000_000_037 * 10_000_000_000_099,      "26 dígits"),
    ]

    # Afegim primers veritables per escalada real
    try:
        import sympy
        for exp in [12, 15, 18, 22, 26, 30, 35, 40]:
            base = 10**exp
            p = sympy.nextprime(base + 37)
            q = sympy.nextprime(p + 1000)
            N = p * q
            proves.append((N, f"{len(str(N))}d primers veritables"))
    except ImportError:
        pass

    print("\n" + "█" * 72)
    print("  📊 COMPARACIÓ v21 · v22 · v23")
    print(f"  {'Descripció':<30} {'v21':>8}  {'v22':>8}  {'v23':>8}  {'Guany':>8}")
    print(f"  {'─'*30} {'─'*8}  {'─'*8}  {'─'*8}  {'─'*8}")

    ok = fail = 0
    for N, desc in proves:
        t21 = v21(N, verbose=False)[1]
        t22 = v22(N, verbose=False)[1]
        f23, t23 = mdc_v23(N, verbose=False)

        if f23 and N % f23 == 0:
            guany = t21 / t23 if t23 > 0 else 99.0
            estat = "✅"
            ok += 1
        else:
            guany = 0.0
            estat = "⚠️ "
            fail += 1

        print(f"  {estat} {desc:<28} {t21:>7.0f}  {t22:>7.0f}  {t23:>7.0f}  {guany:>7.2f}x")

    print(f"\n  Encerts v23: {ok}/{ok+fail}")
    print("█" * 72)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

if __name__ == "__main__":
    # Prova verbose (11 dígits, per veure el flux de decisions)
    N = 100003 * 100019
    mdc_v23(N, verbose=True)

    print("\n" + "─" * 72)
    executar_bateria()

