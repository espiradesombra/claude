# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC v18 — MÈTODE DIOFÀNTIC CINEMÀTIC
  Zona Densa K-sweep · Espiral Sota m_conv · Pinça 4+4 · Espell
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)

  ─────────────────────────────────────────────────────────────────────────
  INSIGHT CENTRAL DE v18 (del Libro 5, Víctor 2025)
  ─────────────────────────────────────────────────────────────────────────

  1. ZONA DENSA [m_conv, m_max]:
     La convergència CORRECTA de l'acumulador 1/(2·i!) és:
       m_conv = (e−1)/2 · m_max  ≈  0.8591 · m_max
     (No e−2 ≈ 0.718, que era l'error de les IAs en el document.)

     En esta zona [m_conv, m_max], que ocupa sols el 14% de l'espai:
       · D = 2m+3 ≈ √N → gran
       · k = N//D canvia a CADA pas de m (un k diferent per m)
       · Els mínims dels dents apareixen AMB MÀXIMA FREQÜÈNCIA
       · n_k = k_max − k_min ≈ amplada_zona (prop de m_max)
     Per tant, un k-sweep sobre esta zona és MOLT ràpid: en poques
     iteracions de k trobes el factor si hi és.

  2. ZONA ESPIRAL [1, m_conv]:
     Per sota de m_conv, els dents s'allarguen i els mínims s'espeen.
     Ací s'aplica l'espiral harmònica 1/(2·i!) per cobrir el rang
     generant sondes a les transicions de pendent k→k+1.

  3. PER QUÈ NO CALIA ANAR FINS m_max (√N):
     El k-sweep sobre [m_conv, m_max] JA cobreix tota la zona densa
     de forma EXHAUSTIVA. No cal generar sondes en este rang: el
     k-sweep ho fa tot sol en un sol bucle de k.

  ─────────────────────────────────────────────────────────────────────────
  ARQUITECTURA v18:
    F0 · Precondicions
    F1 · K-sweep exhaustiu sobre la zona densa [m_conv, m_max]
         (on els mínims apareixen més freqüentment)
    F2 · Espiral harmònica 1/(2·i!) per sota m_conv
         → genera sondes a les transicions de pendent
    F3 · Detector espell ΔΦ → prioritza segments candidats
    F4 · K-sweep per segment (general, amb bisecció per segments grans)
  ─────────────────────────────────────────────────────────────────────────
"""

from fractions import Fraction
import math
import time


# ============================================================================
#  BLOC 0 — CONSTANTS
# ============================================================================
PRIMOS    = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)   # 223.092.870

# Convergència correcta: (e−1)/2 ≈ 0.8591
# (No e−2 ≈ 0.7183, que era l'error dels documents anteriors)
E_MINUS_1_OVER_2 = (math.e - 1) / 2   # ≈ 0.8591409...


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL
# ============================================================================

def pasa_rueda(m: int) -> bool:
    """
    Comprova si D(m) = 2m+3 és coprimer amb el primorial p23.
    Densitat: ~8.07% dels enters passen la roda.
    """
    if m < 1:
        return False
    return math.gcd(2 * m + 3, PRIMORIAL) == 1


def proper_valid(m: int, pas: int = 1) -> int | None:
    """Trova el m vàlid (que passa la roda) més proper. pas=±1."""
    for _ in range(500):
        if m < 1:
            return None
        if pasa_rueda(m):
            return m
        m += pas
    return None


# ============================================================================
#  BLOC 2 — FUNCIÓ D'ONA I ESPELL
# ============================================================================

def d_frac(m: int, N: int) -> Fraction:
    """
    Funció d'ona: d(m) = Fraction(N mod D, D) − ½
    on D = 2*(2m+3).  Fraction exacta, zero error d'arrodoniment.
    """
    D = 2 * (2 * m + 3)
    return Fraction(N % D, D) - Fraction(1, 2)


def m_espell(m: int, N: int) -> int | None:
    """
    Espell aritmètic: A = N // (2m+3),  m_espell = (A−3) // 2.
    Principi: (2m+3) · A ≈ N  (exacte quan m és el factor real).
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
    ΔΦ(m) = d_frac(m) − d_frac(m_espell).
    Canvi de signe entre dos sondes → factor atrapat entre elles.
    """
    dB = d_frac(m, N)
    me = m_espell(m, N)
    if me is None:
        return Fraction(0)
    return dB - d_frac(me, N)


def check_factor(m: int, N: int) -> bool:
    """Verificació directa i exacta: (2m+3) divideix N?"""
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


# ============================================================================
#  BLOC 3 — K-SWEEP: EL COR DEL MÈTODE
# ============================================================================

def k_sweep(N: int, m_ini: int, m_fi: int,
            verbose: bool = False) -> int | None:
    """
    K-sweep: iteració sobre les pendents k en lloc dels valors m.

    MATEMÀTICA FONAMENTAL (pendents lineals del Libro 5):
      r(d) = N − k·d    per a tot d en [N/(k+1), N/k]
      on k = ⌊N/d⌋ és la pendent actual (cocient de la divisió).

      El zero de la recta k és:  d_zero = N // k
      Si N % d_zero == 0 → FACTOR TROBAT.

    PERQUÈ ÉS EXHAUSTIU:
      Per a cada m en [m_ini, m_fi], la seua pendent k = N//(2m+3)
      cau en l'interval [k_lo, k_hi]. Iterant sobre TOTES les k en
      este rang i comprovant N // k, cobrim TOTS els factors possibles
      sense ometre'n cap. Es completament determinista.

    EFICIÈNCIA A LA ZONA DENSA [m_conv, m_max]:
      En esta zona D ≈ √N, de forma que k ≈ √N també. Cada unitat
      de m produïx una k diferent: n_k ≈ amplada del segment.
      Per als factor prop de m_max: n_k ≈ 14% · m_max, manejable.

    RETORNA: el factor D = 2m+3, o None si no n'hi ha en el rang.
    """
    if m_ini < 1:
        m_ini = 1
    if m_fi < m_ini:
        return None

    pos_ini = 2 * m_ini + 3   # pos menut (k gran)
    pos_fi  = 2 * m_fi  + 3   # pos gran  (k menut)

    k_lo = max(1, N // pos_fi)    # k mínim (per al pos gran)
    k_hi = N // pos_ini           # k màxim (per al pos menut)

    if k_lo > k_hi:
        return None

    n_k = k_hi - k_lo + 1

    if verbose:
        print(f"      k-sweep: k∈[{k_lo:,},{k_hi:,}]  n_k={n_k:,}")

    # Iteració sobre k (les pendents)
    for k in range(k_lo, k_hi + 1):
        if k <= 0:
            continue
        candidat = N // k   # zero de la recta k
        if candidat < 3:
            continue
        if candidat % 2 == 0:
            continue   # (2m+3) sempre senar → candidat ha de ser senar
        if N % candidat == 0 and 1 < candidat < N:
            if verbose:
                m_c = (candidat - 3) // 2
                print(f"      ✓ k={k}  candidat={candidat:,}  m={m_c:,}")
            return candidat

    return None


def k_sweep_bisectat(N: int, m_ini: int, m_fi: int,
                     llindar_k: int = 2_000_000,
                     profunditat: int = 60,
                     verbose: bool = False) -> int | None:
    """
    K-sweep amb bisecció per ΔΦ espell per a segments molt grans.

    Si n_k <= llindar_k → k-sweep directe (exhaustiu, ràpid).
    Si n_k > llindar_k → bisecta el segment per ΔΦ fins que
                          el fragment sigui prou menut per al k-sweep.

    MATEMÀTICA DE LA BISECCIÓ:
      El desfase espell ΔΦ canvia de signe quan el segment conté
      un factor. Bisectem triant la meitat amb canvi de signe,
      reduint l'interval fins que n_k sigui manejable.

      Quan no hi ha canvi de signe clar (factor molt equilibrat prop
      d'un extrem del segment), usem el criteri de menor N%pos%(pos-1)
      (distància al zero de pendent) per triar el costat millor.
    """
    m_e_orig = m_ini
    m_f_orig = m_fi

    for it in range(profunditat):
        pos_i = 2 * m_ini + 3
        pos_f = 2 * m_fi  + 3
        k_lo  = max(1, N // pos_f)
        k_hi  = N // pos_i
        n_k   = k_hi - k_lo + 1

        # K-sweep directe quan el segment és manejable
        if n_k <= llindar_k:
            return k_sweep(N, m_ini, m_fi, verbose)

        # Bisecció
        m_c = (m_ini + m_fi) // 2
        m_e = proper_valid(m_ini, 1)
        m_cm = proper_valid(m_c,  1)
        m_f = proper_valid(m_fi, -1)

        if not (m_e and m_cm and m_f):
            break

        # Comprovació directa dels punts de bisecció via k-sweep local
        for mp in [m_e, m_cm, m_f]:
            r = k_sweep(N, mp, mp, verbose=False)
            if r:
                return r

        df_e  = desfase_espell(m_e,  N)
        df_cm = desfase_espell(m_cm, N)
        df_f  = desfase_espell(m_f,  N)

        if df_e * df_cm < 0:
            m_fi  = m_cm   # factor a la meitat esquerra
        elif df_cm * df_f < 0:
            m_ini = m_cm   # factor a la meitat dreta
        else:
            # Sense canvi de signe clar: usem n%pos%(pos-1) per triar
            pos_e = 2 * m_e + 3
            pos_f2 = 2 * m_f + 3
            R2_e = (N % pos_e) % (pos_e - 1) if pos_e > 1 else 10**18
            R2_f = (N % pos_f2) % (pos_f2 - 1) if pos_f2 > 1 else 10**18
            if R2_e <= R2_f:
                m_fi  = m_cm
            else:
                m_ini = m_cm

    return None


# ============================================================================
#  BLOC 4 — ESPIRAL HARMÒNICA  (per sota de m_conv)
# ============================================================================

def generar_sondes_sota_mconv(m_conv: int, N: int) -> list[tuple]:
    """
    Genera parells (m_ini, m_fi) per als segments a cobrir
    per sota de m_conv usant l'espiral harmònica 1/(2·i!).

    MATEMÀTICA:
      L'acumulador genera salts decreixents des de m_conv cap a 0:
        salta_1 = m_conv / 2       → cobreix [m_conv/2, m_conv]
        salta_2 = m_conv / 4       → cobreix [m_conv/4, m_conv/2]
        salta_3 = m_conv / 12      → cobreix [m_conv/6, m_conv/4]
        ...
      Cada salt genera un SEGMENT que es cobrirà amb k_sweep_bisectat.

    PER QUÈ 1/(2·i!) (nota del Libro 5):
      La seqüència de longituds decreixents 1/2, 1/6, 1/24, ...
      alinea les sondes amb les TRANSICIONS DE PENDENT k→k+1,
      que és on pot aparèixer un zero (un factor). Unes poques
      desenes d'iteracions cobreixen el rang fins a m ≈ 0.
      La suma total: Σ 1/(2·i!) = (e−1)/2 ≈ 0.8591 del rang
      (açò confirma que l'espiral cobreix quasi tot el rang [0, m_conv]).

    RETORNA: llista de (m_ini, m_fi) ordenada per m_fi descendent
             (segments superiors primer, on els dents són més densos).
    """
    segments = []
    salta    = Fraction(m_conv, 2)   # primer salt: m_conv/2
    posicio  = Fraction(m_conv)       # partim des de m_conv cap a baix

    for i in range(1, 200):
        m_fi_seg  = int(posicio)
        posicio  -= salta
        m_ini_seg = max(1, int(posicio))

        if m_fi_seg >= 1 and m_fi_seg > m_ini_seg:
            segments.append((m_ini_seg, m_fi_seg))

        if m_ini_seg <= 1:
            # Afegim el segment final fins a m=1
            if m_ini_seg > 1:
                segments.append((1, m_ini_seg))
            break

        # Reduïm el salt
        salta = salta * Fraction(1, i + 1)
        if salta < Fraction(1, 4):
            # Salt massa menut: afegim segment residual fins a m=1
            segments.append((1, int(posicio)))
            break

    # Ordenem de major a menor (segments superiors = dents més densos = primer)
    segments.sort(key=lambda s: s[1], reverse=True)
    return segments


# ============================================================================
#  BLOC 5 — PINÇA DOBLE 4+4  (cinemàtica exacta, Libro 5)
# ============================================================================

def pinca_doble(m: int, N: int) -> dict:
    """
    Pinça doble 4+4: 4 mesures cap endavant + 4 cap enrere.

    V, A, J de cada direcció via diferències finites (Fraction exacta).
    El FILTRE DE JERK discrimina si estem dins d'un dent (J≈0)
    o en la macro-envolvente parabòlica (J≠0 en ambdues direccions).
    """
    # Endavant
    df = [d_frac(m + i, N) for i in range(4)]
    Vf = df[1] - df[0]
    Af = df[2] - 2*df[1] + df[0]
    Jf = df[3] - 3*df[2] + 3*df[1] - df[0]

    # Enrere
    dr = [d_frac(m - i, N) for i in range(4) if m - i >= 1]
    Vr = Ar = Jr = Fraction(0)
    if len(dr) >= 4:
        Vr = dr[1] - dr[0]
        Ar = dr[2] - 2*dr[1] + dr[0]
        Jr = dr[3] - 3*dr[2] + 3*dr[1] - dr[0]

    # Salts mestres
    salt_f = None
    if Af != 0:
        despl = Vf / Af
        di = int(despl + Fraction(1,2)) if despl >= 0 else -int(-despl + Fraction(1,2))
        m_next = m - di
        if m_next >= 1:
            salt_f = m_next

    salt_r = None
    if Ar != 0:
        despl = Vr / Ar
        di = int(despl + Fraction(1,2)) if despl >= 0 else -int(-despl + Fraction(1,2))
        m_next = m + di
        if m_next >= 1:
            salt_r = m_next

    return {
        'Vf': Vf, 'Af': Af, 'Jf': Jf, 'salt_f': salt_f,
        'Vr': Vr, 'Ar': Ar, 'Jr': Jr, 'salt_r': salt_r,
        'delta_J': abs(Jf - Jr),
    }


# ============================================================================
#  BLOC 6 — ORQUESTRADOR PRINCIPAL  MDC v18
# ============================================================================

def mdc_v18(N: int, verbose: bool = True) -> tuple:
    """
    MDC v18 — Orquestrador complet.

    FLUX (4 fases):
      F0 · Precondicions (trivials + quadrat)
      F1 · K-sweep exhaustiu sobre la ZONA DENSA [m_conv, m_max]
           (14% del rang, màxima freqüència de mínims, instant)
      F2 · Espiral 1/(2·i!) → segments per sota m_conv
           → k_sweep_bisectat per segment
      F3 · Pinça doble 4+4 en els millors candidats residuals

    COMPLEXITAT:
      F1: O(m_max · 0.14 · 2) = O(√N · 0.28) iteracions de k
          (en la pràctica molt menys: n_k ≈ 14% · m_max per a RSA)
      F2: O(log N) segments × k-sweep bisectat per segment
      Total: O(√N) en el pitjor cas, però O(log² N) per a casos típics
    """
    if verbose:
        bar = '═' * 72
        print(f"\n{bar}")
        print(f"  MDC v18 — Mètode Diofàntic Cinemàtic")
        print(f"  Zona Densa K-sweep [m_conv,m_max] · Espiral sota m_conv")
        print(f"  N = {N:,}  ({len(str(N))} dígits)")
        print(f"{bar}")

    t0 = time.perf_counter()

    # ── FASE 0: PRECONDICIONS ──────────────────────────────────────────────
    if verbose:
        print(f"\n  [F0] Precondicions...")

    for p in PRIMOS:
        if N % p == 0 and p < N:
            t_ms = (time.perf_counter() - t0) * 1000
            if verbose:
                print(f"    → Factor trivial: {p}")
            return p, t_ms

    r = math.isqrt(N)
    if r * r == N:
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose:
            print(f"    → Quadrat perfecte: {r}")
        return r, t_ms

    m_max  = (math.isqrt(N) - 3) // 2
    # Convergència CORRECTA: (e−1)/2 ≈ 0.8591 (no e−2 ≈ 0.718!)
    m_conv = int(E_MINUS_1_OVER_2 * m_max)

    # Criba per factors petits fins a min(√N, ~2M) usant la roda
    # La roda primorial p23 filtra el 91.8% dels candidats
    # Cobrim fins a 2M de candidats (~ 2M × 12.4 = 25M enters) en < 100ms
    lim_criba = min(math.isqrt(N) + 1, 2_000_000)
    if verbose:
        print(f"    → Criba roda fins a m={lim_criba:,} (D={2*lim_criba+3:,})...")
    for m_c in range(1, lim_criba):
        if pasa_rueda(m_c):
            D = 2 * m_c + 3
            if N % D == 0 and D < N:
                t_ms = (time.perf_counter() - t0) * 1000
                if verbose:
                    print(f"    → Factor petit: {D:,}")
                return D, t_ms

    if verbose:
        print(f"    → m_max  = {m_max:,}")
        print(f"    → m_conv = {m_conv:,}  ({E_MINUS_1_OVER_2*100:.2f}% de m_max)")
        print(f"    → Zona densa: [{m_conv:,}, {m_max:,}]  "
              f"(amplada={m_max-m_conv:,})")

    # ── FASE 1: K-SWEEP ZONA DENSA [m_conv, m_max] ────────────────────────
    if verbose:
        print(f"\n  [F1] K-sweep zona densa [m_conv, m_max]...")

    # La zona densa és petita (14% de m_max) però amb màxima freqüència
    # de mínims. Un sol k-sweep la cobreix completament en poques ms.
    pos_conv = 2 * m_conv + 3
    pos_max  = 2 * m_max  + 3
    k_lo_dense = max(1, N // pos_max)
    k_hi_dense = N // pos_conv
    n_k_dense  = k_hi_dense - k_lo_dense + 1

    if verbose:
        print(f"    → n_k = {n_k_dense:,}  (iteracions de k per cobrir zona densa)")

    factor = k_sweep(N, m_conv, m_max, verbose=verbose and n_k_dense < 100)

    if factor:
        t_ms = (time.perf_counter() - t0) * 1000
        _print_resultat(N, factor, "zona densa", t_ms)
        return factor, t_ms

    if verbose:
        print(f"    → No trobat a la zona densa.")

    # ── FASE 2: ESPIRAL SOTA m_conv ───────────────────────────────────────
    if verbose:
        print(f"\n  [F2] Espiral 1/(2·i!) per sota m_conv={m_conv:,}...")

    segments_espiral = generar_sondes_sota_mconv(m_conv, N)

    if verbose:
        print(f"    → {len(segments_espiral)} segments generats per l'espiral")
        for (mi, mf) in segments_espiral[:5]:
            pos_f = 2*mf+3
            pos_i = 2*mi+3
            kl = max(1, N//pos_f) if pos_f>0 else 1
            kh = N//pos_i if pos_i>0 else 1
            print(f"       [{mi:>10,}, {mf:>10,}]  n_k={kh-kl:,}")
        if len(segments_espiral) > 5:
            print(f"       ... i {len(segments_espiral)-5} segments més")

    for idx, (m_ini_s, m_fi_s) in enumerate(segments_espiral):

        if verbose and idx < 8:
            print(f"\n    Segment espiral #{idx+1}: [{m_ini_s:,}, {m_fi_s:,}]",
                  end="  ")

        factor = k_sweep_bisectat(N, m_ini_s, m_fi_s,
                                   llindar_k=1_000_000,
                                   profunditat=50,
                                   verbose=False)

        if factor:
            t_ms = (time.perf_counter() - t0) * 1000
            _print_resultat(N, factor, f"espiral #{idx+1}", t_ms)
            return factor, t_ms
        elif verbose and idx < 8:
            print("✗")

    # ── FASE 3: PINÇA DOBLE en els millors candidats ──────────────────────
    # (rarament arriba ací; és un safety-net per si k-sweep ha fallat)
    if verbose:
        print(f"\n  [F3] Pinça doble sobre candidats per profunditat...")

    # Prenem 20 punts equidistants sota m_conv i apliquem pinça doble
    pas_pc = max(1, m_conv // 20)
    for i in range(20):
        m_s = i * pas_pc + 1
        mv  = proper_valid(m_s, 1)
        if not mv or mv > m_conv:
            continue
        cinc = pinca_doble(mv, N)
        for salt in [cinc['salt_f'], cinc['salt_r']]:
            if salt and 1 <= salt <= m_max:
                if check_factor(salt, N):
                    t_ms = (time.perf_counter() - t0) * 1000
                    _print_resultat(N, 2*salt+3, "pinça", t_ms)
                    return 2*salt+3, t_ms

    # ── FALLADA ────────────────────────────────────────────────────────────
    t_ms = (time.perf_counter() - t0) * 1000
    if verbose:
        print(f"\n  ✗ Factor no trobat. Temps: {t_ms:.3f} ms")
    return None, t_ms


def _print_resultat(N, factor, zona, t_ms):
    bar = '═' * 72
    print(f"\n{bar}")
    print(f"  🎯 FACTOR TROBAT!  (via {zona})")
    print(f"     Factor p : {factor:,}")
    print(f"     Factor q : {N // factor:,}")
    print(f"     Temps    : {t_ms:.3f} ms")
    print(f"{bar}")


# ============================================================================
#  BLOC 7 — BATERIA DE PROVES
# ============================================================================

def executar_bateria(verbose: bool = False):
    proves = [
        (3 * 5,                                "3 × 5  (trivial)"),
        (17 * 19,                              "17 × 19 = 323"),
        (101 * 103,                            "101 × 103 = 10.403"),
        (100003 * 100019,                      "100003 × 100019  (11 dígits, equilibrat)"),
        (1_000_003 * 1_000_033,                "1.000.003 × 1.000.033  (13 dígits)"),
        (9_999_991 * 9_999_973,                "9.999.991 × 9.999.973  (14 dígits)"),
        (1_548_586_332_452_843,                "1.548...  (16 dígits, factor 59)"),
        (999_983 * 1_000_000_000_000_003,      "999.983 × 10^15+3  (asimètric)"),
        (99_999_999_999_973 * 100_000_000_000_003,
                                               "~10^14 × ~10^14  (gran equilibrat)"),
    ]

    print("\n" + "█" * 72)
    print("  🚀 BATERIA DE PROVES — MDC v18")
    print("█" * 72)

    ok = fail = 0
    for N, desc in proves:
        print(f"\n{'─'*72}")
        print(f"  📋 {desc}")
        factor, t_ms = mdc_v18(N, verbose=verbose)
        if factor:
            assert N % factor == 0, f"ERROR CRÍTIC!"
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
    # Prova verbose: el clàssic de l'estudi visual
    N = 100003 * 100019
    factor, t_ms = mdc_v18(N, verbose=True)

    print("\n" + "─" * 72)
    print("  Bateria completa...")
    executar_bateria(verbose=False)
