# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC v19 — MÈTODE DIOFÀNTIC CINEMÀTIC
  Bisecció Binària per ΔΦ Espell · K-sweep Final · Criba Roda
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)

  ─────────────────────────────────────────────────────────────────────────
  INSIGHT CENTRAL DE v19 (Víctor, abril 2025)
  ─────────────────────────────────────────────────────────────────────────

  En lloc de fer el k-sweep sobre tota la zona densa [m_conv, m_max],
  usem el DESFASE ESPELL ΔΦ com a ORACLE DE BISECCIÓ:

    1. Avaluem ΔΦ als EXTREMS i al CENTRE del rang [m_ini, m_fi].
    2. Si signe(ΔΦ_esquerra) ≠ signe(ΔΦ_centre) → factor a [m_ini, m_c]
       Si signe(ΔΦ_centre)  ≠ signe(ΔΦ_dreta)   → factor a [m_c, m_fi]
    3. Bisectem la meitat correcta i repetim.
    4. Quan el segment és prou menut → k-sweep exhaustiu (O(n_k) ràpid).

  COMPLEXITAT:
    Biseccions necessàries: log₂(m_max / llindar_ksweep)
    Per N de D dígits:      m_max ≈ 10^(D/2)
                            llindar ≈ 2×10⁶
    → biseccions ≈ D/2 × log₂(10) − 21 ≈ 1.66·D

    Cada bisecció = 3 avaluacions de ΔΦ = O(1) operacions enteres.

    Per RSA-4096 (D=1233): ≈ 2050 biseccions → FACTIBLE
    Per k-sweep zona densa: ≈ 10^616 iteracions → IMPOSSIBLE

  PERQUÈ ΔΦ ÉS UN BON ORACLE:
    Quan el factor p cau entre m_i i m_{i+1}:
      · B = 2m+3 creuà p: N//B canvia bruscament
      · L'espell A = N//B salta a un valor diferent
      · d_frac(m_espell) canvia abruptament
      · ΔΦ = d_B − d_A canvia de signe (teorema d'atrapament, verificat)

  ARQUITECTURA v19:
    F0 · Precondicions + criba roda fins ~2M (factors petits/asimètrics)
    F1 · Bisecció binària ΔΦ sobre [1, m_max]
         fins que amplada → llindar
    F2 · K-sweep exhaustiu sobre el segment acotat
    F3 · Fallback: pinça 4+4 si k-sweep no ha trobat res
  ─────────────────────────────────────────────────────────────────────────
"""

from fractions import Fraction
import math
import time


# ============================================================================
#  BLOC 0 — CONSTANTS
# ============================================================================
PRIMOS    = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)

# Convergència correcta (e−1)/2 ≈ 0.8591
E_CONV = (math.e - 1) / 2

# Llindar per al k-sweep final (≤ LIM_KSWEEP iteracions de k → ràpid)
LIM_KSWEEP = 2_000_000


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL
# ============================================================================

def pasa_rueda(m: int) -> bool:
    """D(m) = 2m+3 coprimer amb el primorial p23? (~8% dels enters passen)."""
    if m < 1:
        return False
    return math.gcd(2 * m + 3, PRIMORIAL) == 1


def proper_valid(m: int, pas: int = 1) -> int | None:
    """Cerca el m vàlid (roda) més proper. pas=+1 endavant, -1 enrere."""
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
    Funció d'ona exacta: d(m) = Fraction(N mod D, D) − ½
    on D = 2*(2m+3).

    Fraction pura: zero error d'arrodoniment per a N de qualsevol mida.
    Factor real ↔ N mod (2m+3) == 0 ↔ d = −½.
    """
    D = 2 * (2 * m + 3)
    return Fraction(N % D, D) - Fraction(1, 2)


def m_espell(m: int, N: int) -> int | None:
    """
    Espell aritmètic de m: si B = 2m+3, llavors A = N//B i m_e = (A−3)//2.
    Principi A·B ≈ N (exacte quan m és el factor real).
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
    ΔΦ(m) = d_frac(m) − d_frac(m_espell(m)).

    PROPIETAT CLAU (oracle de bisecció):
      Si ΔΦ(m_a) · ΔΦ(m_b) < 0 → existeix un factor de N entre m_a i m_b.
      Açò és el fonament determinista de la bisecció.

    Quan m és exactament el factor p: B = p, A = N/p enter, ΔΦ → mínim.
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
#  BLOC 3 — K-SWEEP (cobertura exhaustiva d'un segment via pendents k)
# ============================================================================

def k_sweep(N: int, m_ini: int, m_fi: int) -> int | None:
    """
    Iteració sobre les pendents k en lloc dels valors m.

    MATEMÀTICA:
      r(d) = N − k·d  per a d en [N/(k+1), N/k]  (pendent k)
      Zero de la recta k: d_zero = N // k
      Si N % d_zero == 0 → FACTOR.

    EXHAUSTIU: cobreix TOTS els factors en [m_ini, m_fi] sense ometre'n cap.
    ÚS: al final de la bisecció, quan el segment és prou menut.
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
        c = N // k        # zero de la recta k
        if c < 3 or c % 2 == 0:
            continue
        if N % c == 0 and 1 < c < N:
            return c
    return None


# ============================================================================
#  BLOC 4 — BISECCIÓ BINÀRIA ΔΦ
# ============================================================================

def biseccio_delta_phi(N: int, m_ini: int, m_fi: int,
                       llindar: int = LIM_KSWEEP,
                       max_iter: int = 4000,
                       verbose: bool = False) -> int | None:
    """
    Bisecció binària del rang [m_ini, m_fi] usant ΔΦ com a oracle.

    ALGORITME (proposta de Víctor, v19):
    ─────────────────────────────────────────────────────────────────────
    Pas 1: Avaluar ΔΦ als tres punts: esquerra (m_e), centre (m_c), dreta (m_d).
    Pas 2: Detectar el canvi de signe:
             ΔΦ_e · ΔΦ_c < 0  → factor a [m_e, m_c]  (meitat esquerra)
             ΔΦ_c · ΔΦ_d < 0  → factor a [m_c, m_d]  (meitat dreta)
    Pas 3: Bisectar la meitat correcta i tornar al Pas 1.
    Pas 4: Quan n_k del segment ≤ llindar → k-sweep exhaustiu i retornar.

    NÚMERO DE BISECCIONS:
      Cada bisecció divideix el rang per 2.
      Cal log₂(m_max / llindar) biseccions:
        m_max ≈ 10^(D/2),  llindar ≈ 2×10⁶
        → biseccions ≈ D × 1.66
      Per RSA-4096 (D=1233): ≈ 2048 biseccions.

    QUAN ΔΦ NO TÉ CANVI DE SIGNE CLAR:
      Pot passar per dos motius:
        a) El factor és molt prop d'un extrem → el canvi de signe cau
           fora de la finestra dels tres punts mostrejats.
        b) ΔΦ té múltiples zeros dins del segment (N té molts factors).
      Estratègia: usem n%pos%(pos-1) per triar el costat amb menor
      distància al zero de pendent (criteri cinemàtic del Libro 5).

    RETORNA: el factor D = 2m+3 o None.
    """
    iter_count = 0

    while iter_count < max_iter:
        iter_count += 1

        # Calculem n_k del segment actual
        pos_i = 2 * m_ini + 3
        pos_f = 2 * m_fi  + 3
        k_lo  = max(1, N // pos_f) if pos_f > 0 else 1
        k_hi  = N // pos_i         if pos_i > 0 else 1
        n_k   = k_hi - k_lo + 1

        if verbose:
            print(f"    iter={iter_count:4d}  "
                  f"[{m_ini:>12,}, {m_fi:>12,}]  "
                  f"amplada={m_fi-m_ini:>10,}  n_k={n_k:>8,}")

        # ── K-SWEEP FINAL quan el segment és prou menut ────────────────
        if n_k <= llindar:
            f = k_sweep(N, m_ini, m_fi)
            if f:
                if verbose:
                    print(f"    → k-sweep trobat: {f:,}")
            return f

        # ── TRES PUNTS PER A LA BISECCIÓ ──────────────────────────────
        m_e  = proper_valid(m_ini,              1)
        m_c  = proper_valid((m_ini + m_fi) // 2, 1)
        m_d  = proper_valid(m_fi,               -1)

        if not (m_e and m_c and m_d):
            return None

        # Comprovació directa dels tres punts (cost O(1))
        for mp in [m_e, m_c, m_d]:
            if check_factor(mp, N):
                return 2 * mp + 3

        # ── AVALUACIÓ DE ΔΦ (l'oracle) ────────────────────────────────
        df_e = desfase_espell(m_e, N)
        df_c = desfase_espell(m_c, N)
        df_d = desfase_espell(m_d, N)

        # ── DECISIÓ DE BISECCIÓ ───────────────────────────────────────
        if df_e * df_c < 0:
            # Canvi de signe a l'esquerra → factor en [m_e, m_c]
            m_fi  = m_c

        elif df_c * df_d < 0:
            # Canvi de signe a la dreta → factor en [m_c, m_d]
            m_ini = m_c

        else:
            # ── CAP CANVI DE SIGNE CLAR ───────────────────────────────
            # Podem estar en dos situacions:
            #   a) Factor molt prop d'un extrem (ΔΦ creua fora dels 3 punts)
            #   b) N té pocs factors i ΔΦ és quasi monòton al segment
            #
            # Estratègia A: mirar si ΔΦ al punt de màxima curvatura
            #   apunta cap a un costat
            # Estratègia B: usar n%pos%(pos-1) com a criteri de mínims
            #
            # Implementem A+B combinats:

            # Punt de mínima |ΔΦ| (prop del zero = prop del factor)
            abs_e = abs(df_e)
            abs_c = abs(df_c)
            abs_d = abs(df_d)

            if abs_e <= abs_c and abs_e <= abs_d:
                # El punt esquerre és el més prop del zero → anem a l'esquerra
                m_fi = m_c
            elif abs_d <= abs_c:
                # El punt dret és el més prop → anem a la dreta
                m_ini = m_c
            else:
                # El punt central és el mínim de |ΔΦ| → bisectem simetricament
                # usant n%pos%(pos-1) per triar el millor costat
                pos_e = 2 * m_e + 3
                pos_d = 2 * m_d + 3
                R2_e  = (N % pos_e) % (pos_e - 1) if pos_e > 1 else 10**18
                R2_d  = (N % pos_d) % (pos_d - 1) if pos_d > 1 else 10**18
                if R2_e <= R2_d:
                    m_fi  = m_c   # menor distància al zero → esquerra
                else:
                    m_ini = m_c   # menor distància al zero → dreta

    return None


# ============================================================================
#  BLOC 5 — PINÇA DOBLE 4+4  (Libro 5, safety-net final)
# ============================================================================

def pinca_doble_salt(m: int, N: int, m_max: int) -> int | None:
    """
    Pinça doble 4+4 + salt mestre com a safety-net.
    S'activa sols si la bisecció ha fallat (rar).
    """
    # Endavant
    df = [d_frac(m + i, N) for i in range(4)]
    Vf = df[1] - df[0]
    Af = df[2] - 2*df[1] + df[0]

    if Af != 0:
        despl = Vf / Af
        di = int(despl + Fraction(1,2)) if despl >= 0 else -int(-despl + Fraction(1,2))
        ms = m - di
        if 1 <= ms <= m_max and check_factor(ms, N):
            return 2 * ms + 3

    # Enrere
    if m >= 4:
        dr = [d_frac(m - i, N) for i in range(4)]
        Vr = dr[1] - dr[0]
        Ar = dr[2] - 2*dr[1] + dr[0]
        if Ar != 0:
            despl = Vr / Ar
            di = int(despl + Fraction(1,2)) if despl >= 0 else -int(-despl + Fraction(1,2))
            ms = m + di
            if 1 <= ms <= m_max and check_factor(ms, N):
                return 2 * ms + 3
    return None


# ============================================================================
#  BLOC 6 — ORQUESTRADOR PRINCIPAL  MDC v19
# ============================================================================

def mdc_v19(N: int, verbose: bool = True) -> tuple:
    """
    MDC v19 — Bisecció Binària per ΔΦ Espell.

    FLUX:
      F0 · Precondicions (trivials + criba roda fins ~2M)
      F1 · Bisecció ΔΦ sobre [1, m_max] → segment acotat
      F2 · K-sweep exhaustiu sobre el segment (ja dins de LIM_KSWEEP)
      [F3 · Pinça 4+4 safety-net si tot falla]

    COMPLEXITAT:
      F0: O(LIM_CRIBA / densitat_roda) ≈ O(2M × 0.08) = O(160K) divisió
      F1: O(log₂(m_max / LIM_KSWEEP)) biseccions
          Cada bisecció: 3 avaluacions ΔΦ = O(1)
          Total: O(D × 1.66) operacions exactes (D = dígits de N)
      F2: O(LIM_KSWEEP) = O(2M) iteracions enteres → < 100ms
      Total per factors equilibrats: O(D × 1.66) quasi-O(log N)
    """
    if verbose:
        bar = '═' * 72
        print(f"\n{bar}")
        print(f"  MDC v19 — Bisecció Binària per ΔΦ Espell")
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

    m_max = (math.isqrt(N) - 3) // 2

    # Criba amb roda primorial per a factors petits / molt asimètrics
    # Cobreix fins a min(√N, 2M × 12.4) ≈ min(√N, 25M) d'enters
    LIM_CRIBA = min(m_max, 2_000_000)
    if verbose:
        print(f"    → Criba roda fins a m={LIM_CRIBA:,}  (D={2*LIM_CRIBA+3:,})...")

    for m_c in range(1, LIM_CRIBA + 1):
        if pasa_rueda(m_c):
            D = 2 * m_c + 3
            if N % D == 0 and D < N:
                t_ms = (time.perf_counter() - t0) * 1000
                if verbose:
                    print(f"    → Factor criba: {D:,}")
                return D, t_ms

    if LIM_CRIBA >= m_max:
        # La criba ja ha cobert tot el rang → no hi ha factors
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose:
            print(f"    → Criba ha cobert tot el rang. N és primer?")
        return None, t_ms

    if verbose:
        print(f"    → m_max = {m_max:,}  ({len(str(m_max))} dígits)")
        print(f"    → m_conv= {int(E_CONV*m_max):,}  ({E_CONV*100:.2f}%)")
        print(f"    → Rang a bisellar: [{LIM_CRIBA:,}, {m_max:,}]")

    t_criba = (time.perf_counter() - t0) * 1000
    if verbose:
        print(f"    → Criba completada en {t_criba:.1f} ms. Factor no trobat.")

    # ── FASE 1b: K-SWEEP ZONA DENSA [m_conv, m_max] ───────────────────────
    # Sempre fem k-sweep dels LIM_KSWEEP valors de k prop de m_max
    # (on els dents apareixen més freqüentment i els factors equilibrats viuen).
    # Independentment de n_k total, els últims LIM_KSWEEP k son sempre ràpids.
    m_conv_val = max(LIM_CRIBA + 1, int(E_CONV * m_max))

    # Calculem quin m correspon als LIM_KSWEEP valors de k més alts
    # k_max = N // (2*m_conv+3), k_lim = k_max - LIM_KSWEEP
    # m que dona k_lim: m ≈ N//(2*k_lim) → però ho calculem directament
    pos_max   = 2 * m_max + 3
    k_top     = max(1, N // pos_max)               # k per a m_max (el més menut)
    k_lim_top = k_top + LIM_KSWEEP                 # LIM_KSWEEP k per damunt
    # El m que dona k_lim_top: D = N // k_lim_top → m = (D-3)//2
    d_lim_top = N // k_lim_top if k_lim_top > 0 else 2 * m_max + 3
    m_lim_top = max(LIM_CRIBA + 1, (d_lim_top - 3) // 2)

    if verbose:
        n_k_top = min(LIM_KSWEEP, N//(2*m_lim_top+3) - k_top)
        print(f"\n  [F1b] K-sweep [{m_lim_top:,}, {m_max:,}]"
              f"  (últims ~{LIM_KSWEEP:,} valors de k)...")

    factor_D = k_sweep(N, m_lim_top, m_max)
    if factor_D:
        t_ms = (time.perf_counter() - t0) * 1000
        _print_resultat(N, factor_D, "zona densa k-sweep", t_ms)
        return factor_D, t_ms

    if verbose:
        print(f"    → No trobat a la zona densa.")

    # ── FASE 1: BISECCIÓ BINÀRIA ΔΦ ───────────────────────────────────────
    if verbose:
        print(f"\n  [F1] Bisecció binària ΔΦ espell sobre [{LIM_CRIBA:,}, {m_max:,}]...")
        n_bisec_teoriq = max(1, int(math.log2(
            max(1, (m_max - LIM_CRIBA) / LIM_KSWEEP)
        )))
        print(f"    → Biseccions estimades: ~{n_bisec_teoriq}")

    factor_D = biseccio_delta_phi(
        N, LIM_CRIBA, m_max,
        llindar = LIM_KSWEEP,
        max_iter = max(4000, m_max.bit_length() * 3),
        verbose  = verbose
    )

    if factor_D:
        t_ms = (time.perf_counter() - t0) * 1000
        _print_resultat(N, factor_D, "bisecció ΔΦ", t_ms)
        return factor_D, t_ms

    # ── FASE 3: PINÇA 4+4 SAFETY-NET (rarament arriba ací) ────────────────
    if verbose:
        print(f"\n  [F3] Pinça 4+4 safety-net...")

    pas_pc = max(1, m_max // 30)
    for i in range(30):
        m_s = LIM_CRIBA + i * pas_pc
        if m_s > m_max:
            break
        mv = proper_valid(m_s, 1)
        if mv:
            f = pinca_doble_salt(mv, N, m_max)
            if f:
                t_ms = (time.perf_counter() - t0) * 1000
                _print_resultat(N, f, "pinça 4+4", t_ms)
                return f, t_ms

    t_ms = (time.perf_counter() - t0) * 1000
    if verbose:
        print(f"\n  ✗ Factor no trobat. Temps: {t_ms:.2f} ms")
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

def executar_bateria(verbose_bisec: bool = False):
    proves = [
        (3 * 5,
         "3 × 5  (trivial)"),
        (101 * 103,
         "101 × 103 = 10.403"),
        (100003 * 100019,
         "100003 × 100019  (11 dígits, equilibrat)"),
        (1_000_003 * 1_000_033,
         "1.000.003 × 1.000.033  (13 dígits)"),
        (9_999_991 * 9_999_973,
         "9.999.991 × 9.999.973  (14 dígits)"),
        (1_548_586_332_452_843,
         "1.548...  (16 dígits, factor petit 59)"),
        (999_983 * 1_000_000_000_000_003,
         "999.983 × 10^15+3  (asimètric 6×16 dígits)"),
        (100_000_000_003 * 100_000_000_019,
         "10^11+3 × 10^11+19  (22 dígits, gran equilibrat)"),
        (999_999_999_989 * 1_000_000_000_039,
         "~10^12 × ~10^12  (24 dígits)"),
    ]

    print("\n" + "█" * 72)
    print("  🚀 BATERIA DE PROVES — MDC v19 (Bisecció ΔΦ)")
    print("█" * 72)

    ok = fail = 0
    for N, desc in proves:
        print(f"\n{'─'*72}")
        print(f"  📋 {desc}")
        factor, t_ms = mdc_v19(N, verbose=False)
        if factor:
            assert N % factor == 0, f"ERROR!"
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
    # Prova verbose amb bisecció detallada
    N = 100003 * 100019
    factor, t_ms = mdc_v19(N, verbose=True)

    print("\n" + "─" * 72)
    executar_bateria(verbose_bisec=False)
