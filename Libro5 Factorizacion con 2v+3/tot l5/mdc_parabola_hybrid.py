# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC PARABOLA HYBRID
  Algoritme Híbrid de Factorització per Semiprims
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)
  VERSIÓ: Hybrid-1.0
  COL·LABORACIÓ D'ESCRIPTURA: Claude (Anthropic)

  ─────────────────────────────────────────────────────────────────────────
  TRES MÈTODES, UNA ARQUITECTURA
  ─────────────────────────────────────────────────────────────────────────

  MÈTODE A — K-SWEEP PER PENDENTS (de v17-v23)
  ─────────────────────────────────────────────
  Iteració sobre els cocients k = N//D en lloc dels valors m.
  Per cada k en [k_lo, k_hi]:
      candidat = N // k
      Si N % candidat == 0 → FACTOR TROBAT
  Exhaustiu i determinista. O(n_k) on n_k és el nombre de valors de k
  en el rang de cerca. Per a factors molt equilibrats (p ≈ q ≈ √N),
  n_k ≈ LIM_KSWEEP → quasi instantani.

  MÈTODE B — FUNCIÓ D'ONA MDC v15 (cinemàtica)
  ─────────────────────────────────────────────
  Funció d'ona: d(m) = Fraction(N mod D, D) − ½, on D = 2*(2m+3)
  Usada per detectar "zones calentes": valors de m on |d(m)| és petit
  (prop del zero de fase → prop d'un factor).
  La cinemàtica (V, A) permet saltar directament al vèrtex de la paràbola
  local sense iterar tot el rang.

  MÈTODE C — DISCRIMINANT S±k (Libro 5)
  ─────────────────────────────────────
  Si N = (2v+3)(2b+3), llavors:
      S = v + b,   k = v - b
      Δ(S) = S² + 6S − (N − 9) = k²
  Cerca el S que fa Δ(S) un quadrat perfecte.
  Un cop detectat S per MDC (com a pista), el discriminant confirma
  algebraicament els dos factors en O(1) operacions.

  ─────────────────────────────────────────────────────────────────────────
  PER QUÈ COMBINAR-LOS?
  ─────────────────────────────────────────────────────────────────────────

  · K-sweep és el mètode més ràpid per a factors prop de √N (equilibrats).
    Però per a factors lluny de √N necessita massa iteracions.

  · MDC v15 detecta zones calentes en O(n_mostra) avaluacions —
    cadascuna O(1). Però la cinemàtica no sempre salta exactament al factor.

  · Discriminant S±k és EXACTE i O(1) un cop coneixes S aproximadament.
    Però cercar S en [√N, N/6] és massa lent sense una pista.

  FLUX HÍBRID:
    1. K-sweep [m_lim, m_max]: cobreix factors equilibrats → ràpid
    2. MDC detector: per cada segment, genera candidats m (zones calentes)
    3. Discriminant: per cada candidat m, calcula S = m_a_S(m) i comprova
       si Δ(S) és quadrat perfecte → factoritza exactament
    4. Fallback: criba directa per N petits o quan tot falla

  ─────────────────────────────────────────────────────────────────────────
  ARXIUS DE REFERÈNCIA DEL CORPUS:
    mdc_v15.py  → d_frac, check_factor, pasa_rueda, mesurar_cinematica
    mdc_v20.py  → k_sweep per pendents
    mdc_fusio_v3.py → Teorema 1, Teorema 2, salt parabòlic
    mdc_factorizador.py → predecir(), verificar_zona()
    codigos.txt → Δ(S)=S²+6S−(N−9)=k², S±k de Fermat
  ─────────────────────────────────────────────────────────────────────────
"""

from fractions import Fraction
import math
import time


# ============================================================================
#  BLOC 0 — CONSTANTS GLOBALS
# ============================================================================

# Primers petits per a la criba inicial i la roda primorial
PRIMOS    = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)   # 223.092.870

# Nombre màxim de valors de k en el k-sweep inicial (Mètode A)
# 2M valors ≈ 400ms en Python pur → acceptable com a primera fase
LIM_KSWEEP = 2_000_000

# Umbral de "zona calenta" cinemàtica (Mètode B)
# Si |d(m)| < UMBRAL_FASE considerem m com a candidat fort
UMBRAL_FASE = Fraction(1, 8)   # = 0.125 (ajustable)

# Radi de cerca al voltant de la projecció cinemàtica (Mètode B)
RADI_SALT = 3


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL (filtre topològic)
# ============================================================================

def pasa_rueda(m: int) -> bool:
    """
    Comprova si D(m) = 2m+3 és coprimer amb el primorial p23.

    Si mcd(2m+3, PRIMORIAL) > 1 → 2m+3 té un factor petit → no pot ser
    un factor no trivial de N (assumint que N ja no és divisible pels
    primers de la roda).

    Densitat de valors que passen: φ(PRIMORIAL)/PRIMORIAL ≈ 8.07%
    → la roda descarta ~92% dels enters sense perdre cap factor.
    """
    if m < 1:
        return False
    return math.gcd(2 * m + 3, PRIMORIAL) == 1


def proper_valid(m: int, pas: int = 1) -> int | None:
    """
    Cerca el m vàlid (que passa la roda) més proper.
    pas = +1 endavant, -1 enrere.
    """
    for _ in range(500):
        if m < 1:
            return None
        if pasa_rueda(m):
            return m
        m += pas
    return None


# ============================================================================
#  BLOC 2 — FUNCIÓ D'ONA i VERIFICACIÓ (Mètode B, MDC v15)
# ============================================================================

def d_frac(m: int, N: int) -> Fraction:
    """
    Funció d'ona del MDC: d(m) = Fraction(N mod D, D) − ½
    on D = 2*(2m+3).

    · d(m) → 0    quan 2m+3 ≈ factor de N (zona calenta)
    · N % (2m+3) == 0 exactament → factor trobat

    Fraction exacta: evita errors d'arrodoniment per a N de molts dígits.
    """
    D = 2 * (2 * m + 3)
    R = N % D
    return Fraction(R, D) - Fraction(1, 2)


def check_factor(m: int, N: int) -> bool:
    """Verificació directa: (2m+3) divideix N?"""
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


def mesurar_cinematica(m: int, N: int) -> tuple:
    """
    Mesura cinemàtica de la funció d'ona en el punt m.

    Retorna (d0, V, A) amb Fraction exacta:
      d0 = posició actual
      V  = d(m+1) − d(m)        velocitat (1a diferència)
      A  = d(m+2) − 2d(m+1) + d(m)  acceleració (2a diferència)

    El vèrtex de la paràbola local és a: m* = m − V/A
    """
    d0 = d_frac(m,     N)
    d1 = d_frac(m + 1, N)
    d2 = d_frac(m + 2, N)
    V  = d1 - d0
    A  = d2 - 2 * d1 + d0
    return d0, V, A


# ============================================================================
#  BLOC 3 — DISCRIMINANT S±k (Mètode C, Libro 5)
# ============================================================================

def m_a_S(m: int, N: int) -> int:
    """
    Converteix un candidat m a la variable S de Fermat.

    Si m és el factor p = 2m+3, aleshores el cofactor és q = N // p.
    La suma d'índexs és:
        v = (p − 3) // 2 = m
        b = (q − 3) // 2
        S = v + b

    RETORNA: S = m + (N // (2m+3) − 3) // 2
    """
    p = 2 * m + 3
    if p <= 1 or p >= N:
        return -1
    q = N // p
    if q < 3:
        return -1
    b = (q - 3) // 2
    return m + b   # S = v + b


def discriminant_exacte(S: int, N: int) -> tuple:
    """
    Comprova si Δ(S) = S² + 6S − (N − 9) és un quadrat perfecte.

    MATEMÀTICA (Libro 5, Teorema 1):
        N = (2v+3)(2b+3) = 4vb + 6(v+b) + 9
        S = v + b,  k = v − b
        Δ(S) = S² + 6S − (N−9) = (v−b)² = k²

    Si Δ(S) ≥ 0 i és quadrat perfecte k²:
        v = (S + k) // 2
        b = (S − k) // 2
        p = 2v + 3,  q = 2b + 3
        p × q = N  (Teorema 1)

    RETORNA: (p, q) si factoritza, o (None, None) en cas contrari.
    """
    if S < 0:
        return None, None

    delta = S * S + 6 * S - (N - 9)

    if delta < 0:
        return None, None

    k = math.isqrt(delta)
    if k * k != delta:
        return None, None   # Δ no és quadrat perfecte

    # Extraiem v i b
    Sp_k = S + k
    Sm_k = S - k

    if Sp_k % 2 != 0 or Sm_k % 2 != 0:
        return None, None   # v i b han de ser enters

    v = Sp_k // 2
    b = Sm_k // 2

    if b < 0:
        return None, None   # factor trivial (B = 2b+3 < 3)

    p = 2 * v + 3
    q = 2 * b + 3

    if p > 1 and q > 1 and p != N and q != N and p * q == N:
        return min(p, q), max(p, q)

    return None, None


def cerca_discriminant_zona(m_centre: int, N: int,
                             radi_S: int = 50) -> tuple:
    """
    Cerca el quadrat perfecte en Δ(S) per a S properes a m_centre.

    Donada una zona calenta m detectada per MDC (Mètode B),
    calculem S_0 = m_a_S(m) i explorem S ∈ [S_0 − radi_S, S_0 + radi_S].

    Açò combina la potència de detecció de zones del MDC amb la precisió
    algebraica del discriminant.

    RETORNA: (p, q) si factoritza, o (None, None).
    """
    S_base = m_a_S(m_centre, N)
    if S_base < 0:
        return None, None

    # Explorar S properes a S_base (en ambdues direccions)
    for delta_S in range(radi_S + 1):
        for signe in [0, 1] if delta_S == 0 else [1, -1]:
            S = S_base + signe * delta_S
            if S < 0:
                continue
            p, q = discriminant_exacte(S, N)
            if p is not None:
                return p, q

    return None, None


def cerca_discriminant_rang(N: int, S_ini: int, S_fi: int) -> tuple:
    """
    Cerca exhaustiva del quadrat perfecte en Δ(S) per a S en [S_ini, S_fi].

    Usat com a fallback quan el k-sweep i el MDC no han trobat res.
    Equivalent a la "trayectoria S" del codigos.txt (Mati/Gemini).

    MATEMÀTICA:
      S creix de 1 en 1 (o de 2 en 2 si la paritat ho permet).
      Per a cada S: Δ(S) = S² + 6S − (N−9). Si és k² → factoritza.

    NOTA: Per a N gran, este rang és massa llarg per a recórrer complet.
    S'usa sols en el rang estret [√N − marge, √N + marge].

    RETORNA: (p, q) o (None, None).
    """
    for S in range(S_ini, S_fi + 1):
        p, q = discriminant_exacte(S, N)
        if p is not None:
            return p, q
    return None, None


# ============================================================================
#  BLOC 4 — K-SWEEP PER PENDENTS (Mètode A, de v17-v23)
# ============================================================================

def k_sweep(N: int, m_ini: int, m_fi: int) -> int | None:
    """
    Iteració sobre les pendents k = N//D en lloc dels valors m.

    MATEMÀTICA (pendents lineals del Libro 5):
      La funció d(m) és lineal per trams: per a m ∈ [N/(k+1), N/k],
      el cocient k = N//(2m+3) és constant. El zero de cada pendent
      (on d(m)=0.5, és a dir, on hi podria haver un factor) és:
          candidat = N // k
      Si N % candidat == 0 → FACTOR TROBAT.

    PER QUÈ ÉS EXHAUSTIU:
      Per cada m ∈ [m_ini, m_fi], la seua pendent k = N//(2m+3) cau en
      [k_lo, k_hi]. Iterant sobre totes les k i comprovant N // k,
      cobrim tots els factors sense ometre'n cap.

    AVANTATGE per a factors equilibrats:
      Quan p ≈ q ≈ √N: k ≈ √N i el rang [k_lo, k_hi] és molt estret
      → l'algoritme és quasi instantani.

    RETORNA: el factor (enter), o None si no n'hi ha en [m_ini, m_fi].
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
        candidat = N // k   # zero de la recta k
        if candidat < 3 or candidat % 2 == 0:
            continue        # (2m+3) és sempre senar
        if N % candidat == 0 and 1 < candidat < N:
            return candidat

    return None


# ============================================================================
#  BLOC 5 — DETECTOR MDC (zones calentes cinemàtiques)
# ============================================================================

def detectar_zones_calentes(N: int, m_ini: int, m_fi: int,
                             n_mostra: int = 60,
                             verbose: bool = False) -> list[int]:
    """
    Detecta "zones calentes" en [m_ini, m_fi]: valors de m on la funció
    d'ona d(m) és propera a zero (indicadors cinemàtics de factors propers).

    ESTRATÈGIA:
      1. Mostrar n_mostra punts equidistants en [m_ini, m_fi].
      2. Per cada punt m que passa la roda primorial:
         a. Calcular d0 = d_frac(m, N) (Fraction exacta).
         b. Si |d0| < UMBRAL_FASE → zona calenta → afegir a la llista.
         c. Mesura cinemàtica: si A ≠ 0, projectar salt m* = m − V/A.
            Si m* ∈ [m_ini, m_fi] → afegir m* com a zona calenta secundària.
      3. Ordenar per |d(m)| ascendent (zones més calentes primer).

    RETORNA: llista de m candidats (zones calentes), ordenada per calor.
    """
    if m_ini >= m_fi or m_ini < 1:
        return []

    amplada = m_fi - m_ini
    pas = max(1, amplada // n_mostra)

    zonas = []   # (|d0|, m)

    m = m_ini
    while m <= m_fi:
        if pasa_rueda(m):

            # ── Verificació directa (sempre) ──────────────────────────────
            if check_factor(m, N):
                zonas.insert(0, (Fraction(0), m))   # zona absolutament calenta
                if verbose:
                    print(f"    [Directe] Factor immediat a m={m:,}")
                return [m]

            # ── Mesura d(m) exacta ─────────────────────────────────────────
            d0, V, A = mesurar_cinematica(m, N)

            if abs(d0) < UMBRAL_FASE:
                zonas.append((abs(d0), m))
                if verbose:
                    print(f"    [Zona calenta] m={m:,}  |d|={float(abs(d0)):.4f}")

            # ── Salt cinemàtic: projectar vèrtex de la paràbola ───────────
            if A != Fraction(0):
                despl = V / A
                if despl >= 0:
                    despl_int = int(despl + Fraction(1, 2))
                else:
                    despl_int = -int(-despl + Fraction(1, 2))
                m_proj = m - despl_int
                if m_ini <= m_proj <= m_fi and pasa_rueda(m_proj):
                    dp = abs(d_frac(m_proj, N))
                    if dp < UMBRAL_FASE * 2:   # llindar més permissiu per a projeccions
                        zonas.append((dp, m_proj))
                        if verbose:
                            print(f"    [Salt cinem.] Projecció m_proj={m_proj:,}"
                                  f"  |d|={float(dp):.4f}")

        m += pas

    # Ordenar per |d| ascendent, eliminar duplicats
    vistos = set()
    resultat = []
    for _, m_cand in sorted(zonas, key=lambda x: x[0]):
        if m_cand not in vistos:
            vistos.add(m_cand)
            resultat.append(m_cand)

    return resultat


# ============================================================================
#  BLOC 6 — ORQUESTRADOR HÍBRID PRINCIPAL
# ============================================================================

def mdc_parabola_hybrid(N: int,
                        n_mostra: int = 60,
                        verbose: bool = True) -> tuple:
    """
    MDC PARABOLA HYBRID — Algoritme híbrid de factorització.

    Combina tres mètodes complementaris:

    ┌─────────────────────────────────────────────────────────────────────┐
    │  FASE 0: Precondicions (trivials, quadrat perfecte)                 │
    │  FASE 1: K-sweep [m_lim, m_max]     ← Mètode A (ràpid per p≈q)    │
    │  FASE 2: MDC + Discriminant S±k     ← Mètodes B+C (precís)        │
    │    · MDC detecta zones calentes cinemàticament                     │
    │    · Per cada zona: discriminant confirma algebraicament            │
    │  FASE 3: Discriminant rang estret   ← Mètode C pur (exacte)       │
    │  FASE 4: Fallback directe (N petits)← Criba de seguretat           │
    └─────────────────────────────────────────────────────────────────────┘

    ARGUMENTS:
      N        : enter a factoritzar (senar, ≥ 9)
      n_mostra : punts de mostra per segment en el detector MDC
      verbose  : si True, imprimeix el flux detallat

    RETORNA: ((p, q), temps_ms) on p ≤ q o ((None, None), temps_ms)
    """

    def log(msg):
        if verbose:
            print(msg)

    def retornar(p, q, etapa):
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose:
            bar = '═' * 70
            print(f"\n{bar}")
            print(f"  🎯 FACTOR TROBAT! [{etapa}]")
            print(f"     Factor p : {p:,}")
            print(f"     Factor q : {q:,}")
            print(f"     Temps    : {t_ms:.3f} ms")
            print(f"{bar}")
        return (min(p, q), max(p, q)), t_ms

    # ── CAPÇALERA ────────────────────────────────────────────────────────────
    if verbose:
        bar = '═' * 70
        print(f"\n{bar}")
        print(f"  MDC PARABOLA HYBRID  (MDC + Cinemàtica + Discriminant S±k)")
        print(f"  N = {N:,}  ({len(str(N))} dígits)")
        print(f"{bar}")

    t0 = time.perf_counter()

    # ────────────────────────────────────────────────────────────────────────
    # FASE 0 — PRECONDICIONS
    # ────────────────────────────────────────────────────────────────────────
    log("\n  [FASE 0] Precondicions...")

    # 0a. Factors trivials (primers de la roda)
    for p in PRIMOS:
        if N % p == 0 and p < N:
            log(f"    → Factor trivial: {p}")
            return retornar(p, N // p, "F0 trivial")

    # 0b. Quadrat perfecte: N = r²
    r = math.isqrt(N)
    if r * r == N:
        log(f"    → Quadrat perfecte: √N = {r}")
        return retornar(r, r, "F0 quadrat")

    # 0c. Criba ràpida amb roda primorial per a factors petits (<4M)
    # Cobreix factors molt asimètrics (p << q) sense costar massa temps
    m_max_tmp = (math.isqrt(N) - 3) // 2
    LIM_CRIBA = min(m_max_tmp + 1, 2_000_000)
    log(f"    → Criba roda fins a D={2*LIM_CRIBA+3:,}...")
    for m_c in range(1, LIM_CRIBA + 1):
        if pasa_rueda(m_c):
            D = 2 * m_c + 3
            if N % D == 0 and D < N:
                log(f"    → Factor criba: {D:,}")
                return retornar(D, N // D, "F0 criba")

    # 0d. Paràmetres bàsics
    m_max = m_max_tmp
    log(f"    → m_max = {m_max:,}  (√N ≈ {2*m_max+3:,})")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 1 — K-SWEEP (Mètode A): factors molt equilibrats prop de √N
    # ────────────────────────────────────────────────────────────────────────
    log(f"\n  [FASE 1] K-sweep Mètode A: {LIM_KSWEEP:,} valors de k prop de √N...")

    # Calculem el rang de m que correspon a LIM_KSWEEP valors de k
    # k_top = N // (2*m_max + 3)  (la k mínima, per al m màxim)
    # Cobrim els LIM_KSWEEP valors de k immediats per sobre de k_top
    k_top = max(1, N // (2 * m_max + 3))
    k_lim = k_top + LIM_KSWEEP
    d_lim = N // k_lim if k_lim > 0 else 2 * m_max + 3
    m_lim = max(1, (d_lim - 3) // 2)

    log(f"    → k_top = {k_top:,}  |  m_lim = {m_lim:,}")

    factor_A = k_sweep(N, m_lim, m_max)

    if factor_A:
        return retornar(factor_A, N // factor_A, "F1 k-sweep")

    log(f"    → Cap factor en la zona k-sweep. Continuant...")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 2 — MDC DETECTOR + DISCRIMINANT (Mètodes B + C)
    # ────────────────────────────────────────────────────────────────────────
    log(f"\n  [FASE 2] MDC detector (zones calentes) + Discriminant S±k...")

    # Segmentació armònica del rang [1, m_max] (com en v15)
    # Prioritzem la zona alta (equilibrada) on els factors solen ser
    segments = _generar_segments(m_max, m_lim)

    log(f"    → {len(segments)} segments MDC a escanejar")

    for seg_idx, (m_ini_s, m_fi_s) in enumerate(segments):

        if verbose:
            print(f"\n    Segment {seg_idx+1}/{len(segments)}: "
                  f"[{m_ini_s:,}, {m_fi_s:,}]  "
                  f"(amplada={m_fi_s-m_ini_s:,})", end="  ")

        # ── MÈTODE B: Detector MDC (zones calentes) ───────────────────────
        zones = detectar_zones_calentes(N, m_ini_s, m_fi_s,
                                        n_mostra=n_mostra,
                                        verbose=False)

        if not zones:
            log("→ cap zona calenta")
            continue

        log(f"→ {len(zones)} zones calentes")

        # ── MÈTODES B+C: Per cada zona calenta, aplicar discriminant ─────
        for m_zona in zones:

            # B1. Verificació directa del MDC
            if check_factor(m_zona, N):
                return retornar(2*m_zona+3, N//(2*m_zona+3), "F2 MDC directe")

            # B2. Cerca en el veïnat ±RADI_SALT
            for dm in range(-RADI_SALT, RADI_SALT + 1):
                m_t = m_zona + dm
                if m_t >= 1 and pasa_rueda(m_t) and check_factor(m_t, N):
                    return retornar(2*m_t+3, N//(2*m_t+3), "F2 MDC veïnat")

            # C. Discriminant S±k (Mètode C) des de la zona calenta com a pista
            # La zona m_zona ens dona una estimació de S = v + b
            p_d, q_d = cerca_discriminant_zona(m_zona, N, radi_S=40)
            if p_d is not None:
                return retornar(p_d, q_d, "F2 MDC+Discriminant")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 3 — DISCRIMINANT RANG ESTRET (Mètode C pur)
    # ────────────────────────────────────────────────────────────────────────
    log(f"\n  [FASE 3] Discriminant S±k rang estret (zona equilibrada)...")

    # El discriminant és més eficient per a factors equilibrats (S ≈ 2√N − 3)
    # Explorar en un rang estret de ±1000 al voltant de S_init = ⌈2√N⌉ − 3
    S_init = max(0, 2 * math.isqrt(N) - 3)
    S_rang = min(1000, S_init // 2)   # rang de cerca (ajustable)

    log(f"    → Rang S: [{S_init - S_rang:,}, {S_init + S_rang:,}]")

    p_d3, q_d3 = cerca_discriminant_rang(N, S_init - S_rang, S_init + S_rang)
    if p_d3 is not None:
        return retornar(p_d3, q_d3, "F3 Discriminant pur")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 4 — FALLBACK DIRECTE (criba de seguretat per N menuts)
    # ────────────────────────────────────────────────────────────────────────
    log(f"\n  [FASE 4] Fallback: criba directa L1 (N ≤ LIM en precaucions)...")

    # Únicament recorrem tot el rang si N és prou menut
    # Per a N grans (> 10^14) açò seria massa lent → informem i sortim
    LIM_FALLBACK = 10_000_000   # ≈ 10M candidats L1 → < 1s en Python

    if m_max <= LIM_FALLBACK:
        log(f"    → m_max={m_max:,} ≤ LIM_FALLBACK → criba completa")
        for m_f in range(1, m_max + 1):
            if pasa_rueda(m_f) and check_factor(m_f, N):
                return retornar(2*m_f+3, N//(2*m_f+3), "F4 criba completa")
    else:
        log(f"    → m_max={m_max:,} > LIM_FALLBACK → N podria ser primer")

    # Cap fase ha trobat un factor → N és primer (o cal augmentar n_mostra)
    t_ms = (time.perf_counter() - t0) * 1000
    log(f"\n  ✗ Cap factor trobat. N és primer o cal augmentar n_mostra.")
    log(f"    Temps total: {t_ms:.3f} ms")
    return (None, None), t_ms


# ============================================================================
#  AUXILIAR: GENERACIÓ DE SEGMENTS
# ============================================================================

def _generar_segments(m_max: int, m_lim: int) -> list[tuple]:
    """
    Genera segments harmònics en [1, m_lim] (la zona NO coberta pel k-sweep).

    Segmentació: [m_max // (i+1), m_max // i] per i = 1, 2, 3, ...
    Filtra els segments que ja estan coberts per k-sweep (per sobre de m_lim).

    Prioritat: segments superiors (factors equilibrats) primer.
    """
    segments = []
    i = 1
    while True:
        m_fi_s  = min(m_lim - 1, m_max // i)        # limitar per sota de m_lim
        m_ini_s = m_max // (i + 1)
        if m_fi_s < 1:
            break
        if m_ini_s < 1:
            m_ini_s = 1
        if m_fi_s > m_ini_s:
            segments.append((m_ini_s, m_fi_s))
        i += 1
        if i > 150 or m_max // i < 1:
            break

    # Ordenar per m_fi descendent (segments superiors primer → factors equilibrats)
    segments.sort(key=lambda s: s[1], reverse=True)
    return segments


# ============================================================================
#  BLOC 7 — BATERIA DE PROVES
# ============================================================================

def executar_bateria(verbose_seg: bool = False):
    """
    Bateria de proves per verificar el híbrid.

    Per cada cas:
      · Executa mdc_parabola_hybrid()
      · Verifica que p × q = N
      · Registra el temps i la fase on s'ha trobat el factor
    """
    proves = [
        # (N,  descripció)
        (3 * 5,
         "Trivial: 3×5=15"),
        (17 * 19,
         "Petit: 17×19=323"),
        (101 * 103,
         "Clàssic MDC: 101×103=10.403"),
        (100003 * 100019,
         "Target original: 100.003×100.019 (10 dígits)"),  # ← el del prompt
        (1_000_003 * 1_000_033,
         "Equilibrat 7d×7d"),
        (9_999_991 * 9_999_973,
         "Equilibrat 8d×8d"),
        (1_548_586_332_452_843,
         "Factor petit 59 (16 dígits)"),
        (999_983 * 1_000_000_000_000_003,
         "Asimètric 6d×16d"),
        (100_000_000_003 * 100_000_000_019,
         "Equilibrat 12d×12d"),
    ]

    print("\n" + "█" * 70)
    print("  🚀 BATERIA DE PROVES — MDC PARABOLA HYBRID")
    print("  [MDC + Cinemàtica + Discriminant S±k]")
    print("█" * 70)

    ok_total = fail_total = 0

    for N, desc in proves:
        print(f"\n{'─'*70}")
        print(f"  📋 {desc}")
        print(f"  N = {N:,}  ({len(str(N))} dígits)")

        (p, q), t_ms = mdc_parabola_hybrid(N, n_mostra=60, verbose=verbose_seg)

        if p is not None:
            assert p * q == N, f"ERROR! {p}×{q} ≠ {N}"
            print(f"  ✅  {p:,} × {q:,}   ({t_ms:.2f} ms)")
            ok_total += 1
        else:
            print(f"  ⚠️   No trobat  ({t_ms:.2f} ms)")
            fail_total += 1

    print(f"\n{'█'*70}")
    print(f"  Encerts: {ok_total}/{ok_total+fail_total}")
    print("█" * 70)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

if __name__ == "__main__":

    # ── Prova individual verbose (el cas del prompt) ───────────────────────
    N_PROVA = 100003 * 100019
    print(f"Prova individual: N = {N_PROVA:,} = 100.003 × 100.019")
    (p, q), t_ms = mdc_parabola_hybrid(N_PROVA, n_mostra=60, verbose=True)

    if p:
        print(f"\nResultat: {p:,} × {q:,}  ({t_ms:.3f} ms)")
        assert p * q == N_PROVA
    else:
        print("\nNo trobat (N podria ser primer)")

    # ── Bateria completa (sense verbose detallat) ──────────────────────────
    print()
    executar_bateria(verbose_seg=False)
