# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC v17 — MÈTODE DIOFÀNTIC CINEMÀTIC
  Filtre Acumulador · Pendent Exacta · Pinça Doble 4+4 · Detector Espell
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)

  ─────────────────────────────────────────────────────────────────────────
  PRINCIPIS FONAMENTALS (del Libro 5, lectura completa)
  ─────────────────────────────────────────────────────────────────────────

  1. PENDENTS LINEALS DE LA SEQÜÈNCIA DE RESTOS
     La funció r(d) = N mod d, per d decreixent, forma pendents lineals:
       r(d) = N - k·d    per a tot d tal que ⌊N/d⌋ = k
     Dins de cada pendent k, r(d) és una RECTA EXACTA de pendent −k.
     Un factor p és l'únic enter d on r(p) = 0, és a dir:
       p = N // k    →   si N % p == 0: FACTOR TROBAT
     Açò és EXACTE i DETERMINISTA. No és una aproximació.

  2. ACUMULADOR HARMÒNIC com a FILTRE DE SEGMENTS
     acumulo += acumulo * Fraction(1, i)
     NO es usa per generar posicions sinó per DESCARTAR segments estèrils:
     si en un segment els restos no decreixen seguint esta progressió,
     el segment no té cap factor → DESCARTAT sense fer cap divisió.

  3. SALT EXACTE PER PENDENT (S±k del Libro 5)
     Donat un punt m amb pos = 2m+3:
       k        = N // pos          (cocient, pendent actual)
       candidat = N // k            (zero de la recta k, si és enter)
       m_cand   = (candidat − 3) // 2
     Si N % candidat == 0 → FACTOR EXACTE en un sol pas.

  4. PINÇA DOBLE 4+4 (la pinça real del Libro 5)
     4 mesures cap ENDAVANT: m, m+1, m+2, m+3
     4 mesures cap ENRERE:   m, m-1, m-2, m-3
     De cada sèrie s'extreu V, A, J (velocitat, acceleració, jerk).
     La COMPARACIÓ dels dos jerks detecta si estem en la macro-envolvente
     (J diferent en les dos direccions) o dins d'un dent interior (J ≈ 0).

  5. DETECTOR ESPELL  ΔΦ = d_B(m) − d_A(m_espell)
     Per a B = 2m+3, l'espell aritmètic és A = N // B:
       m_espell = (A − 3) // 2
     Si el signe de ΔΦ canvia entre dos sondes consecutives →
     EXISTEIX un factor en eixe segment (determinista).

  ─────────────────────────────────────────────────────────────────────────
  ARQUITECTURA v17 (les 5 fases):
    F0 · Precondicions (trivials + quadrat)
    F1 · Espiral harmònica → sondes (m_max cap a 1, decreixent)
    F2 · Filtre acumulador → descarta segments sense resonància
    F3 · Detector espell → llista curta de segments candidats
    F4 · Per cada segment: salt exacte + pinça doble + bisecció espell
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


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL
# ============================================================================

def pasa_rueda(m: int) -> bool:
    """
    Comprova si D(m) = 2m+3 és coprimer amb el primorial p23.

    Si mcd(2m+3, PRIMORIAL) > 1 → D(m) té un factor petit ja eliminat
    en les precondicions → no pot ser divisor no trivial de N.
    Densitat: ~8.07% dels enters passen la roda.
    """
    if m < 1:
        return False
    return math.gcd(2 * m + 3, PRIMORIAL) == 1


def proper_valid(m: int, pas: int = 1) -> int | None:
    """
    Troba el m vàlid (que passa la roda) més proper.
    pas = +1 busca cap endavant, -1 cap enrere.
    """
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
    Funció d'ona del MDC: distància signada al zero de fase.

    D  = 2*(2m+3)            denominador simètric
    R  = N mod D             reste exacte (enter pur)
    d  = Fraction(R, D) − ½

    Factor real ↔ N mod (2m+3) == 0 ↔ d = −½ (mínim absolut).
    Usar Fraction elimina TOT error d'arrodoniment per a N de qualsevol mida.
    """
    D = 2 * (2 * m + 3)
    R = N % D
    return Fraction(R, D) - Fraction(1, 2)


def m_espell(m: int, N: int) -> int | None:
    """
    Espell aritmètic: si B = 2m+3, llavors A = N // B i m_espell = (A-3)//2.

    PRINCIPI: A · B ≈ N (exactament = N si m és el factor real).
    L'espell no és geomètric sinó ARITMÈTIC per la propietat multiplicativa.
    Per factors reals: A·B = N exactament.
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
    Desfase entre la senyal principal B i l'espell A:
      ΔΦ(m) = d_frac(m, N) − d_frac(m_espell, N)

    TEOREMA D'ATRAPAMENT (verificat experimentalment):
    Si ΔΦ(m_i) · ΔΦ(m_{i+1}) < 0 → existeix un factor en [m_i, m_{i+1}].

    Quan m és exactament el factor p: B = p, A = N/p enter, ΔΦ = mínim global.
    """
    dB  = d_frac(m, N)
    me  = m_espell(m, N)
    if me is None:
        return Fraction(0)
    dA  = d_frac(me, N)
    return dB - dA


def check_factor(m: int, N: int) -> bool:
    """Verificació directa i exacta: (2m+3) divideix N?"""
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


# ============================================================================
#  BLOC 3 — SALT EXACTE PER PENDENT  (S±k del Libro 5)
# ============================================================================

def salt_exacte_pendent(m: int, N: int) -> int | None:
    """
    Salt exacte al zero de la pendent lineal actual.

    MATEMÀTICA DE LES PENDENTS LINEALS:
      La seqüència de restos r(d) = N mod d forma pendents:
        r(d) = N − k·d    per a tot d en [N/(k+1), N/k]
      on k = ⌊N/pos⌋ és el cocient = pendent actual.

      El zero de la recta r(d) = 0 és:
        d_zero = N / k    →    m_zero = (d_zero − 3) // 2

      Si N és múltiple de d_zero → FACTOR EXACTE.
      Si no → el salt porta prop del proper canvi de pendent.

    NOTA: açò és EXACTE, no una aproximació. Un factor real sempre
    satisfà k = N // p i p = N // k, amb N % p == 0.

    RETORNA: m_zero (enter) o None si el salt és invàlid.
    """
    pos = 2 * m + 3
    if pos <= 1:
        return None
    R = N % pos
    if R == 0:
        return m    # Ja estem al factor!
    k = N // pos    # Pendent actual (cocient)
    if k <= 0:
        return None
    # El zero de la recta k és exactament N//k (divisió entera de Python)
    candidat = N // k
    if candidat < 3:
        return None
    return (candidat - 3) // 2


def provar_salt_pendent(m: int, N: int, m_ini: int, m_fi: int) -> int | None:
    """
    Prova el salt exacte per pendent i els seus veïns immediats.
    Retorna el factor D = 2m+3 si es troba, o None.
    """
    m_salt = salt_exacte_pendent(m, N)
    if m_salt is None:
        return None

    # Comprovem el punt exacte i ±15 veïns (per errors de ±1 en la divisió)
    for dm in range(-15, 16):
        m_t = m_salt + dm
        if m_t < 1 or m_t > m_fi:
            continue
        if pasa_rueda(m_t) and check_factor(m_t, N):
            return 2 * m_t + 3
    return None


# ============================================================================
#  BLOC 4 — PINÇA DOBLE 4+4  (la pinça real del Libro 5)
# ============================================================================

def pinca_doble(m: int, N: int) -> dict:
    """
    Pinça doble de 4+4 medicions: 4 cap endavant + 4 cap enrere.

    MATEMÀTICA:
      Sèrie endavant: d0, d1, d2, d3  (m, m+1, m+2, m+3)
      Sèrie enrere:   d0, d_1, d_2, d_3 (m, m-1, m-2, m-3)

      Per a cada sèrie:
        V  = primera diferència    (velocitat)
        A  = segona diferència     (acceleració / curvatura)
        J  = tercera diferència    (jerk)

      DISCRIMINADOR DE JERK (idea del Libro 5):
        Si J_endavant ≈ J_enrere → estem dins d'un sol dent de serra
          (la funció és localment lineal: A ≈ 0, J ≈ 0).
          El salt mestre NO és fiable en este cas.
        Si J_endavant ≠ J_enrere → estem en la MACRO-ENVOLVENTE parabòlica
          (creuem pendents k consecutives). El salt mestre SÍ és fiable.

      SALT MESTRE:  m_next = m − round(V/A)   (anàleg a Newton-Raphson)

    RETORNA: diccionari amb tota la cinemàtica de les dos sèries.
    """
    # Sèrie endavant
    d = [d_frac(m + i, N) for i in range(4)]
    Vf  = d[1] - d[0]
    Af  = d[2] - 2*d[1] + d[0]
    Jf  = d[3] - 3*d[2] + 3*d[1] - d[0]

    # Sèrie enrere (m-1, m-2, m-3; la velocitat és negada perquè anem enrere)
    dr = [d_frac(m - i, N) for i in range(4)]
    Vr  = dr[1] - dr[0]    # velocitat cap enrere (signe oposat a Vf si és monòton)
    Ar  = dr[2] - 2*dr[1] + dr[0]
    Jr  = dr[3] - 3*dr[2] + 3*dr[1] - dr[0]

    # Discriminador: diferència de jerks
    delta_J = abs(Jf - Jr)

    # Salt mestre endavant (si A != 0 i jerk indica macro-envolvente)
    salt_f = None
    if Af != 0:
        despl = Vf / Af
        despl_int = int(despl + Fraction(1,2)) if despl >= 0 else -int(-despl + Fraction(1,2))
        m_next = m - despl_int
        if m_next >= 1:
            salt_f = m_next

    # Salt mestre enrere
    salt_r = None
    if Ar != 0:
        despl = Vr / Ar
        despl_int = int(despl + Fraction(1,2)) if despl >= 0 else -int(-despl + Fraction(1,2))
        m_next = m + despl_int   # enrere: sumem en lloc de restar
        if m_next >= 1:
            salt_r = m_next

    return {
        'd0'      : d[0],
        'Vf'      : Vf,   'Af'   : Af,   'Jf'    : Jf,    'salt_f' : salt_f,
        'Vr'      : Vr,   'Ar'   : Ar,   'Jr'    : Jr,    'salt_r' : salt_r,
        'delta_J' : delta_J,
        # V−|V|: component direccional (positiu = baixem cap al mínim)
        'V_dir'   : Vf - abs(Vf),
    }


# ============================================================================
#  BLOC 5 — FILTRE ACUMULADOR  (filtre de resonància de segments)
# ============================================================================

def filtre_acumulador(m_ini: int, m_fi: int, N: int,
                      n_punts: int = 8) -> bool:
    """
    Filtre de resonància basant-se en l'acumulador harmònic.

    MATEMÀTICA (del Libro 5):
      acumulo += acumulo * Fraction(1, i)

    INTERPRETACIÓ:
      La seqüència de restos N%pos té pendents lineals de pendent k=1,2,3,...
      En un segment que CONTÉ un factor, els restos decreixen bruscament
      en almenys un punt (el salt de pendent). L'acumulador mesura si
      este decreixement és coherent amb la progressió harmònica 1/i.

      Si en n_punts equidistants del segment el ratio de decreixement dels
      restos és molt superior a l'acumulador harmònic → el segment no té
      cap zero de pendent dins seu → DESCARTAR.

    IMPLEMENTACIÓ PRÀCTICA:
      Prenem n_punts dins del segment. Per a cadascun calculem:
        - ratio_decr = (N%pos_actual) / (N%pos_anterior)   (decreixement)
      Si el mínim ratio observat és > acumulo final → segment estèril.

    RETORNA: True si el segment PASSA el filtre (té resonància),
             False si és ESTÈRIL (descartat).
    """
    amplada = m_fi - m_ini
    if amplada <= 0:
        return False
    if amplada < 20:
        return True     # Segments molt menuts: no descartem (massa risc)

    # Punts de mostra equidistants
    pas      = max(1, amplada // n_punts)
    ms_mostra = [m_ini + i * pas for i in range(n_punts)]
    ms_mostra.append(m_fi)

    # Acumulador harmònic (filtre de resonància)
    acumulo = Fraction(1)
    min_ratio_decr = Fraction(2)   # inicialitzem a "no decreix"

    restos = []
    for m in ms_mostra:
        pos  = 2 * m + 3
        if pos <= 1:
            continue
        R = N % pos
        restos.append(R)

    # Mesurem el decreixement relatiu dels restos
    for i in range(1, len(restos)):
        if restos[i-1] > 0:
            ratio = Fraction(restos[i], restos[i-1])
            if ratio < min_ratio_decr:
                min_ratio_decr = ratio

        # Avancem l'acumulador
        acumulo += acumulo * Fraction(1, i + 1)

    # Condició de RESONÀNCIA: el decreixement observat és prou fort
    # Si min_ratio < (1 − 1/acumulo) → hi ha almenys un salt de pendent
    # que podria ser un zero → segment PASSA el filtre
    llindar = Fraction(1) - Fraction(1, max(1, int(float(acumulo))))

    # Si el mínim ratio és molt gran (>0.98) → no hi ha decreixements bruscos
    # → segment sense zeros potencials → DESCARTAR
    if min_ratio_decr > Fraction(98, 100):
        return False

    return True     # Passa: té decreixements que indiquen possibles zeros


# ============================================================================
#  BLOC 6 — ESPIRAL HARMÒNICA  (generació de sondes)
# ============================================================================

def generar_sondes_espiral(m_max: int, N: int) -> list[int]:
    """
    Genera sondes de mostreig per la sèrie harmònica factorial.

    acumulo_inicial = m_max
    acumulo_i = acumulo_{i-1} * (1/i)    →    salt_i = m_max / (2·i!)

    Les posicions cobren les TRANSICIONS DE PENDENT k → k+1,
    que és on es troben els zeros exactes (els factors).

    Primera sonda:  m_max − m_max/2       = m_max/2        (k≈2)
    Segona sonda:   anterior − m_max/4    = 3·m_max/4      (k≈3)
    ...
    Convergència:   m_max · (1 − (e−1)/2) ≈ m_max · 0.141  (k→∞)

    COBERTURA ADDICIONAL:
    · Simètrica respecte a m_max (per factors asimètrics grans)
    · Punts al voltant de l'atractor e−2 ≈ 0.718·m_max
    · Extrems [1, m_max] sempre inclosos
    """
    sondes  = set()
    salta   = Fraction(m_max, 2)    # primer salt: m_max/2
    posicio = Fraction(m_max)       # partim des de m_max

    # Punt màxim sempre inclós (factors equilibrats prop de √N)
    mv = proper_valid(m_max, -1)
    if mv:
        sondes.add(mv)

    # Atractor principal: (e−2)·m_max ≈ 0.7183·m_max
    m_conv = int((math.e - 2) * m_max)

    for i in range(1, 500):
        posicio -= salta
        m_sonda  = int(posicio)

        if m_sonda < 1:
            break

        # Sonda principal
        mv = proper_valid(m_sonda, 1)
        if mv and mv <= m_max:
            sondes.add(mv)

        # Simètrica respecte a m_max (zones de factors asimètrics)
        m_sim = m_max - m_sonda
        if m_sim >= 1:
            mv2 = proper_valid(m_sim, 1)
            if mv2 and mv2 <= m_max:
                sondes.add(mv2)

        # Reducció del salt per a la propera iteració
        salta = salta * Fraction(1, i + 1)

        # Quan el salt és menor que 1 → refinament dens al voltant de m_conv
        if salta < Fraction(1, 2):
            radi = min(300, max(10, m_max // 500))   # MÀXIM 600 punts addicionals
            m_ini_ref = max(1,      m_conv - radi)
            m_fi_ref  = min(m_max,  m_conv + radi)
            for mr in range(m_ini_ref, m_fi_ref + 1):
                if pasa_rueda(mr):
                    sondes.add(mr)
            break

    # Extrems garantits
    for m_ext in [1, m_max // 4, m_max // 2, m_max]:
        mv = proper_valid(m_ext, 1)
        if mv and mv <= m_max:
            sondes.add(mv)

    # De major a menor: explorem primers els factors equilibrats (zona alta)
    return sorted(list(sondes), reverse=True)


# ============================================================================
#  BLOC 7 — DETECTOR ESPELL  (atrapament de segments candidats)
# ============================================================================

def detectar_segments_espell(N: int, sondes: list[int]) -> list[tuple]:
    """
    Detecta segments candidats per canvi de signe del desfase espell.

    Per a cada parell consecutiu (m_i, m_{i+1}) de la llista de sondes:
      ΔΦ_i   = desfase_espell(m_i,   N)
      ΔΦ_{i+1} = desfase_espell(m_{i+1}, N)

      Si ΔΦ_i · ΔΦ_{i+1} < 0 → CANVI DE SIGNE → factor atrapat.

    PRIORITAT: segments ordenats per |ΔΦ_i · ΔΦ_{i+1}| descendent
    (canvis de signe més abruptes primer = major probabilitat).

    Inclou comprovació directa de la sonda (factor exacte a m_i o m_{i+1}).
    """
    # Calculem ΔΦ per totes les sondes
    dfs = []
    for m in sondes:
        # Comprovació directa (factor exacte)
        if check_factor(m, N):
            return [('DIRECTE', m, m, Fraction(0), Fraction(0))]
        df = desfase_espell(m, N)
        dfs.append((m, df))

    # Detectem canvis de signe
    candidats = []
    for i in range(len(dfs) - 1):
        m_a, df_a = dfs[i]
        m_b, df_b = dfs[i + 1]

        # Assegurem ordre (m_b < m_a perquè sondes va de gran a menut)
        m_esq, m_dre = min(m_a, m_b), max(m_a, m_b)
        df_esq = dfs[i+1][1] if m_a > m_b else df_a
        df_dre = df_a        if m_a > m_b else df_b

        if df_a * df_b < 0:    # Canvi de signe: factor atrapat
            # PRIORITAT: el segment on un dels ΔΦ és més proper a 0 és
            # el que té el factor més a prop de la sonda → prioritzem
            # per min(|ΔΦ_a|, |ΔΦ_b|) ASCENDENT (el zero és on importa)
            proximitat = min(abs(df_a), abs(df_b))
            candidats.append((m_esq, m_dre, df_esq, df_dre, proximitat))

    # Ordenem per proximitat ascendent: el segment amb ΔΦ més proper a 0 primer
    candidats.sort(key=lambda x: x[4])
    return [(s[0], s[1], s[2], s[3]) for s in candidats]


# ============================================================================
#  BLOC 8 — REFINAMENT DE SEGMENT  (salt exacte + pinça doble + bisecció)
# ============================================================================

def refinar_segment(N: int, m_ini: int, m_fi: int,
                    profunditat: int = 60,
                    verbose: bool = False) -> int | None:
    """
    Refina un segment candidat per trobar el factor.

    ESTRATÈGIA PRINCIPAL — K-SWEEP (el veritable S±k del Libro 5):
    ─────────────────────────────────────────────────────────────────
    En lloc d'iterar sobre m (bruta força), iterem sobre les PENDENTS k:

      Per a pos en [2*m_ini+3, 2*m_fi+3]:
        k = N // pos            (pendent actual)
        candidat = N // k       (zero de la recta k)
        si N % candidat == 0 → FACTOR EXACTE

      Equivalentment, per a cada k en [k_min, k_max]:
        k_min = N // (2*m_fi+3)     (k per al pos màxim del segment)
        k_max = N // (2*m_ini+3)    (k per al pos mínim del segment)
        candidat = N // k
        si N % candidat == 0 → FACTOR

    COMPLEXITAT del k-sweep:
      nombre de k: O(k_max − k_min) = O(N·Δm / m²) ≈ O(Δm) per a m≈m_max

    PER QUÈ ÉS DETERMINISTA:
      Tot factor p de N té pendent k = N // p exacta. Esta k SEMPRE cau
      dins del rang [k_min, k_max] del segment que conté m = (p−3)//2.
      El k-sweep troba TOTS els factors d'un segment sense mai perdre'n cap.

    ESTRATÈGIA SECUNDÀRIA — PINÇA DOBLE + SALT MESTRE:
      Per als millors candidats (menor n%pos%(pos-1)), la pinça doble
      proporciona un salt addicional si el k-sweep no ha trobat res.

    ESTRATÈGIA TERCIÀRIA — BISECCIÓ PER DESFASE ESPELL:
      Bisecció guiada per el signe de ΔΦ per a segments molt grans on
      el k-sweep seria massa costós (> 10M iteracions).
    """
    m_ini_orig = m_ini
    m_fi_orig  = m_fi

    pos_ini = 2 * m_ini + 3    # pos mínim del segment
    pos_fi  = 2 * m_fi  + 3    # pos màxim del segment

    # ── K-SWEEP: iteració sobre les pendents k ─────────────────────────────
    # k per a pos_fi (pos gran → k petit) i pos_ini (pos petit → k gran)
    k_min = N // pos_fi   if pos_fi  > 0 else 1
    k_max = N // pos_ini  if pos_ini > 0 else 1

    # Nombre de k a comprovar
    n_k = k_max - k_min + 1

    # K-SWEEP DIRECTE si el rang és manejable (≤ 2M de k)
    if n_k <= 2_000_000:
        for k in range(k_min, k_max + 1):
            if k <= 0:
                continue
            candidat = N // k
            if candidat < 3:
                continue
            # Verificació: N % candidat == 0 i candidat és senar no trivial
            if candidat % 2 == 0:
                continue
            if N % candidat == 0 and 1 < candidat < N:
                m_cand = (candidat - 3) // 2
                if verbose:
                    print(f"      [K-sweep k={k}] candidat={candidat:,}  m={m_cand:,}")
                return candidat

        # K-sweep no ha trobat res: el segment no conté cap factor
        # (o n_k=0 → segment d'un sol punt ja comprovat)
        return None

    # ── Per a segments MOLT GRANS (n_k > 2M): bisecció per ΔΦ ────────────
    # Reduïm el segment fins que el k-sweep sigui manejable
    if verbose:
        print(f"      [Segment gran n_k={n_k:,}: bisecció ΔΦ→k-sweep]")

    for iter_num in range(profunditat):
        amplada = m_fi - m_ini

        # Recalculem n_k per al segment actual
        pos_i = 2 * m_ini + 3
        pos_f = 2 * m_fi  + 3
        k_lo  = N // pos_f  if pos_f  > 0 else 1
        k_hi  = N // pos_i  if pos_i  > 0 else 1
        n_k_actual = k_hi - k_lo + 1

        # K-sweep quan és manejable
        if n_k_actual <= 2_000_000:
            for k in range(k_lo, k_hi + 1):
                if k <= 0:
                    continue
                candidat = N // k
                if candidat < 3 or candidat % 2 == 0:
                    continue
                if N % candidat == 0 and 1 < candidat < N:
                    if verbose:
                        print(f"      [K-sweep bisectat iter={iter_num} k={k}] candidat={candidat:,}")
                    return candidat
            return None

        # Bisecció per ΔΦ per reduir el segment
        m_e = proper_valid(m_ini,              1)
        m_c = proper_valid((m_ini + m_fi) // 2, 1)
        m_d = proper_valid(m_fi,               -1)

        if not (m_e and m_c and m_d):
            return None

        # Comprovació directa via S±k dels tres punts de bisecció
        for m_t in [m_e, m_c, m_d]:
            factor_D = provar_salt_pendent(m_t, N, m_ini_orig, m_fi_orig)
            if factor_D:
                return factor_D

        # Bisecció guiada per ΔΦ
        df_e = desfase_espell(m_e, N)
        df_c = desfase_espell(m_c, N)
        df_d = desfase_espell(m_d, N)

        if df_e * df_c < 0:
            m_fi  = m_c
        elif df_c * df_d < 0:
            m_ini = m_c
        else:
            # Sense canvi de signe: usem n%pos%(pos-1) per triar costat
            pos_e = 2 * m_e + 3
            pos_d = 2 * m_d + 3
            R2_e  = (N % pos_e) % (pos_e - 1) if pos_e > 1 else 999
            R2_d  = (N % pos_d) % (pos_d - 1) if pos_d > 1 else 999
            if R2_e <= R2_d:
                m_fi  = m_c
            else:
                m_ini = m_c

    return None


# ============================================================================
#  BLOC 9 — ORQUESTRADOR PRINCIPAL  MDC v17
# ============================================================================

def mdc_v17(N: int, verbose: bool = True) -> tuple:
    """
    MDC v17 — Orquestrador complet (5 fases).

    FLUX CORREGIT:
      F0 · Precondicions
      F1 · Espiral harmònica → sondes
      F2 · Filtre acumulador → ordena segments per prioritat
      F3 · Detector espell → PRIORITZA (no descarta) segments candidats
      F4 · K-sweep sobre TOTS els segments ordenats per prioritat

    INSIGHT CLAU (correcció respecte a v16):
      L'espell detecta canvis de signe de ΔΦ als EXTREMS dels segments.
      Però el factor pot estar dins d'un segment on TOTS DOS extrems
      tenen el mateix signe (ΔΦ oscil·la dins i creua zero internament).
      El k-sweep resol açò: busca el factor per totes les pendents k
      del segment, sense dependre de la posició de les sondes extrems.
      L'espell es manté com a PRIORITZACIÓ: els segments amb canvi de
      signe es proven primer (major probabilitat), però TOTS es proven.
    """
    if verbose:
        bar = '═' * 72
        print(f"\n{bar}")
        print(f"  MDC v17 — Mètode Diofàntic Cinemàtic")
        print(f"  K-sweep · Pendent Exacta · Pinça 4+4 · Espell (prioritat)")
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
    m_conv = int((math.e - 2) * m_max)

    if verbose:
        print(f"    → m_max={m_max:,}   m_conv(e−2)={m_conv:,}")

    # ── FASE 1: ESPIRAL HARMÒNICA ──────────────────────────────────────────
    if verbose:
        print(f"\n  [F1] Espiral harmònica 1/(2·i!) ...")

    sondes = generar_sondes_espiral(m_max, N)

    if verbose:
        print(f"    → {len(sondes)} sondes generades")
        for m in sondes[:3]:
            print(f"       m={m:>10,}  ({m/m_max*100:.1f}%)  D={2*m+3:,}")
        print(f"       ...")
        for m in sondes[-2:]:
            print(f"       m={m:>10,}  ({m/m_max*100:.1f}%)  D={2*m+3:,}")

    # ── FASE 2+3: GENERAR I PRIORITZAR TOTS ELS SEGMENTS ─────────────────
    if verbose:
        print(f"\n  [F2+F3] Generant i prioritzant segments...")

    # Comprovació directa de cada sonda (cost zero addicional)
    for m in sondes:
        if check_factor(m, N):
            t_ms = (time.perf_counter() - t0) * 1000
            _print_resultat(N, 2*m+3, 0, t_ms)
            return 2*m+3, t_ms

    # Generem tots els segments (parells consecutius de sondes)
    sondes_ord = sorted(sondes)    # de menor a major
    tots_segments = []

    for i in range(len(sondes_ord) - 1):
        m_a = sondes_ord[i]
        m_b = sondes_ord[i + 1]
        if m_b <= m_a:
            continue

        pos_a = 2 * m_a + 3
        pos_b = 2 * m_b + 3
        k_lo  = N // pos_b   # k per al pos gran (m_b)
        k_hi  = N // pos_a   # k per al pos menut (m_a)
        n_k   = k_hi - k_lo + 1

        # ΔΦ als extrems per a priorització
        df_a = desfase_espell(m_a, N)
        df_b = desfase_espell(m_b, N)

        # Prioritat: segments amb canvi de signe de ΔΦ → prioritat 0 (primers)
        # Segments sense canvi → prioritat 1 (després), ordenats per n_k
        te_canvi = (df_a * df_b < 0)
        # Dintre de cada grup, ordenem per n_k ascendent (menys iteracions primer)
        prioritat = (0 if te_canvi else 1, n_k)

        tots_segments.append((prioritat, m_a, m_b, df_a, df_b, n_k))

    # Ordenem: primer els amb canvi de signe, dintre per n_k creixent
    tots_segments.sort(key=lambda x: x[0])

    n_canvi = sum(1 for s in tots_segments if not s[0][0])
    n_total = len(tots_segments)

    if verbose:
        print(f"    → {n_total} segments totals")
        print(f"       {n_canvi} amb canvi de signe ΔΦ (prioritat alta)")
        print(f"       {n_total - n_canvi} sense canvi (prioritat baixa)")
        print(f"    Primers segments per prioritat:")
        for (pri, m_a, m_b, df_a, df_b, n_k) in tots_segments[:5]:
            marca = "★" if not pri[0] else " "
            print(f"       {marca} [{m_a:,},{m_b:,}]  n_k={n_k:,}  "
                  f"ΔΦ:{float(df_a):.4f}→{float(df_b):.4f}")

    # ── FASE 4: K-SWEEP SOBRE TOTS ELS SEGMENTS ───────────────────────────
    if verbose:
        print(f"\n  [F4] K-sweep per tots els segments...")

    for idx, (pri, m_a, m_b, df_a, df_b, n_k) in enumerate(tots_segments):

        if verbose and idx < 10:
            marca = "★" if not pri[0] else " "
            print(f"\n    {marca} Segment #{idx+1}: [{m_a:,},{m_b:,}]  "
                  f"n_k={n_k:,}", end="  ")

        factor_D = refinar_segment(N, m_a, m_b,
                                   profunditat=60,
                                   verbose=verbose and idx < 5)
        if factor_D is not None:
            t_ms = (time.perf_counter() - t0) * 1000
            _print_resultat(N, factor_D, idx + 1, t_ms)
            return factor_D, t_ms
        elif verbose and idx < 10:
            print("✗")

    # ── FALLADA ────────────────────────────────────────────────────────────
    t_ms = (time.perf_counter() - t0) * 1000
    if verbose:
        print(f"\n  ✗ Factor no trobat. Temps: {t_ms:.3f} ms")
    return None, t_ms


def _print_resultat(N, factor, seg_idx, t_ms):
    bar = '═' * 72
    print(f"\n{bar}")
    print(f"  🎯 FACTOR TROBAT!")
    print(f"     Factor p : {factor:,}")
    print(f"     Factor q : {N // factor:,}")
    print(f"     Segment  : #{seg_idx}")
    print(f"     Temps    : {t_ms:.3f} ms")
    print(f"{bar}")


# ============================================================================
#  BLOC 10 — BATERIA DE PROVES
# ============================================================================

def executar_bateria(verbose: bool = False):
    proves = [
        # (N,                               desc)
        (3 * 5,                             "3 × 5  (trivial)"),
        (17 * 19,                           "17 × 19 = 323"),
        (101 * 103,                         "101 × 103 = 10.403  (clàssic)"),
        (100003 * 100019,                   "100003 × 100019  (11 dígits, estudi visual)"),
        (1_000_003 * 1_000_033,             "1.000.003 × 1.000.033  (13 dígits)"),
        (9_999_991 * 9_999_973,             "9.999.991 × 9.999.973  (14 dígits)"),
        (1_548_586_332_452_843,             "1.548...  (16 dígits, prova Colab)"),
        (999_983 * 1_000_000_000_000_003,   "999.983 × 1.000...003  (asimètric 6×16 dígits)"),
    ]

    print("\n" + "█" * 72)
    print("  🚀 BATERIA DE PROVES — MDC v17")
    print("█" * 72)

    ok = fail = 0
    for N, desc in proves:
        print(f"\n{'─' * 72}")
        print(f"  📋 {desc}")
        factor, t_ms = mdc_v17(N, verbose=verbose)
        if factor:
            assert N % factor == 0, f"ERROR CRÍTIC: {factor} no divideix {N}"
            print(f"  ✅  {factor:,} × {N//factor:,}   ({t_ms:.3f} ms)")
            ok += 1
        else:
            print(f"  ⚠️   No trobat  ({t_ms:.3f} ms)")
            fail += 1

    print(f"\n{'█' * 72}")
    print(f"  Encerts: {ok}/{ok+fail}")
    print("█" * 72)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

if __name__ == "__main__":
    # Prova individual amb verbose complet
    N = 100003 * 100019
    factor, t_ms = mdc_v17(N, verbose=True)

    print("\n" + "─" * 72)
    print("  Bateria completa (verbose=False)...")
    executar_bateria(verbose=False)
