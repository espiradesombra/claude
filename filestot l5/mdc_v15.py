# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC v15 — MÈTODE DIOFÀNTIC CINEMÀTIC
  Versió Didàctica i d'Aritmètica Exacta
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)
  VERSIÓ: 15

  NOVETATS RESPECTE v14:
  ─────────────────────────────────────────────────────────────────────────
  1. ARITMÈTICA EXACTA (Fraction):
     Tots els càlculs crítics (funció d'ona, filtre de mínims, offset
     vectorial) usen `fractions.Fraction` en lloc de `float`. Açò elimina
     completament els errors d'arrodoniment que apareixien amb N de molts
     dígits.

  2. SEGMENTACIÓ ARMÒNICA AMB COMPARACIÓ DE VÀLIDS:
     El rang de cerca es divideix en segments usant la progressió harmònica
         acumulo(i) = acumulo(i-1) + acumulo(i-1) * Fraction(1, i)
     Per a cadascun d'estos segments es fa una COMPARACIÓ explícita:
     ¿hi ha valors m que passen la roda primorial? (= "i vàlids a la
     seqüència"). Si la resposta és NO → el segment s'omet completament,
     estalviant feina inútil.

  3. CODI DIDÀCTIC:
     Cada funció explica la matemàtica que implementa. L'objectiu és que
     el codi siga llegible com un article tècnic.

  PRINCIPIS MATEMÀTICS CENTRALS:
  ─────────────────────────────────────────────────────────────────────────
  · Tot divisor senar de N pren la forma D = 2k+1. Paramètritzem com
    D(m) = 2m+3 (forma que evita trivials) → factors no trivials ↔ m ≥ 1.

  · La FUNCIÓ D'ONA mesura la distància al "zero de fase":
      d(m) = N mod D(m) / D(m) − 1/2    (Fraction exacta)
    Quan d(m) → 0: N mod D(m) ≈ D(m)/2, indici fort que D(m) | N.
    Quan N % (2m+3) == 0 exactament: hem trobat el factor.

  · La CINEMÀTICA (V, A) permet saltar al vèrtex de la paràbola local:
      m_next = m − V/A    (Salt Mestre, exacte amb Fraction)
    en lloc d'iterar m a m, que seria força bruta.

  · L'ACUMULADOR HARMÒNIC 1/i defineix les fronteres de segment i actua
    com a filtre macroscòpic: segments sense candidats vàlids s'ometen.
==============================================================================
"""

from fractions import Fraction
import math
import time


# ============================================================================
#  BLOC 0 — CONSTANTS GLOBALS
# ============================================================================
# Roda primorial fins a p=23. El producte (primorial) és 223.092.870.
# La roda garanteix que D(m)=2m+3 no és múltiple de cap d'estos primers,
# descartant ~92 % del espai de cerca sense perdre cap factor (suposant
# que N ja ha sigut dividit pels primers petits en la fase de precondicions).

PRIMOS   = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)   # 223092870


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL  (Filtre topològic)
# ============================================================================

def pasa_rueda(m: int) -> bool:
    """
    Comprova si D(m) = 2m+3 és COPRIMER amb el primorial.

    MATEMÀTICA:
      Si mcd(2m+3, PRIMORIAL) > 1, aleshores 2m+3 és divisible per algun
      primer petit. Per tant NO pot ser un factor no trivial de N (assumint
      que la fase de precondicions ja ha eliminat estos primers de N).

    RETORNA: True  → m és un candidat vàlid (val la pena analitzar-lo)
             False → m és un candidat descartable
    """
    if m < 1:
        return False
    return math.gcd(2 * m + 3, PRIMORIAL) == 1


def alinear_rueda(m: int) -> int | None:
    """
    Ajusta m cap avall fins al primer valor vàlid de la roda.

    Ús típic: quan generem un pilar en una posició harmònica exacta,
    és possible que eixa posició no passe la roda. Baixem d'un en un
    fins a trobar el candidat vàlid més proper.

    RETORNA: m vàlid ≥ 1, o None si no n'hi ha cap en [1, m_original].
    """
    while m >= 1 and not pasa_rueda(m):
        m -= 1
    return m if m >= 1 else None


# ============================================================================
#  BLOC 2 — FUNCIÓ D'ONA  (El cor del MDC)
# ============================================================================

def d_frac(m: int, N: int) -> Fraction:
    """
    Funció d'ona del MDC: distància signada al zero de fase.

    MATEMÀTICA:
      D  = 2*(2m+3)         denominador de la forma simètrica
      R  = N mod D          reste exacte (operació entera pura)
      d  = Fraction(R, D) − Fraction(1, 2)

    La tria D = 2*(2m+3) (no D = 2m+3 directament) centra la "fase" a 0.5:
      · Quan R = 0   → d = −1/2  (N és múltiple de D, però D és parell)
      · Quan R = D/2 → d = 0     → PUNT D'INTERÈS MÀXIM
      · Quan R = D   → d = +1/2

    La verificació final de factor usa sempre N % (2m+3) == 0, que és
    l'aritmètica directa. La funció d'ona ens dona la direcció i curvatura
    per navegar fins a eixe punt sense iterar tots els m.

    AVANTATGE de Fraction: amb float, per a N de 50+ dígits apareixien
    errors d'arrodoniment que confondrien la cinemàtica. Fraction és exacta.
    """
    D = 2 * (2 * m + 3)   # Denominador simètric (sempre parell)
    R = N % D              # Reste exacte (operació entera, sense error)
    return Fraction(R, D) - Fraction(1, 2)


def check_factor(m: int, N: int) -> bool:
    """
    Verificació directa i EXACTA de si (2m+3) divideix N.

    Esta és la condició definitiva de factorització. La cinemàtica ens
    guia fins ací, però la paraula final l'té sempre esta divisió entera.
    """
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


# ============================================================================
#  BLOC 3 — CINEMÀTICA  (Velocitat, Acceleració, Salt Mestre)
# ============================================================================

def mesurar_cinematica(m: int, N: int) -> tuple[Fraction, Fraction, Fraction]:
    """
    Mesura velocitat V i acceleració A de la funció d'ona en el punt m.

    MATEMÀTICA (Diferències Finites Exactes):
      d0 = d_frac(m)       posició actual
      d1 = d_frac(m+1)     posició al pas m+1
      d2 = d_frac(m+2)     posició al pas m+2

      V  = d1 − d0                  1a diferència  (velocitat)
      A  = d2 − 2·d1 + d0           2a diferència  (acceleració / curvatura)

    Amb tres punts suposem localment una paràbola:
      d(m+k) ≈ d0 + V·k + (A/2)·k²

    El vèrtex d'esta paràbola (on la funció és mínima en valor absolut) és:
      k* = −V/A   →   m* = m + k* = m − V/A

    RETORNA: (d0, V, A) com a Fraction exactes (sense arrodoniment).
    """
    d0 = d_frac(m,     N)
    d1 = d_frac(m + 1, N)
    d2 = d_frac(m + 2, N)

    V = d1 - d0                             # Velocitat (Fraction)
    A = d2 - Fraction(2) * d1 + d0         # Acceleració (Fraction)

    return d0, V, A


def salt_mestre(m: int, V: Fraction, A: Fraction) -> int | None:
    """
    Salt Mestre: projecta m cap al vèrtex de la paràbola local.

    MATEMÀTICA:
      Si A ≠ 0: m_next = m − round(V/A)
      Si A = 0: la funció és lineal en este interval, no hi ha vèrtex →
                retornem None (no hi ha salt útil).

    L'arrodoniment round() es fa amb Fraction per evitar biaixos de float.
    Resultat: el m enter més proper al vèrtex analític.

    RETORNA: m_next (enter ≥ 1) o None si A == 0.
    """
    if A == Fraction(0):
        return None   # Curvatura zero: no hi ha vèrtex, no hi ha salt

    # Fracció exacta del desplaçament fins al vèrtex
    despl = V / A

    # Arrodoniment exacte amb Fraction: equivalent a round() però sense float
    if despl >= 0:
        despl_int = int(despl + Fraction(1, 2))
    else:
        despl_int = -int(-despl + Fraction(1, 2))

    m_next = m - despl_int
    return m_next if m_next >= 1 else None


# ============================================================================
#  BLOC 4 — FILTRE DE MÍNIMS  (Profunditat del Reste)
# ============================================================================

def profunditat_minimo(m: int, N: int) -> Fraction:
    """
    Calcula la "profunditat" del mínim de reste per al candidat m.

    MATEMÀTICA (del Libro 5 / MDC v4.2):
      pos   = 2m+3           el candidat a divisor
      R     = N mod pos      reste de la divisió
      R2    = R mod (pos-1)  "mínim de la recta" (doble aritmètica modular)
      prof  = Fraction(R2, pos)   distància normalitzada al zero

    INTERPRETACIÓ:
      Quan pos divideix N → R = 0 → R2 = 0 → prof = 0 (mínim absolut).
      Com més prop de 0 estigui prof, més "interessant" és este candidat.
      Ordenar per prof ascendent equival a ordenar per "probabilitat de
      ser un factor" descendent.

    AVANTATGE de Fraction vs float:
      En v14, `profundidad = minimo_recta / pos` usava divisió entera de
      Python (que retorna float). Per a pos de 30+ dígits, la mantissa de
      64 bits de float no tenia prou precisió. Fraction és exacta.
    """
    pos = 2 * m + 3
    if pos <= 1:
        return Fraction(1)          # Cas degeneratiu: màxima "distància"
    R  = N % pos
    R2 = R % (pos - 1)
    return Fraction(R2, pos)        # Fraction exacta


def offset_vectorial(m: int, N: int, variacio: Fraction) -> int | None:
    """
    Calcula l'offset vectorial per al "Disparo Vectorial" (v14 → v15).

    MATEMÀTICA (versió Fraction exacta de v14):
      pos       = 2m+3
      offset    = (pos − N + variacio) / 2    (Fraction exacta)
      m_impacte = m − round(offset)

    En v14 esta operació usava `float(...)`, que introduïa error per a
    N grans. Ací usem Fraction per tindre el resultat exacte.

    RETORNA: m_impacte (enter ≥ 1) o None si el resultat és invàlid.
    """
    pos      = 2 * m + 3
    # Operació completament en Fraction: sense conversió a float
    offset_f = (Fraction(pos) - Fraction(N) + variacio) / Fraction(2)

    # Arrodoniment exacte
    if offset_f >= 0:
        offset_i = int(offset_f + Fraction(1, 2))
    else:
        offset_i = -int(-offset_f + Fraction(1, 2))

    m_imp = m - offset_i
    return m_imp if m_imp >= 1 else None


# ============================================================================
#  BLOC 5 — SEGMENTACIÓ ARMÒNICA  (La novetat central de v15)
# ============================================================================

def generar_segments_armonics(N: int) -> list[tuple]:
    """
    Divideix el rang de cerca [1, m_max] en segments usant la progressió
    harmònica de l'acumulador del Libro 5:

        acumulo(2)  = m_max
        acumulo(i)  = acumulo(i-1) + acumulo(i-1) * Fraction(1, i)
                    = acumulo(i-1) * Fraction(i+1, i)

    Desenrotllant:
        acumulo(i) = m_max * Fraction(3,2) * Fraction(4,3) * ... * Fraction(i+1,i)
                   = m_max * Fraction(i+1, 2)       [producte telescòpic]

    Les FRONTERES de segment s'obtenen invertint: on acumulo(i) / m_max = k/2
    per a k enter. En pràctica, generem les fronteres harmòniques simples:
        frontera(i) = m_max // i     (harmònica: m_max/1, m_max/2, m_max/3, ...)

    El segment S_i és l'interval:
        S_i = [m_max // (i+1) , m_max // i]

    PROPIETAT CLAU: Els segments amb i petit (S_1, S_2...) cobrixen la part
    alta del rang [√N/2, √N], on es troben els factors "equilibrats" (els dos
    factors de N similars en magnitud). Els segments amb i gran cobrixen la
    part baixa, on es troben factors molt asimètrics.

    RETORNA: llista de tuples (m_inici, m_fi, fraccio_rang, i_segment)
             ordenada per prioritat descendent (S_1 primer).
    """
    m_max = (math.isqrt(N) - 3) // 2
    if m_max < 1:
        return []

    segments = []
    i = 1
    while True:
        m_fi    = m_max // i          # Frontera superior (inclosa)
        m_inici = m_max // (i + 1)    # Frontera inferior (inclosa)

        if m_fi < 1:
            break
        if m_inici < 1:
            m_inici = 1               # Truncar al límit inferior del rang

        if m_fi > m_inici:            # Segment no buit en termes d'enters
            # Fracció exacta (Fraction) de la longitud d'este segment
            fraccio = Fraction(m_fi - m_inici, m_max)
            segments.append((m_inici, m_fi, fraccio, i))

        i += 1
        if i > 200 or (m_max // i) < 1:
            break                     # Segments massa petits: no aporten res

    return segments


def comparar_valids_en_segment(m_inici: int, m_fi: int,
                               n_mostra: int = 40) -> tuple[int, list[int]]:
    """
    COMPARACIÓ EXPLÍCITA: ¿hi ha valors m vàlids (de la roda) en [m_inici, m_fi]?

    Esta és la novetat central de v15. En v14, la xarxa intersticial
    generava punts sense comprovar si el segment tenia candidats vàlids.
    Ací fem la comprovació EXPLÍCITA per a cada segment:

      1. Prenem n_mostra punts equidistants dins de [m_inici, m_fi].
      2. Comprovem quants passen pasa_rueda().
      3. Si cap passa → el segment és "buit" → s'omet.
      4. Si algun passa → el segment és "vàlid" → es processa.

    PER QUÈ UNA MOSTRA EN LLOC DE TOTS ELS PUNTS?
      La roda primorial (p23) té cicle de longitud PRIMORIAL ≈ 2.2·10⁸.
      Dins d'un cicle, exactament φ(PRIMORIAL) valors passen la roda.
      La densitat és ~8.07 %, és a dir 1 de cada 12.4 enters aprox. Per
      tant, en un segment de longitud ≥ 50 quasi sempre hi ha candidats.
      En segments molt estrets (longitud < 13), pot no haver-n'hi cap.

    RETORNA: (n_valids_en_mostra, llista_de_m_valids)
    """
    amplada = m_fi - m_inici
    if amplada <= 0:
        return 0, []

    # Pas entre punts de mostra (com a mínim 1)
    pas = max(1, amplada // n_mostra)

    valids = []
    m = m_inici
    while m <= m_fi:
        if pasa_rueda(m):
            valids.append(m)
        m += pas

    # Assegurar que el límit superior sempre es comprova
    if m_fi not in valids and pasa_rueda(m_fi):
        valids.append(m_fi)

    return len(valids), valids


def filtrar_segments_amb_valids(segments: list[tuple],
                                n_mostra: int = 40) -> tuple[list, list]:
    """
    Aplica la COMPARACIÓ DE VÀLIDS a tots els segments i els classifica.

    PER A CADA SEGMENT:
      → Crida comparar_valids_en_segment()
      → Si té ≥ 1 candidat vàlid: va a 'segments_valids'
      → Si té 0 candidats vàlids: va a 'segments_omesos'

    RETORNA:
      segments_valids : llista de (m_inici, m_fi, fraccio, i, llista_m_valids)
      segments_omesos : llista de (m_inici, m_fi, fraccio, i)
    """
    segments_valids = []
    segments_omesos = []

    for (m_inici, m_fi, fraccio, i) in segments:
        n_v, llista_m = comparar_valids_en_segment(m_inici, m_fi, n_mostra)

        if n_v > 0:
            segments_valids.append((m_inici, m_fi, fraccio, i, llista_m))
        else:
            # Segment buit: s'omet completament (estalvi computacional)
            segments_omesos.append((m_inici, m_fi, fraccio, i))

    return segments_valids, segments_omesos


# ============================================================================
#  BLOC 6 — ESCANEIG PER SEGMENT  (Cinemàtica + Disparo Dual)
# ============================================================================

def escanejar_segment(N: int, m_inici: int, m_fi: int,
                      llista_m_valids: list[int],
                      verbose: bool = True) -> int | None:
    """
    Escandeix un segment buscant un factor de N, usant la cinemàtica MDC.

    ESTRATÈGIA PER A CADA PUNT m DE LA LLISTA DE VÀLIDS:

      ETAPA 1 — Verificació directa:
        Si N % (2m+3) == 0 → FACTOR TROBAT. Fi immediat.

      ETAPA 2 — Mesura cinemàtica (Fraction exacta):
        Calcular (d0, V, A) per al punt m.
        Si d0 és "quasi zero" (|d0| < umbral) → candidat fort.

      ETAPA 3 — Salt Mestre:
        Si A ≠ 0: calcular m_next = m − round(V/A).
        Si m_next és dins del segment → explorar [m_next−2, m_next+2].

      ETAPA 4 — Disparo Vectorial (Fraction):
        Calcular l'offset vectorial (versió exacta de v14) i verificar.

    RETORNA: el factor trobat (enter > 1) o None.
    """
    # Umbral de "quasi zero" per identificar candidats forts cinemàticament
    # Usem Fraction(1, 10) = 0.1 com a llindar (ajustable)
    UMBRAL_FASE = Fraction(1, 10)

    for m in llista_m_valids:

        # ── ETAPA 1: Verificació directa ──────────────────────────────────
        if check_factor(m, N):
            if verbose:
                print(f"      ✓ [Directe] Factor a m={m:,}  →  D={2*m+3:,}")
            return 2 * m + 3

        # ── ETAPA 2: Mesura cinemàtica (Fraction exacta) ──────────────────
        d0, V, A = mesurar_cinematica(m, N)

        # Candidat fort si la fase és molt propera a zero
        if abs(d0) < UMBRAL_FASE:
            # Comprovem també m+1 i m+2 (els tres punts usats en la mesura)
            for dm in range(3):
                m_t = m + dm
                if pasa_rueda(m_t) and check_factor(m_t, N):
                    if verbose:
                        print(f"      ✓ [Fase~0] Factor a m={m_t:,}  →  D={2*m_t+3:,}")
                    return 2 * m_t + 3

        # ── ETAPA 3: Salt Mestre ──────────────────────────────────────────
        m_proj = salt_mestre(m, V, A)

        if m_proj is not None and m_inici <= m_proj <= m_fi:
            # Exploració local al voltant de la projecció (±2 passos)
            for dm in range(-2, 3):
                m_t = m_proj + dm
                if m_t >= 1 and pasa_rueda(m_t):
                    if check_factor(m_t, N):
                        if verbose:
                            print(f"      ✓ [Salt Mestre] Factor a m={m_t:,}"
                                  f"  (projecció des de m={m:,})")
                        return 2 * m_t + 3

        # ── ETAPA 4: Disparo Vectorial (Fraction exacta) ──────────────────
        # Calculem la variació com Fraction per a l'offset vectorial
        pos      = 2 * m + 3
        R        = N % pos
        variacio = Fraction(R, pos - 1) if pos > 1 else Fraction(0)

        m_imp = offset_vectorial(m, N, variacio)
        if m_imp is not None and pasa_rueda(m_imp):
            if check_factor(m_imp, N):
                if verbose:
                    print(f"      ✓ [Vectorial] Factor a m={m_imp:,}"
                          f"  (offset des de m={m:,})")
                return 2 * m_imp + 3

    return None   # Cap factor trobat en este segment amb esta mostra


# ============================================================================
#  BLOC 7 — FILTRE DE ZONES CALENTES  (Priorització per Profunditat)
# ============================================================================

def prioritzar_candidats(N: int, llista_m: list[int],
                         top_k: int = 50) -> list[tuple]:
    """
    Ordena els candidats de la llista per "profunditat de mínim" (Fraction).

    MATEMÀTICA:
      profunditat(m) = Fraction( (N % (2m+3)) % (2m+2), 2m+3 )

    Candidats amb profunditat menor → restes més properes a zero → primers
    de la llista d'exploració.

    NOTA v15: En v14 esta operació usava divisió float. Ací tot és Fraction.

    RETORNA: llista de (profunditat_Fraction, m, variacio_Fraction) ordenada
             per profunditat ascendent, limitada a top_k elements.
    """
    zonas = []
    for m in llista_m:
        pos = 2 * m + 3
        if pos <= 1:
            continue
        R        = N % pos
        R2       = R % (pos - 1)
        prof     = Fraction(R2, pos)        # Profunditat exacta (Fraction)
        variacio = Fraction(R, pos - 1)     # Variació exacta (Fraction)
        zonas.append((prof, m, variacio))

    zonas.sort(key=lambda x: x[0])         # Menor profunditat = més prioritari
    return zonas[:top_k]


# ============================================================================
#  BLOC 8 — ORQUESTRADOR PRINCIPAL  MDC v15
# ============================================================================

def mdc_v15(N: int, n_mostra: int = 40, verbose: bool = True) -> tuple:
    """
    MDC v15 — Orquestrador principal.

    FLUX COMPLET:
    ─────────────────────────────────────────────────────────────────────────
    FASE 0: Precondicions
      · Factors trivials (primers petits de la roda)
      · Quadrat perfecte

    FASE 1: Segmentació Armònica (Fraction exacta)
      · Genera les fronteres harmòniques m_max//i
      · Cada segment representa un interval [m_max//(i+1), m_max//i]

    FASE 2: Comparació de Vàlids (la novetat v15)
      · Per a cada segment: ¿hi ha m que passen la roda primorial?
      · Segments buits → OMESOS (estalvi computacional)
      · Segments vàlids → PROCESSATS

    FASE 3: Priorització per Profunditat (Fraction)
      · Dins de cada segment vàlid, ordena els candidats per profunditat
        de mínim (criteri cinemàtic del MDC)

    FASE 4: Escaneig Cinemàtic (Fraction + Salt Mestre + Disparo Vectorial)
      · Per cada candidat prioritzat: Directe → Fase~0 → Salt Mestre → Vec.

    ARGUMENTS:
      N        : enter a factoritzar (ha de ser compost, senar, > 8)
      n_mostra : nombre de punts de mostra per segment (per a la comparació)
      verbose  : si True, imprimeix tot el procés pas a pas

    RETORNA: (factor, temps_ms) on factor és un divisor no trivial de N,
             o (None, temps_ms) si no s'ha trobat.
    """

    if verbose:
        ample = 72
        print(f"\n{'═'*ample}")
        print(f"  MDC v15 — Mètode Diofàntic Cinemàtic (Aritmètica Exacta)")
        print(f"  N = {N:,}")
        print(f"  ({len(str(N))} dígits)")
        print(f"{'═'*ample}")

    t0 = time.perf_counter()

    # ────────────────────────────────────────────────────────────────────────
    # FASE 0: PRECONDICIONS
    # ────────────────────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [FASE 0] Precondicions...")

    # 0a. Factors trivials: primers de la roda
    for p in PRIMOS:
        if N % p == 0 and p < N:
            t_ms = (time.perf_counter() - t0) * 1000
            if verbose:
                print(f"    → Factor trivial: {p}  (N = {p} × {N//p})")
            return p, t_ms

    # 0b. Quadrat perfecte: N = r² → factor = r
    r = math.isqrt(N)
    if r * r == N:
        t_ms = (time.perf_counter() - t0) * 1000
        if verbose:
            print(f"    → N és quadrat perfecte: N = {r}²")
        return r, t_ms

    m_max = (math.isqrt(N) - 3) // 2
    if verbose:
        print(f"    → Rang MDC: m ∈ [1, {m_max:,}]")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 1: SEGMENTACIÓ ARMÒNICA
    # ────────────────────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [FASE 1] Segmentació armònica (acumulo += acumulo·1/i)...")

    segments = generar_segments_armonics(N)

    if verbose:
        print(f"    → {len(segments)} segments generats")
        # Mostrar els primers 6 segments per il·lustrar l'estructura
        for (mi, mf, frac, i) in segments[:6]:
            pct = float(frac) * 100
            print(f"       S_{i:>2}: [{mi:>12,} , {mf:>12,}]  "
                  f"({pct:.2f}% del rang)  longitud={mf-mi:,}")
        if len(segments) > 6:
            print(f"       ... i {len(segments)-6} segments més de menor longitud")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 2: COMPARACIÓ — ¿HI HA i VÀLIDS EN CADA SEGMENT?
    # ────────────────────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [FASE 2] Comparació de candidats vàlids (roda p23)...")
        print(f"    (mostra de {n_mostra} punts per segment)")

    segs_valids, segs_omesos = filtrar_segments_amb_valids(segments, n_mostra)

    if verbose:
        print(f"    → Segments AMB candidats vàlids : {len(segs_valids)}")
        print(f"    → Segments OMESOS (sense vàlids): {len(segs_omesos)}")
        if segs_omesos:
            ids_omesos = [s[3] for s in segs_omesos]
            print(f"       IDs omesos: {ids_omesos}")

    # ────────────────────────────────────────────────────────────────────────
    # FASES 3+4: PRIORITZACIÓ I ESCANEIG PER SEGMENT
    # ────────────────────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [FASES 3+4] Priorització i escaneig cinemàtic per segment...")

    for (m_inici, m_fi, fraccio, i, llista_m_valids) in segs_valids:

        if verbose:
            print(f"\n  ── Segment S_{i}: [{m_inici:,}, {m_fi:,}] "
                  f"({len(llista_m_valids)} candidats a la mostra) ──")

        # FASE 3: Prioritzar per profunditat de mínim (Fraction)
        candidats_prioritzats = prioritzar_candidats(N, llista_m_valids, top_k=50)

        if verbose and candidats_prioritzats:
            m_top = candidats_prioritzats[0][1]
            p_top = float(candidats_prioritzats[0][0])
            print(f"    → Millor candidat: m={m_top:,}  profunditat={p_top:.6f}")

        # FASE 4: Escaneig cinemàtic sobre els candidats prioritzats
        llista_m_ordenats = [c[1] for c in candidats_prioritzats]
        factor = escanejar_segment(N, m_inici, m_fi, llista_m_ordenats, verbose)

        if factor is not None:
            t_ms = (time.perf_counter() - t0) * 1000
            if verbose:
                print(f"\n{'═'*72}")
                print(f"  🎯 FACTOR TROBAT!")
                print(f"     Factor p : {factor:,}")
                print(f"     Factor q : {N // factor:,}")
                print(f"     Segment  : S_{i}")
                print(f"     Temps    : {t_ms:.3f} ms")
                print(f"{'═'*72}")
            return factor, t_ms

    # ────────────────────────────────────────────────────────────────────────
    # FALLADA: cap segment ha lliurat el factor
    # ────────────────────────────────────────────────────────────────────────
    t_ms = (time.perf_counter() - t0) * 1000
    if verbose:
        print(f"\n  ✗ Factor no trobat amb la mida de mostra actual (n_mostra={n_mostra}).")
        print(f"    Prova a augmentar n_mostra o usar la bateria de verificació.")
        print(f"    Temps total: {t_ms:.3f} ms")
    return None, t_ms


# ============================================================================
#  BLOC 9 — BATERIA DE PROVES
# ============================================================================

def executar_bateria(verbose_global: bool = True):
    """
    Bateria de proves per verificar la correcció i robustesa del MDC v15.

    S'inclouen casos de dificultat creixent:
      · Semiprimers trivials i petits → verifiquen la correcció bàsica
      · Semiprimers mitjans           → comproven la cinemàtica
      · Semiprimers grans             → estic el rendiment

    Per a cada cas, SEMPRE es verifica que factor × (N/factor) == N.
    """

    proves = [
        # (N,                             n_mostra, descripció)
        (3 * 5,                              30,  "Trivial: 3 × 5 = 15"),
        (17 * 19,                            30,  "Petit: 17 × 19"),
        (101 * 103,                          30,  "Clàssic MDC: 101 × 103 = 10.403"),
        (1_000_003 * 1_000_033,              40,  "Mig: 7 dígits × 7 dígits"),
        (999_983 * 1_000_003,                60,  "Mig asimètric: 6 × 7 dígits"),
        (1_548_586_332_452_843,             100,  "Prova Colab (16 dígits)"),
        (1_000_000_000_039 * 100_000_000_3, 120,  "Gran: 13 × 10 dígits"),
    ]

    print("\n" + "█"*72)
    print("  🚀 BATERIA DE PROVES — MDC v15 (Aritmètica Exacta)")
    print("█"*72)

    resultats = []
    for (N, n_mostra, desc) in proves:
        print(f"\n{'─'*72}")
        print(f"  📋 {desc}")
        print(f"  N = {N:,}  ({len(str(N))} dígits)")
        print(f"  Mostra: {n_mostra} punts/segment")

        factor, t_ms = mdc_v15(N, n_mostra=n_mostra, verbose=verbose_global)

        if factor is not None:
            assert N % factor == 0, f"ERROR DE VERIFICACIÓ: {factor} no divideix {N}!"
            ok = "✅"
            info = f"{factor:,} × {N//factor:,}"
        else:
            ok = "⚠️ "
            info = "No trobat (augmenta n_mostra)"

        print(f"\n  {ok} Resultat : {info}")
        print(f"  ⏱️  Temps   : {t_ms:.3f} ms")
        resultats.append((desc, factor is not None, t_ms))

    # Resum final
    print(f"\n{'█'*72}")
    print(f"  RESUM FINAL:")
    n_ok   = sum(1 for r in resultats if r[1])
    n_fail = len(resultats) - n_ok
    print(f"  ✅ Encerts: {n_ok}/{len(resultats)}")
    if n_fail > 0:
        print(f"  ⚠️  Fallades: {n_fail}  (augmenta n_mostra en eixos casos)")
    print("█"*72)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

if __name__ == "__main__":
    # ── Prova ràpida individual ────────────────────────────────────────────
    N_PROVA = 10403          # 101 × 103  (el clàssic de les sessions anteriors)
    factor, t_ms = mdc_v15(N_PROVA, n_mostra=40, verbose=True)

    # ── Bateria completa ───────────────────────────────────────────────────
    # Descomenta la línia següent per executar totes les proves:
    # executar_bateria(verbose_global=False)
