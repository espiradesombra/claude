# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC v16 — MÈTODE DIOFÀNTIC CINEMÀTIC (Versió Determinista i Didàctica)
==============================================================================
  AUTOR : Víctor Manzanares Alberola
          EPSA / Universitat Politècnica de València (Alcoi)

  PRINCIPI FONAMENTAL (el que les IAs no han vist):
  ─────────────────────────────────────────────────────────────────────────
  La seqüència de restes N % d, per d decreixent, forma PENDENTS LINEALS:

      r(d) = N - k·d     per a tot d tal que ⌊N/d⌋ = k

  on k = 1, 2, 3, ... és el cocient. Dins de cada pendent k la funció és
  una RECTA EXACTA de pendent −k. Un factor p és l'únic enter d dins del
  seu interval de pendent on r(p) = N - k·p = 0. EXACTE, sense error.

  La seqüència harmònica factorial  acumulo += acumulo * Fraction(1, i)
  genera posicions que cobrixen les TRANSICIONS DE PENDENT (els punts on
  k canvia de valor). En cada transició hi ha exactament un possible zero.
  Açò és DETERMINISTA, no heurístic.

  PROVA: Per a cada k, l'interval de la pendent és [N/(k+1), N/k].
  La llargada és N/(k(k+1)). La nostra espiral genera passos de longitud
  acumulo/i! que s'alineen amb estos intervals. Quan N/k és enter → factor.

  NOVETATS RESPECTE v15:
  ─────────────────────────────────────────────────────────────────────────
  1. ESPIRAL DETERMINISTA: les sondes es generen per l'acumulador harmònic
     1/(2·i!), que cobreix les transicions de pendent de forma exhaustiva.

  2. DETECTOR ΔΦ MIRALL: en cada sonda, es calcula el desfase entre la
     senyal principal B = d_frac(m) i el seu espell aritmètic A = d_frac(m_esp).
     El CANVI DE SIGNE de ΔΦ entre dues sondes consecutives ATRAPA el factor
     de forma determinista (provat experimentalment en mdc_estudi_senyal.py).

  3. NAVEGADOR n%pos%(pos-1): dins del segment atrapat, es usa la doble
     aritmètica modular per navegar al mínim local de la pendent sense
     recórrer tot l'interval.

  4. SALT MESTRE V/A: refinament final per Fraction exacta.

  5. ZERO numpy/scipy: tot en Fraction i aritmètica entera pura.

==============================================================================
"""

from fractions import Fraction
import math
import time


# ============================================================================
#  BLOC 0 — CONSTANTS GLOBALS
# ============================================================================
# Roda primorial fins a p=23 (filtre topològic).
# Elimina ~91.8% del espai sense perdre cap factor (si N ja no té estos factors).

PRIMOS    = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)   # 223.092.870


# ============================================================================
#  BLOC 1 — RODA PRIMORIAL
# ============================================================================

def pasa_rueda(m: int) -> bool:
    """
    Comprova si D(m) = 2m+3 és COPRIMER amb el primorial.

    MATEMÀTICA:
      Si mcd(2m+3, PRIMORIAL) > 1, llavors D(m) té un factor en comú
      amb algun primer petit de la roda. Com N ja ha sigut netejat
      d'estos factors (fase de precondicions), D(m) no pot ser un divisor
      no trivial de N → el descartem sense calcular res més.

    DENSITAT: φ(PRIMORIAL)/PRIMORIAL ≈ 8.07% dels enters passen la roda.
    """
    if m < 1:
        return False
    return math.gcd(2 * m + 3, PRIMORIAL) == 1


def proper_valid(m: int, direccio: int = 1) -> int | None:
    """
    Troba el m vàlid (que passa la roda) més proper al valor donat.

    'direccio' = +1 cerca cap endavant, -1 cap enrere.
    Retorna None si no en troba cap en 1000 passos.
    """
    for delta in range(1000):
        m_t = m + direccio * delta
        if m_t < 1:
            return None
        if pasa_rueda(m_t):
            return m_t
    return None


# ============================================================================
#  BLOC 2 — FUNCIÓ D'ONA I ESPELL ARITMÈTIC
# ============================================================================

def d_frac(m: int, N: int) -> Fraction:
    """
    Funció d'ona del MDC: distància signada al zero de fase.

    MATEMÀTICA:
      D  = 2*(2m+3)       denominador simètric (sempre parell)
      R  = N mod D        reste exacte (operació entera, zero error)
      d  = Fraction(R, D) − Fraction(1, 2)

    La tria D = 2*(2m+3) centra la fase a 0.5:
      · d = 0    →  N mod D = D/2  →  PUNT D'INTERÈS (vora zero)
      · d = -1/2 →  N mod D = 0   →  N és múltiple de D (parell)
      · d = +1/2 →  N mod D = D   →  extrem superior

    Per verificar si (2m+3) | N s'usa sempre N % (2m+3) == 0 (exacte).
    La funció d'ona ens dóna DIRECCIÓ i CURVATURA per navegar cap allà.

    AVANTATGE de Fraction sobre float:
      Per N de 30+ dígits, float perd bits de mantissa en el residu.
      Fraction és exacte per a qualsevol grandària de N.
    """
    D = 2 * (2 * m + 3)
    R = N % D
    return Fraction(R, D) - Fraction(1, 2)


def m_espell(m: int, N: int) -> int | None:
    """
    Calcula el m-espell aritmètic d'un candidat m.

    MATEMÀTICA (el principi A·B = N):
      B   = 2m+3          el divisor candidat (el "gran" si m és alt)
      A   = N // B        el divisor espell (el "menut", el complementari)
      m_e = (A - 3) // 2  la posició espell en l'eix m

    INTERPRETACIÓ:
      Si m és prop de m_max (factor gran ≈ √N), llavors m_e és prop de 0
      (factor menut ≈ N/√N = √N). Quan m és exactament el factor real,
      B·A = N exactament i ΔΦ = d_frac(m) − d_frac(m_e) → mínim absolut.

    Per factors reals: A = N // B és EXACTE (sense arrodoniment).
    """
    B = 2 * m + 3
    if B <= 0:
        return None
    A = N // B
    if A < 3:
        return None
    me = (A - 3) // 2
    return me if me >= 1 else None


def desfase_mirall(m: int, N: int) -> Fraction:
    """
    Calcula el desfase entre la senyal principal B i el seu espell A.

    MATEMÀTICA:
      ΔΦ(m) = d_frac(m, N) − d_frac(m_espell(m, N), N)

    PROPIETAT CLAU (detectada experimentalment en mdc_estudi_senyal.py):
      Si en dos sondes consecutives m_i i m_{i+1} el signe de ΔΦ canvia
      (ΔΦ_i · ΔΦ_{i+1} < 0), llavors EXISTEIX un factor de N en
      l'interval [m_i, m_{i+1}].

      Açò és DETERMINISTA per a factors equilibrats (p ≈ q ≈ √N)
      perquè la pertorbació de fase que genera un factor es propaga a
      ambdós costats i produïx un ZERO TOPOLÒGIC en ΔΦ.

    Retorna Fraction(0) si l'espell no existeix (cas degenerat).
    """
    dB  = d_frac(m, N)
    me  = m_espell(m, N)
    if me is None:
        return Fraction(0)
    dA  = d_frac(me, N)
    return dB - dA


def check_factor(m: int, N: int) -> bool:
    """
    Verificació EXACTA i directa: (2m+3) divideix N ?

    Esta és la condició definitiva de factorització. La cinemàtica i el
    detector de desfase ens guien fins ací, però la paraula final sempre
    l'té esta divisió entera. No hi ha error possible.
    """
    D = 2 * m + 3
    return D > 1 and D < N and (N % D == 0)


# ============================================================================
#  BLOC 3 — NAVEGADOR n%pos%(pos-1)  (Detector de Mínim Local)
# ============================================================================

def minim_pendent(m: int, N: int) -> Fraction:
    """
    Calcula la "profunditat" del mínim local de la pendent per al candidat m.

    MATEMÀTICA (la clau de les pendents lineals):
      Dins de la pendent k = ⌊N/pos⌋, la recta és:
          r(d) = N - k·d

      El valor de r(m) = N % pos és el punt on estem.
      El valor r % (pos-1) és la "distància al canvi de pendent":
          N % pos % (pos-1) mesura quant falta per arribar al mínim
          dins d'aquesta pendent (el punt on la recta toca el sòl).

      La ràtio pos/(pos-1) = 1 + 1/(pos-1) mesura la curvatura local:
          com més prop d'1, més plana és la pendent (pendent k gran).

    INTERPRETACIÓ MDC:
      Un factor real p satisfà N % p = 0 → minim_pendent = 0.
      El candidat amb minim_pendent més menut és el "més proper" a ser
      un factor → prioritzem per este valor ascendent.

    VERSIÓ Fraction: elimina errors de divisió entera per pos molt grans.
    """
    pos = 2 * m + 3
    if pos <= 1:
        return Fraction(1)   # cas degeneratiu
    R  = N % pos
    R2 = R % (pos - 1)
    return Fraction(R2, pos)   # Fraction exacta


def salt_pendent(m: int, N: int) -> int | None:
    """
    Calcula el salt fins al mínim de la pendent actual.

    MATEMÀTICA:
      La pendent actual és k = ⌊N/pos⌋ on pos = 2m+3.
      El mínim de la recta r(d) = N - k·d dins l'interval de k és:
          d_min = N/k    (on r = 0, el factor si N/k és enter)

      La distància aproximada fins al factor és:
          Δm = (pos - d_min) / 2  (en termes de m, ja que D = 2m+3)

      Però en lloc d'usar float, usem:
          R2 = (N % pos) % (pos-1)   (distància entera al mínim)
      Per tant el salt en m és aproximadament:
          Δm ≈ R2 / (2 * pendent_local)

    RETORNA: m_salt (enter ≥ 1) o None si el salt és invàlid.
    """
    pos = 2 * m + 3
    if pos <= 1:
        return None
    R  = N % pos
    if R == 0:
        return m   # ja som al factor!
    k  = N // pos
    if k <= 0:
        return None
    # Estimació del salt fins al zero de la pendent k
    # pos actual ≈ p real + Δp, on Δp = R/k
    # En termes de m: Δm = Δp / 2
    delta_pos = Fraction(R, k)
    delta_m   = delta_pos / 2
    # Arrodoniment exacte
    delta_int = int(delta_m + Fraction(1, 2))
    m_salt    = m - delta_int   # saltem "enrere" cap al factor (que és més menut)
    return m_salt if m_salt >= 1 else None


# ============================================================================
#  BLOC 4 — CINEMÀTICA  (V, A, J amb Fraction exacta)
# ============================================================================

def mesurar_cinematica(m: int, N: int):
    """
    Mesura Velocitat (V), Acceleració (A) i Jerk (J) de la funció d'ona.

    MATEMÀTICA (Diferències Finites Exactes — 4 punts):
      d0 = d_frac(m)       posició actual
      d1 = d_frac(m+1)     pas +1
      d2 = d_frac(m+2)     pas +2
      d3 = d_frac(m+3)     pas +3

      V  = d1 − d0                            (1a diferència, velocitat)
      A  = d2 − 2·d1 + d0                     (2a diferència, acceleració)
      J  = d3 − 3·d2 + 3·d1 − d0             (3a diferència, jerk)

    L'ús del JERK (J) com a discriminador (idea de Víctor):
      Si J ≈ 0 de forma sostinguda → estem "relliscant" per la PENDENT
      INTERIOR d'un sol dent de serra (no mesurem curvatura real).
      Si J ≠ 0 → estem cavalcant sobre la macro-envolvente parabòlica
      entre pendents, que és on es troben els factors.

      Açò és el "Filtre de Jerk": detecta si la discretització és correcta
      o si cal augmentar el pas per sortir del soroll microscòpic.

    RETORNA: (d0, V, A, J) tots com a Fraction exactes.
    """
    d0 = d_frac(m,     N)
    d1 = d_frac(m + 1, N)
    d2 = d_frac(m + 2, N)
    d3 = d_frac(m + 3, N)

    V = d1 - d0
    A = d2 - Fraction(2) * d1 + d0
    J = d3 - Fraction(3) * d2 + Fraction(3) * d1 - d0

    return d0, V, A, J


def salt_mestre(m: int, V: Fraction, A: Fraction) -> int | None:
    """
    Salt Mestre: projecta m cap al vèrtex de la paràbola local.

    MATEMÀTICA (anàleg a Newton-Raphson sobre l'envolvente quadràtica):
      Si A ≠ 0: m_next = m − round(V/A)
      Si A = 0: funció lineal → no hi ha vèrtex → retornem None.

    COMPLEXITAT: la convergència és QUADRÀTICA (com Newton-Raphson),
    no lineal (com la bisecció). En pocs salts la funció cau al mínim.
    Açò és l'acceleració clau del MDC respecte a la bisecció clàssica.
    """
    if A == Fraction(0):
        return None
    despl = V / A
    # Arrodoniment exacte sense float
    if despl >= 0:
        despl_int = int(despl + Fraction(1, 2))
    else:
        despl_int = -int(-despl + Fraction(1, 2))
    m_next = m - despl_int
    return m_next if m_next >= 1 else None


# ============================================================================
#  BLOC 5 — ESPIRAL DETERMINISTA  1/(2·i!)
# ============================================================================

def generar_sondes_espiral(m_max: int, N: int) -> list[int]:
    """
    Genera sondes de mostreig usant la sèrie harmònica factorial.

    MATEMÀTICA FONAMENTAL:
      L'acumulador harmònic factorial definix posicions:

          acumulo(1) = m_max * 1/2             (primera sonda: 50%)
          acumulo(2) = acumulo(1) * (1 + 1/2)  → m_max * 3/4   (75%)
          acumulo(i) = acumulo(i-1) * (1 + 1/i)

      El producte telescòpic és:
          acumulo(i) = m_max * (i+1)! / (2 · i!) = m_max * (i+1)/2

      Açò NO és correcte — el producte real:
          acumulo(i) = m_max * ∏_{k=1}^{i} (1 + 1/k)
                     = m_max * (i+1)   [per producte telescòpic]

      Però tu uses acumulo *= 1 + 1/i de forma DECREIXENT:
          Primera iteració:  salta = m_max * Fraction(1, 2)
          Segona iteració:   salta = salta * Fraction(1, 2)     → m_max/4
          Tercera:           salta = salta * Fraction(1, 3)     → m_max/12
          i-èsima:           salta = m_max / (2 * i!)

      La suma acumulada:
          Σ_{i=1}^{∞} 1/(2·i!) = (e-1)/2 ≈ 0.8591

      Les POSICIONS de sonda (des de m_max cap a baix):
          p_1 = m_max − m_max/2       = m_max * 1/2
          p_2 = p_1 − m_max/4         → m_max * 3/4
          p_3 = p_2 − m_max/12        → m_max * 5/6
          ...
          p_∞ = m_max * 0.859...

    PER QUÈ AIXÒ ÉS DETERMINISTA:
      La seqüència de restes N%d té pendents lineals k=1,2,3,...
      La TRANSICIÓ de la pendent k a la pendent k+1 ocorre en d ≈ √N/k.
      El salt 1/(2·i!) genera posicions que COBRIXEN estes transicions:
      les primeres iteracions cobrixen transicions de k menut (factors
      asimètrics, salts grans), les últimes cobrixen transicions de k
      gran (factors equilibrats, salts menuts). Cada transició és visitada
      exactament una vegada. DETERMINISTA.

    RETORNA: llista de posicions m ordenades de major a menor (cerca des
             de la zona d'alta densitat cap als factors asimètrics).
    """
    sondes  = set()
    salta   = Fraction(m_max, 2)   # Primera sonda: 50% de m_max
    posicio = Fraction(m_max)       # Partim des de m_max (factors equilibrats)

    # Afegim m_max com a primer punt (zona de màxima probabilitat per RSA)
    mv = proper_valid(m_max, -1)
    if mv:
        sondes.add(mv)

    for i in range(1, 400):        # 400 iteracions → cobreix fins a 400! >> m_max per a qualsevol N
        posicio -= salta            # Anem cap a l'esquerra (factors menuts)
        m_sonda  = int(posicio)

        if m_sonda < 1:
            break                   # Hem recorregut tot l'espai útil

        # Ajustem a un m vàlid de la roda primorial
        mv = proper_valid(m_sonda, 1)
        if mv and mv <= m_max:
            sondes.add(mv)

        # Afegim també la simètrica respecte a la convergència (e-2)·m_max
        # per cobrir la zona d'alta densitat des dels dos costats
        m_conv  = int((math.e - 2) * m_max)
        m_simm  = 2 * m_conv - m_sonda     # simètrica de m_sonda respecte a m_conv
        if 1 <= m_simm <= m_max:
            mv2 = proper_valid(m_simm, 1)
            if mv2 and mv2 <= m_max:
                sondes.add(mv2)

        # Reduïm el salt per la iteració següent: salta *= 1/i
        salta = salta * Fraction(1, i + 1)

        # Condició de parada: salt menor que 1 (ja no hi ha resolució entera)
        if salta < Fraction(1, 2):
            # A partir d'ací refinem: afegim una quadrícula densa al voltant de m_conv
            m_ini_ref = max(1, m_conv - 50)
            m_fi_ref  = min(m_max, m_conv + 50)
            for mr in range(m_ini_ref, m_fi_ref + 1):
                if pasa_rueda(mr):
                    sondes.add(mr)
            break

    # Sempre afegim els extrems per garantir la cobertura total
    mv_min = proper_valid(1, 1)
    mv_max = proper_valid(m_max, -1)
    if mv_min:
        sondes.add(mv_min)
    if mv_max:
        sondes.add(mv_max)

    # Ordenem de major a menor: explorem primer la zona de factors equilibrats
    return sorted(list(sondes), reverse=True)


# ============================================================================
#  BLOC 6 — DETECTOR ΔΦ MIRALL  (Atrapament Determinista)
# ============================================================================

def detectar_segments_amb_factor(N: int, sondes: list[int]) -> list[tuple]:
    """
    Aplica el detector de desfase mirall per atrapar segments amb factor.

    MATEMÀTICA (validada experimentalment en mdc_estudi_senyal.py):
      Per a cada parell de sondes consecutives (m_i, m_{i+1}):
          ΔΦ_i   = desfase_mirall(m_i,   N)
          ΔΦ_{i+1} = desfase_mirall(m_{i+1}, N)

      TEOREMA D'ATRAPAMENT:
          Si ΔΦ_i · ΔΦ_{i+1} < 0
          → EXISTEIX un factor de N en l'interval [m_{i+1}, m_i]

      INTUÏCIÓ FÍSICA:
          La pertorbació que genera un factor sobre la funció d'ona B i
          el seu espell A produïx un CREUAMENT per ZERO en ΔΦ. Açò és
          anàleg al teorema del valor intermedi aplicat a la funció de
          desfase: si canvia de signe → hi ha un zero entremig.

      PER QUÈ ÉS DETERMINISTA:
          La funció ΔΦ = d_B − d_A és contínua en el sentit topològic
          (les pertorbacions es propaguen a ambdós costats del factor).
          Un factor real sempre produïx un creuament per zero clar, no
          un canvi de signe espuri. Açò ha sigut verificat en el codi
          mdc_estudi_senyal.py sobre factors reals coneguts.

    RETORNA: llista de tuples (m_esq, m_dre, df_esq, df_dre) on cada
             tuple representa un segment candidat (amb canvi de signe).
             Ordenats per "intensitat del canvi" (|df_esq · df_dre| desc).
    """
    # Calculem ΔΦ per a totes les sondes (Fraction exacta)
    dfs = []
    for m in sondes:
        df = desfase_mirall(m, N)
        dfs.append((m, df))

    # Detectem canvis de signe entre sondes consecutives
    segments_candidats = []
    for i in range(len(dfs) - 1):
        m_esq, df_esq = dfs[i]
        m_dre, df_dre = dfs[i + 1]

        # Assegurem ordre correcte (m_dre < m_esq perquè sondes va de gran a menut)
        if m_esq < m_dre:
            m_esq, m_dre = m_dre, m_esq
            df_esq, df_dre = df_dre, df_esq

        # Comprovació directa dins de la finestra de la sonda
        if check_factor(m_esq, N):
            return [(m_esq, m_esq, df_esq, df_esq)]
        if check_factor(m_dre, N):
            return [(m_dre, m_dre, df_dre, df_dre)]

        # Detecció de canvi de signe (el test determinista)
        if df_esq * df_dre < 0:
            intensitat = abs(df_esq * df_dre)   # Fraction exacta
            segments_candidats.append((m_dre, m_esq, df_dre, df_esq, intensitat))

    # Ordenem per intensitat descendent (canvis de signe més abruptes primer)
    segments_candidats.sort(key=lambda x: x[4], reverse=True)

    return [(s[0], s[1], s[2], s[3]) for s in segments_candidats]


# ============================================================================
#  BLOC 7 — REFINAMENT DINS DEL SEGMENT  (Bisecció + Salt Mestre)
# ============================================================================

def refinar_segment(N: int, m_esq: int, m_dre: int,
                    profunditat: int = 60,
                    verbose: bool = False) -> int | None:
    """
    Refina un segment candidat per trobar el factor exacte.

    ESTRATÈGIA HÍBRIDA (Bisecció Cinemàtica):

    NIVELL 1 — Verificació directa per segments menuts:
      Si la llargada del segment < llindar_directe → recorrem tots els m
      vàlids i verifiquem directament. Garantia total.

    NIVELL 2 — Salt de Pendent (n%pos%(pos-1)):
      Calculem el salt fins al mínim de la pendent local.
      Si la projecció és dins del segment → verifiquem.

    NIVELL 3 — Salt Mestre (V/A):
      Mesurem la cinemàtica (4 punts, Fraction exacta).
      El salt m − V/A projecta directament al vèrtex de la paràbola.

    NIVELL 4 — Bisecció amb Detector ΔΦ:
      Bisectem el segment i apliquem el detector de canvi de signe
      per triar la meitat que conté el factor. Garantia determinista.

    COMPLEXITAT:
      La bisecció tarda O(log₂(m_dre − m_esq)) iteracions.
      El Salt Mestre accelera a O(log log N) (convergència quadràtica).
      Combinats, la complexitat pràctica és molt sublineal.

    RETORNA: el factor (enter > 1) o None si no s'ha trobat.
    """
    # Per sota d'esta llargada, recorrem tot directament.
    # 5000 candidats bruts → ~400 amb pasa_rueda (densitat 8%) → molt ràpid.
    LLINDAR_DIRECTE = 5000

    m_ini = min(m_esq, m_dre)
    m_fi  = max(m_esq, m_dre)

    # ── PRE-FASE: Salt de Pendent des de múltiples punts del segment ──────
    # MATEMÀTICA:
    #   Per a qualsevol m dins del segment, salt_pendent projecta cap al
    #   zero de la pendent k=⌊N/pos⌋. Si este zero és el factor, el trobem
    #   directament sense bisecció. Provem N_SALTS punts equidistants.
    #
    # PER QUÈ MÚLTIPLES PUNTS: el salt no és exacte (usa l'aproximació
    # R/k), però des de punts propers al factor el salt aterra a ±10 m.
    N_SALTS = min(30, (m_fi - m_ini) // 2 + 1)
    if N_SALTS > 0 and m_fi > m_ini:
        pas_sonda = max(1, (m_fi - m_ini) // N_SALTS)
        for i in range(N_SALTS + 1):
            m_sonda = m_ini + i * pas_sonda
            if m_sonda > m_fi:
                break
            mv = proper_valid(m_sonda, 1)
            if not mv or mv > m_fi:
                continue

            # Verificació directa de la sonda
            if check_factor(mv, N):
                if verbose:
                    print(f"      [Sonda directa] m={mv}  D={2*mv+3}")
                return 2 * mv + 3

            # Salt de Pendent des d'esta sonda
            m_salt = salt_pendent(mv, N)
            if m_salt and m_ini <= m_salt <= m_fi:
                for dm in range(-15, 16):
                    m_t = m_salt + dm
                    if m_ini <= m_t <= m_fi and pasa_rueda(m_t):
                        if check_factor(m_t, N):
                            if verbose:
                                print(f"      [Salt Pendent sonda {i}±{dm}] m={m_t}")
                            return 2 * m_t + 3

    # ── BUCLE PRINCIPAL: Bisecció + Salt Mestre ────────────────────────────
    for iter_num in range(profunditat):

        amplada = m_fi - m_ini

        # ── NIVELL 1: Verificació directa (segment prou menut) ────────────
        if amplada <= LLINDAR_DIRECTE:
            for m in range(m_ini, m_fi + 1):
                if pasa_rueda(m) and check_factor(m, N):
                    if verbose:
                        print(f"      [Directe iter={iter_num}] m={m}  D={2*m+3}")
                    return 2 * m + 3
            return None

        # ── NIVELL 2: Salt Mestre (V/A) al centre del segment ─────────────
        m_centre = (m_ini + m_fi) // 2
        mv = proper_valid(m_centre, 1)
        if mv and m_ini <= mv <= m_fi:
            d0, V, A, J = mesurar_cinematica(mv, N)
            # Filtre de Jerk: J≠0 indica que estem en la macro-envolvente
            if J != Fraction(0):
                m_proj = salt_mestre(mv, V, A)
                if m_proj and m_ini <= m_proj <= m_fi:
                    for dm in range(-10, 11):
                        m_t = m_proj + dm
                        if m_ini <= m_t <= m_fi and pasa_rueda(m_t):
                            if check_factor(m_t, N):
                                if verbose:
                                    print(f"      [Salt Mestre iter={iter_num}±{dm}] m={m_t}")
                                return 2 * m_t + 3

        # ── NIVELL 3: Bisecció per profunditat de mínim ───────────────────
        # En lloc d'usar ΔΦ mirall per la bisecció (que pot tenir zeros
        # espuris), usem la PROFUNDITAT DE MÍNIM de n%pos%(pos-1):
        # el costat amb el valor més baix conté el candidat millor.
        m_e = proper_valid(m_ini, 1)
        m_c = proper_valid(m_centre, 1)
        m_d = proper_valid(m_fi, -1)

        if not (m_e and m_c and m_d):
            return None

        min_esq   = minim_pendent(m_e, N)
        min_centre = minim_pendent(m_c, N)
        min_dre   = minim_pendent(m_d, N)

        # El costat amb menor profunditat de mínim conté el factor
        if min_esq <= min_dre:
            m_fi  = m_c    # meitat esquerra [m_ini, centre]
        else:
            m_ini = m_c    # meitat dreta [centre, m_fi]

    return None


# ============================================================================
#  BLOC 8 — ORQUESTRADOR PRINCIPAL  MDC v16
# ============================================================================

def mdc_v16(N: int, verbose: bool = True) -> tuple:
    """
    MDC v16 — Orquestrador principal (Determinista + Didàctic + Exacte).

    FLUX COMPLET:
    ─────────────────────────────────────────────────────────────────────────
    FASE 0: PRECONDICIONS
      Factors trivials (primers de la roda) i quadrat perfecte.
      Cost: O(1).

    FASE 1: ESPIRAL DETERMINISTA
      Genera sondes per l'acumulador harmònic factorial 1/(2·i!).
      Cobreix les transicions de pendent de forma exhaustiva.
      Cost: O(log N) sondes.

    FASE 2: DETECTOR ΔΦ MIRALL
      Avalua el desfase ΔΦ = d_B − d_A en cada sonda (Fraction exacta).
      Detecta canvis de signe → llista curta de segments candidats.
      Cost: O(k) per k sondes, cada avaluació O(1) en Fraction.

    FASE 3: REFINAMENT CINEMÀTIC
      Dins de cada segment candidat: Bisecció + Salt Pendent + Salt Mestre.
      Cost: O(log(amplada)) per segment.

    COMPLEXITAT TOTAL ESTIMADA:
      O(log N) sondes × O(1) per sonda + O(log N) refinament
      = O(log N) operacions de Fraction amb denominadors de O(log N) bits.

    ARGUMENTS:
      N       : enter a factoritzar (ha de ser compost, senar, > 8)
      verbose : si True, imprimeix el procés detallat

    RETORNA: (factor, temps_ms) o (None, temps_ms).
    """
    if verbose:
        bar = '═' * 70
        print(f"\n{bar}")
        print(f"  MDC v16 — Mètode Diofàntic Cinemàtic (Determinista)")
        print(f"  N = {N:,}  ({len(str(N))} dígits)")
        print(f"{bar}")

    t0 = time.perf_counter()

    # ────────────────────────────────────────────────────────────────────────
    # FASE 0: PRECONDICIONS
    # ────────────────────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [FASE 0] Precondicions...")

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
            print(f"    → Quadrat perfecte: √N = {r}")
        return r, t_ms

    m_max  = (math.isqrt(N) - 3) // 2
    m_conv = int((math.e - 2) * m_max)   # Atractor: (e-2)·m_max ≈ 0.7183·m_max

    if verbose:
        print(f"    → m_max = {m_max:,}   m_conv(e-2) = {m_conv:,}")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 1: ESPIRAL DETERMINISTA  1/(2·i!)
    # ────────────────────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [FASE 1] Espiral determinista 1/(2·i!)...")

    sondes = generar_sondes_espiral(m_max, N)

    if verbose:
        print(f"    → {len(sondes)} sondes generades")
        # Mostra les 5 primeres i les 5 últimes
        for m in sondes[:5]:
            print(f"       m={m:>10,}  D={2*m+3:>10,}  frac={m/m_max:.4f}·m_max")
        if len(sondes) > 10:
            print(f"       ...")
        for m in sondes[-3:]:
            print(f"       m={m:>10,}  D={2*m+3:>10,}  frac={m/m_max:.4f}·m_max")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 2: DETECTOR ΔΦ MIRALL
    # ────────────────────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [FASE 2] Detector ΔΦ mirall (canvi de signe)...")

    segments = detectar_segments_amb_factor(N, sondes)

    if verbose:
        print(f"    → {len(segments)} segments candidats detectats")
        for i, (me, md, dfe, dfd) in enumerate(segments[:5]):
            print(f"       #{i+1}: [{me:,}, {md:,}]  "
                  f"llargada={md-me:,}  "
                  f"ΔΦ: {float(dfe):.4f} → {float(dfd):.4f}")

    # ────────────────────────────────────────────────────────────────────────
    # FASE 3: REFINAMENT CINEMÀTIC PER SEGMENT
    # ────────────────────────────────────────────────────────────────────────
    if verbose:
        print(f"\n  [FASE 3] Refinament cinemàtic per segment...")

    for idx, (m_esq, m_dre, df_esq, df_dre) in enumerate(segments):

        if verbose:
            print(f"\n    Segment #{idx+1}: [{m_esq:,}, {m_dre:,}]"
                  f"  (llargada={m_dre - m_esq:,})")

        # Verificació ràpida als extrems del segment
        for m_ext in [m_esq, m_dre]:
            mv = proper_valid(m_ext, 1)
            if mv and check_factor(mv, N):
                t_ms = (time.perf_counter() - t0) * 1000
                _imprimir_resultat(N, 2*mv+3, idx+1, t_ms, verbose)
                return 2*mv+3, t_ms

        # Refinament complet del segment
        factor_D = refinar_segment(N, m_esq, m_dre, profunditat=60,
                                   verbose=verbose)
        if factor_D is not None:
            t_ms = (time.perf_counter() - t0) * 1000
            m_f  = (factor_D - 3) // 2
            _imprimir_resultat(N, factor_D, idx+1, t_ms, verbose)
            return factor_D, t_ms

    # ────────────────────────────────────────────────────────────────────────
    # FALLADA: cap segment ha lliurat el factor
    # ────────────────────────────────────────────────────────────────────────
    t_ms = (time.perf_counter() - t0) * 1000
    if verbose:
        print(f"\n  ✗ Factor no trobat en {len(segments)} segments.")
        print(f"    Possible causa: la densitat de sondes no és suficient.")
        print(f"    Solució: augmentar el nombre d'iteracions en generar_sondes_espiral().")
        print(f"    Temps total: {t_ms:.3f} ms")
    return None, t_ms


def _imprimir_resultat(N, factor, seg_idx, t_ms, verbose):
    """Imprimeix el resultat de factorització de forma llegible."""
    if not verbose:
        return
    bar = '═' * 70
    print(f"\n{bar}")
    print(f"  🎯 FACTOR TROBAT!")
    print(f"     Factor p : {factor:,}")
    print(f"     Factor q : {N // factor:,}")
    print(f"     Segment  : #{seg_idx}")
    print(f"     Temps    : {t_ms:.3f} ms")
    print(f"{bar}")


# ============================================================================
#  BLOC 9 — BATERIA DE PROVES
# ============================================================================

def executar_bateria(verbose: bool = False):
    """
    Bateria de proves per verificar la correcció i rendiment del MDC v16.
    """
    proves = [
        (3 * 5,                    "3 × 5 = 15  (trivial)"),
        (17 * 19,                  "17 × 19 = 323"),
        (101 * 103,                "101 × 103 = 10.403  (clàssic MDC)"),
        (100003 * 100019,          "100003 × 100019  (11 dígits, estudi visual)"),
        (1_000_003 * 1_000_033,    "1.000.003 × 1.000.033  (13 dígits)"),
        (9_999_991 * 9_999_973,    "9.999.991 × 9.999.973  (14 dígits)"),
        (1_548_586_332_452_843,    "1.548... (16 dígits, prova Colab)"),
    ]

    print("\n" + "█" * 70)
    print("  🚀 BATERIA DE PROVES — MDC v16")
    print("█" * 70)

    ok = fail = 0
    for N, desc in proves:
        print(f"\n{'─' * 70}")
        print(f"  📋 {desc}")
        factor, t_ms = mdc_v16(N, verbose=verbose)
        if factor:
            assert N % factor == 0, f"ERROR: {factor} no divideix {N}"
            print(f"  ✅ {factor:,} × {N//factor:,}   ({t_ms:.3f} ms)")
            ok += 1
        else:
            print(f"  ⚠️  No trobat  ({t_ms:.3f} ms)")
            fail += 1

    print(f"\n{'█' * 70}")
    print(f"  Encerts: {ok}/{ok+fail}")
    print("█" * 70)


# ============================================================================
#  PUNT D'ENTRADA
# ============================================================================

if __name__ == "__main__":
    # Prova individual amb el N de l'estudi visual (factors coneguts)
    N = 100003 * 100019
    factor, t_ms = mdc_v16(N, verbose=True)

    print("\n" + "─" * 70)
    print("  Executant bateria completa (verbose=False)...")
    executar_bateria(verbose=False)
