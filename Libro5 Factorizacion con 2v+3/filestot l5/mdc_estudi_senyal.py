# -*- coding: utf-8 -*-
"""
==============================================================================
  MDC — ESTUDI VISUAL DE LA SENYAL I EL SEU ESPELL
  Pre-estudi per a v16: Detector de Desfase de Freqüència
==============================================================================
  Propòsit: Visualitzar com es veuen les dues senyals (B i A-espell),
  detectar els seus màxims i mínims, calcular cinemàtica (V, A, jerk)
  a partir de 4 pics consecutius, i mesurar el desfase de freqüència
  entre la senyal i el seu espell. Quan el signe del desfase canvia →
  hi ha un factor atrapat en eixe interval.

  Nomenclatura:
    N       : el nombre a factoritzar
    m       : paràmetre de cerca  (D = 2m+3)
    d_B(m)  : N mod D / D − 1/2         (senyal principal, divisor gran)
    m_esp   : (N // D − 3) // 2         (m del divisor espell)
    d_A(m)  : d_frac(m_esp)             (senyal espell, divisor menut)
    desfase : freqüència_B − freqüència_A a cada finestra de 4 pics
==============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.signal import argrelextrema

# ─────────────────────────────────────────────────────────────────────────────
#  PARÀMETRES D'ESTUDI  (canvia ací per explorar)
# ─────────────────────────────────────────────────────────────────────────────

# N petit i verificable: factors coneguts per poder comprovar la teoria
# 100003 × 100019  (dos primers propers → factors "quasi equilibrats")
N = 100003 * 100019       # = 10001900057

# Finestra de visualització al voltant de la convergència e-2
# m_max ≈ 50003  →  convergència ≈ 0.7183 × 50003 ≈ 35919
# m dels factors reals (necessaris per centrar la finestra)
m_factor_A = (100003 - 3) // 2    # 50000
m_factor_B = (100019 - 3) // 2    # 50008

M_MAX     = (int(N**0.5) - 3) // 2
M_CONV    = int((np.e - 2) * M_MAX)    # convergència e-2 (límit inferior de cerca)
FINESTRA  = 300                         # ± punts al voltant del centre triat

# Per factors EQUILIBRATS (tots dos ≈ √N): la zona interessant és prop de m_max
# Per factors ASIMÈTRICS: la zona és prop de m_conv (e-2 · m_max)
# Ací usem m_max com a centre (factors 100003 ≈ 100019 estan prop de m_max)
M_CENTRE  = m_factor_A                  # centres sobre el factor real conegut
M_INI = max(1, M_CENTRE - FINESTRA)
M_FI  = min(M_MAX, M_CENTRE + FINESTRA)

print(f"N = {N:,}  ({len(str(N))} dígits)")
print(f"Factors reals: 100003 × 100019")
print(f"m_max        = {M_MAX:,}")
print(f"m_conv (e-2) = {M_CONV:,}")
print(f"Finestra     : [{M_INI:,}, {M_FI:,}]")

# ─────────────────────────────────────────────────────────────────────────────
#  FUNCIONS DE LA SENYAL (float per velocitat en visualització)
# ─────────────────────────────────────────────────────────────────────────────

def d_B(m):
    """Senyal principal: distància de fase del divisor gran D=2m+3."""
    D = 2 * (2 * m + 3)
    return (N % D) / D - 0.5

def d_A(m):
    """
    Senyal espell: A·B = N, B = 2m+3  →  A = N // B  →  m_esp = (A-3)//2
    Mesura la fase del factor complementari (el factor menut si B és gran).
    """
    B = 2 * m + 3
    if B == 0:
        return 0.0
    A = N // B
    if A < 3:
        return 0.0
    m_esp = (A - 3) // 2
    if m_esp < 1:
        return 0.0
    D_esp = 2 * (2 * m_esp + 3)
    return (N % D_esp) / D_esp - 0.5

def m_espell(m):
    """Retorna m_espell corresponent a m (per marcar-lo a l'eix)."""
    B = 2 * m + 3
    A = N // B
    return (A - 3) // 2 if A >= 3 else None

# ─────────────────────────────────────────────────────────────────────────────
#  GENERACIÓ DE DADES
# ─────────────────────────────────────────────────────────────────────────────

ms      = np.arange(M_INI, M_FI + 1)
yB      = np.array([d_B(m) for m in ms])
yA      = np.array([d_A(m) for m in ms])
desfase = yB - yA     # diferència de fase B−A punt a punt

# ─────────────────────────────────────────────────────────────────────────────
#  DETECCIÓ DE MÀXIMS I MÍNIMS  (finestra de 3 punts)
# ─────────────────────────────────────────────────────────────────────────────

def trobar_extrems(y, ordre=2):
    """
    Retorna índexs dels màxims locals i mínims locals.
    'ordre' = quants veïns a cada costat es comparen.
    """
    idx_max = argrelextrema(y, np.greater, order=ordre)[0]
    idx_min = argrelextrema(y, np.less,    order=ordre)[0]
    return idx_max, idx_min

idx_max_B, idx_min_B = trobar_extrems(yB, ordre=3)
idx_max_A, idx_min_A = trobar_extrems(yA, ordre=3)

# ─────────────────────────────────────────────────────────────────────────────
#  CINEMÀTICA A PARTIR DE 4 PICS CONSECUTIUS
#  (velocitat, acceleració, jerk usant diferències finites)
# ─────────────────────────────────────────────────────────────────────────────

def cinematica_4punts(vals):
    """
    Donat un array de valors (altituds de pics/valls consecutius),
    calcula amb diferències finites:
      V  = primera diferència   (velocitat de canvi d'amplitud)
      A  = segona diferència    (acceleració)
      J  = tercera diferència   (jerk = variació de l'acceleració)

    Amb 4 punts p0,p1,p2,p3 obtenim 3 V, 2 A i 1 J.
    Retorna arrays de longitud 3, 2, 1 respectivament.
    """
    V = np.diff(vals)           # [p1-p0, p2-p1, p3-p2]
    A = np.diff(V)              # [V1-V0, V2-V1]
    J = np.diff(A)              # [A1-A0]  → el jerk
    return V, A, J

def freqüència_local(posicions):
    """
    Freqüència local = 1 / (distància mitja entre pics consecutius).
    'posicions' = array de m-valors dels pics.
    Retorna un array de freqüències (una per cada parell de pics).
    """
    if len(posicions) < 2:
        return np.array([])
    intervals = np.diff(posicions)   # distàncies entre pics consecutius
    return 1.0 / intervals           # freqüència = inversa del període

# Agafem fins als primers 12 màxims de cada senyal per la cinemàtica
N_PICS = min(12, len(idx_max_B), len(idx_max_A))

pos_pics_B = ms[idx_max_B[:N_PICS]]      # posicions m dels màxims de B
pos_pics_A = ms[idx_max_A[:N_PICS]]      # posicions m dels màxims de A

amp_pics_B = yB[idx_max_B[:N_PICS]]      # amplituds dels màxims de B
amp_pics_A = yA[idx_max_A[:N_PICS]]      # amplituds dels màxims de A

# Cinemàtica sobre amplituds dels màxims
if len(amp_pics_B) >= 4:
    V_B, A_B, J_B = cinematica_4punts(amp_pics_B[:4])
    print(f"\nCinemàtica Senyal B (4 màxims):")
    print(f"  Velocitat V : {V_B}")
    print(f"  Acceleració A: {A_B}")
    print(f"  Jerk J      : {J_B}")

if len(amp_pics_A) >= 4:
    V_A, A_A, J_A = cinematica_4punts(amp_pics_A[:4])
    print(f"\nCinemàtica Espell A (4 màxims):")
    print(f"  Velocitat V : {V_A}")
    print(f"  Acceleració A: {A_A}")
    print(f"  Jerk J      : {J_A}")

# Freqüències locals
freq_B = freqüència_local(pos_pics_B)
freq_A = freqüència_local(pos_pics_A)

# ─────────────────────────────────────────────────────────────────────────────
#  DESFASE DE FREQÜÈNCIA I CANVIS DE SIGNE
#  (la clau del detector de factors de v16)
# ─────────────────────────────────────────────────────────────────────────────

# Desfase punt a punt (yB - yA a cada m de la finestra)
signes = np.sign(desfase)

# Detecció de canvis de signe: on signes[i] ≠ signes[i+1]
canvis_signe = np.where(np.diff(signes) != 0)[0]
m_canvis     = ms[canvis_signe]    # posicions m dels canvis de signe

print(f"\nCanvis de signe del desfase (segments candidats):")
for mc in m_canvis:
    print(f"  m = {mc:,}   D = {2*mc+3:,}   N//D = {N//(2*mc+3):,}")

# ─────────────────────────────────────────────────────────────────────────────
#  VISUALITZACIÓ  (4 subgràfics)
# ─────────────────────────────────────────────────────────────────────────────

plt.style.use('dark_background')
fig = plt.figure(figsize=(16, 13), facecolor='#0d0d0d')
fig.suptitle(
    f'MDC — Anàlisi de Senyal i Espell  |  N = {N:,}  |  '
    f'Finestra m ∈ [{M_INI:,}, {M_FI:,}]  |  m_conv(e-2) = {M_CONV:,}',
    fontsize=12, color='#e0e0e0', y=0.98
)

gs = gridspec.GridSpec(4, 1, hspace=0.55,
                       left=0.07, right=0.97, top=0.93, bottom=0.05)

# Colors del tema
COL_B      = '#00aaff'    # blau   → senyal B (divisor gran)
COL_A      = '#ff6600'    # taronja → senyal A-espell (divisor menut)
COL_MAX    = '#00ff88'    # verd   → màxims
COL_MIN    = '#ff3366'    # roig   → mínims
COL_CANVI  = '#ffff00'    # groc   → canvi de signe (alerta de factor)
COL_FACTOR = '#ff00ff'    # magenta → posició factor real

def marca_factors(ax):
    """Afig línies verticals on estan els factors reals."""
    for mf, etiq in [(m_factor_A, 'p=100003'), (m_factor_B, 'q=100019')]:
        if M_INI <= mf <= M_FI:
            ax.axvline(mf, color=COL_FACTOR, lw=1.5, ls='--', alpha=0.8)
            ax.text(mf + 1, ax.get_ylim()[1] * 0.85,
                    etiq, color=COL_FACTOR, fontsize=7.5, va='top')

def marca_convergencia(ax):
    """Afig línia vertical a m_conv = (e-2)·m_max."""
    ax.axvline(M_CONV, color='#888888', lw=1, ls=':', alpha=0.6)
    ax.text(M_CONV + 1, ax.get_ylim()[1] * 0.95,
            f'm_conv={M_CONV:,}', color='#aaaaaa', fontsize=7)

# ── Subgràfic 1: Senyals B i A (ona principal + espell) ────────────────────
ax1 = fig.add_subplot(gs[0])
ax1.set_facecolor('#111111')
ax1.plot(ms, yB, color=COL_B,   lw=0.9, alpha=0.9, label='d_B(m)  [divisor gran D=2m+3]')
ax1.plot(ms, yA, color=COL_A,   lw=0.9, alpha=0.9, label='d_A(m_espell)  [espell A=N//B]')
ax1.axhline(0, color='white', lw=0.5, ls='--', alpha=0.4)

# Màxims i mínims de la senyal B
if len(idx_max_B):
    ax1.scatter(ms[idx_max_B], yB[idx_max_B],
                color=COL_MAX, s=18, zorder=5, label='Màxims B')
if len(idx_min_B):
    ax1.scatter(ms[idx_min_B], yB[idx_min_B],
                color=COL_MIN, s=18, zorder=5, label='Mínims B')

# Màxims i mínims de la senyal A-espell
if len(idx_max_A):
    ax1.scatter(ms[idx_max_A], yA[idx_max_A],
                color=COL_MAX, s=10, marker='^', zorder=5, alpha=0.6)
if len(idx_min_A):
    ax1.scatter(ms[idx_min_A], yA[idx_min_A],
                color=COL_MIN, s=10, marker='v', zorder=5, alpha=0.6)

ax1.set_ylabel('d(m) = N%D/D − ½', color='#cccccc', fontsize=9)
ax1.set_title('Senyal B (blau) i Espell A (taronja)  —  '
              'Màxims ● / Mínims ●', color='#dddddd', fontsize=10)
ax1.legend(fontsize=7.5, loc='upper right', facecolor='#222222',
           edgecolor='#444444', labelcolor='#dddddd')
ax1.tick_params(colors='#888888', labelsize=8)
ax1.spines[:].set_color('#333333')
ax1.set_xlim(M_INI, M_FI)
marca_convergencia(ax1)
marca_factors(ax1)

# ── Subgràfic 2: Desfase punt a punt (yB − yA) ────────────────────────────
ax2 = fig.add_subplot(gs[1], sharex=ax1)
ax2.set_facecolor('#111111')
ax2.plot(ms, desfase, color='#cc88ff', lw=0.9, alpha=0.95,
         label='Desfase(m) = d_B − d_A')
ax2.axhline(0, color='white', lw=0.8, ls='-', alpha=0.5)
ax2.fill_between(ms, desfase, 0,
                 where=(desfase > 0), color='#cc88ff', alpha=0.15)
ax2.fill_between(ms, desfase, 0,
                 where=(desfase < 0), color='#ff6644', alpha=0.15)

# Canvis de signe → segments candidats
for mc in m_canvis:
    ax2.axvline(mc, color=COL_CANVI, lw=1.2, alpha=0.8, zorder=6)

# Marcador especial: canvis que encerclen el factor real
for mc in m_canvis:
    if abs(mc - m_factor_A) < 20 or abs(mc - m_factor_B) < 20:
        ax2.axvline(mc, color=COL_CANVI, lw=2.5, alpha=1.0, zorder=7)
        ax2.text(mc + 1, ax2.get_ylim()[0] * 0.8,
                 '⚡', color=COL_CANVI, fontsize=11, zorder=8)

ax2.set_ylabel('d_B − d_A', color='#cccccc', fontsize=9)
ax2.set_title(f'Desfase Punt a Punt  |  {len(m_canvis)} canvis de signe '
              f'(línies ⚡ grogues = segments candidats)',
              color='#dddddd', fontsize=10)
ax2.legend(fontsize=7.5, loc='upper right', facecolor='#222222',
           edgecolor='#444444', labelcolor='#dddddd')
ax2.tick_params(colors='#888888', labelsize=8)
ax2.spines[:].set_color('#333333')
marca_convergencia(ax2)
marca_factors(ax2)

# ── Subgràfic 3: Freqüències locals de B i A ──────────────────────────────
ax3 = fig.add_subplot(gs[2], sharex=ax1)
ax3.set_facecolor('#111111')

if len(freq_B) > 0 and len(pos_pics_B) >= 2:
    m_freq_B = (pos_pics_B[:-1] + pos_pics_B[1:]) / 2   # punt mig entre pics
    ax3.stem(m_freq_B, freq_B, linefmt=COL_B+'88',
             markerfmt='C0o', basefmt='none', label='Freq B')

if len(freq_A) > 0 and len(pos_pics_A) >= 2:
    m_freq_A = (pos_pics_A[:-1] + pos_pics_A[1:]) / 2
    ax3.stem(m_freq_A, freq_A, linefmt=COL_A+'88',
             markerfmt='C1o', basefmt='none', label='Freq A')

# Desfase de freqüència: si tenim prou pics, interpolem i restem
if len(freq_B) >= 2 and len(freq_A) >= 2:
    # Interpolació al rang comú de m
    m_comuns = np.linspace(max(m_freq_B[0], m_freq_A[0]),
                           min(m_freq_B[-1], m_freq_A[-1]), 50)
    if len(m_freq_B) >= 2 and len(m_freq_A) >= 2:
        fB_interp = np.interp(m_comuns, m_freq_B, freq_B)
        fA_interp = np.interp(m_comuns, m_freq_A, freq_A)
        desfase_freq = fB_interp - fA_interp
        ax3_twin = ax3.twinx()
        ax3_twin.plot(m_comuns, desfase_freq, color='#ffee44',
                      lw=1.2, ls='--', alpha=0.7, label='Δfreq B−A')
        ax3_twin.axhline(0, color='#ffee44', lw=0.5, alpha=0.4)
        ax3_twin.set_ylabel('Δfreqüència', color='#ffee44', fontsize=8)
        ax3_twin.tick_params(colors='#ffee44', labelsize=7)
        ax3_twin.legend(fontsize=7.5, loc='lower right',
                        facecolor='#222222', edgecolor='#444444',
                        labelcolor='#dddddd')

ax3.set_ylabel('Freqüència local\n(1/interval entre pics)', color='#cccccc', fontsize=8)
ax3.set_title('Freqüència Local dels Pics  |  Desfase de Freqüència (línia groga)',
              color='#dddddd', fontsize=10)
ax3.legend(fontsize=7.5, loc='upper right', facecolor='#222222',
           edgecolor='#444444', labelcolor='#dddddd')
ax3.tick_params(colors='#888888', labelsize=8)
ax3.spines[:].set_color('#333333')
marca_factors(ax3)

# ── Subgràfic 4: Cinemàtica (V, A, J) dels màxims de B ───────────────────
ax4 = fig.add_subplot(gs[3], sharex=ax1)
ax4.set_facecolor('#111111')

if len(amp_pics_B) >= 4:
    # Calculem V, A, J per a finestres lliscants de 4 pics consecutius
    Vs, As, Js, mids = [], [], [], []
    for k in range(len(amp_pics_B) - 3):
        seg    = amp_pics_B[k:k+4]
        V_, A_, J_ = cinematica_4punts(seg)
        m_mid  = np.mean(pos_pics_B[k:k+4])
        Vs.append(V_[1])    # velocitat central del grup
        As.append(A_[0])    # acceleració del grup
        Js.append(J_[0])    # jerk del grup
        mids.append(m_mid)

    mids = np.array(mids)
    Vs   = np.array(Vs)
    As   = np.array(As)
    Js   = np.array(Js)

    ax4.plot(mids, Vs, color='#44ffaa', lw=1.2, marker='o',
             ms=5, label='V (velocitat d\'amplitud)')
    ax4.plot(mids, As, color='#ff8844', lw=1.2, marker='s',
             ms=5, label='A (acceleració)')
    ax4.plot(mids, Js, color='#ff44ff', lw=1.2, marker='^',
             ms=5, label='J (jerk)')
    ax4.axhline(0, color='white', lw=0.5, ls='--', alpha=0.4)

    # V−|V| (el terme directional del MDC, mencionat en el document)
    v_directional = Vs - np.abs(Vs)
    ax4.fill_between(mids, v_directional, 0,
                     color='#44ffaa', alpha=0.12,
                     label='V−|V| (component direccional)')

ax4.set_xlabel('m  (paràmetre de cerca,  D = 2m+3)', color='#cccccc', fontsize=9)
ax4.set_ylabel('Magnitud', color='#cccccc', fontsize=9)
ax4.set_title('Cinemàtica dels Màxims de B:  V · A · J  i  V−|V|  (directional)',
              color='#dddddd', fontsize=10)
ax4.legend(fontsize=7.5, loc='upper right', facecolor='#222222',
           edgecolor='#444444', labelcolor='#dddddd')
ax4.tick_params(colors='#888888', labelsize=8)
ax4.spines[:].set_color('#333333')
marca_convergencia(ax4)
marca_factors(ax4)

# ─────────────────────────────────────────────────────────────────────────────
#  GUARDAR I MOSTRAR
# ─────────────────────────────────────────────────────────────────────────────
plt.savefig('/mnt/user-data/outputs/mdc_estudi_senyal.png',
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
print(f"\nGràfic guardat: mdc_estudi_senyal.png")
plt.show()
print("Done.")
