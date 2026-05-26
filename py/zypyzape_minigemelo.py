"""
╔══════════════════════════════════════════════════════════════════╗
║   ZYPYZAPE — MiniGemelo Digital v1.0                            ║
║   Simulación dinámica de 5 aerogeneradores con batería cinética  ║
║   Víctor Manzanares Alberola · UPV-EPSA Alcoy · 2026            ║
╚══════════════════════════════════════════════════════════════════╝

Requisitos: Python 3.8+, matplotlib, numpy
Ejecutar:   python3 zypyzape_minigemelo.py
Controles:
  - Slider VIENTO   → velocidad media del viento [m/s]
  - Slider REPARTO  → % turbinas en captación vs batería
  - Botón PERT.     → inyecta perturbación de ±200 MW en red
  - Botón PAUSA     → congela la simulación
  - Botón RESET     → reinicia el sistema
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════
# FÍSICA DEL SISTEMA
# ═══════════════════════════════════════════════════════════════════

# Parámetros físicos
J        = 5.0e6    # Momento de inercia del rotor [kg·m²]
R        = 40.0     # Radio del rotor [m]
RHO      = 1.225    # Densidad del aire [kg/m³]
OMEGA_NOM = 1.57    # Velocidad angular nominal [rad/s]
OMEGA_MIN = 1.20    # Límite inferior seguro [rad/s]
OMEGA_MAX = 1.85    # Límite superior seguro [rad/s]
P_NOM    = 2.5e6    # Potencia nominal por turbina [W]
F0       = 50.0     # Frecuencia nominal de red [Hz]
H_SIS    = 2.5      # Constante de inercia del sistema [s]
S_TOTAL  = 32.0e9   # Potencia total sistema [W]
DT       = 0.05     # Paso de simulación [s]
N_TURB   = 5        # Número de turbinas

ROLES = {
    'captacion':  {'label': 'CAPT',  'color': '#02C39A'},
    'bat_acelera':{'label': 'ACEL',  'color': '#F4A261'},
    'bat_frena':  {'label': 'FREN',  'color': '#E74C3C'},
}

def Cp_curva(lam) -> float:
    """Curva Cp(λ) paramétrica — perfil HAWT típico 2.5 MW."""
    lam = np.clip(float(lam), 2.0, 14.0)
    lam0 = 8.0
    return float(0.48 * np.exp(-0.5 * ((lam - lam0) / 2.5) ** 2))

def T_aerodinamico(v: float, omega: float) -> float:
    """Par aerodinámico [N·m] a partir del viento y velocidad angular."""
    if omega < 0.1 or v < 1.0:
        return 0.0
    lam = omega * R / v
    cp  = Cp_curva(lam)
    P   = 0.5 * RHO * np.pi * R**2 * v**3 * cp
    return P / omega

def T_MPPT(omega: float) -> float:
    """Par óptimo MPPT: T = K_opt · ω² """
    K_opt = 0.5 * RHO * np.pi * R**5 * 0.48 / (8.0**3)
    return K_opt * omega**2

# ═══════════════════════════════════════════════════════════════════
# ESTADO DEL SISTEMA
# ═══════════════════════════════════════════════════════════════════

class SistemaTurbinas:
    def __init__(self, n=N_TURB):
        self.n = n
        self.nombres = [f'T{i+1}' for i in range(n)]
        self.omega   = np.full(n, OMEGA_NOM) + np.random.uniform(-0.1, 0.1, n)
        self.roles   = ['captacion', 'captacion', 'captacion',
                        'bat_acelera', 'bat_frena'][:n]
        self.SoE     = np.random.uniform(0.4, 0.6, n)
        self.lambda_ = np.zeros(n)
        self.cp      = np.zeros(n)
        self.P_elec  = np.zeros(n)   # [MW]
        self.freq    = F0
        self.hist_f  = np.full(200, F0)
        self.t       = 0.0
        self.viento  = 9.0
        self.pert_activa = False
        self.pert_mag    = 200.0     # MW
        self.pert_restante = 0.0
        self.log     = ["[00:00] Sistema inicializado."]
        self.rotacion_timer = 0.0
        self._actualizar_metricas()

    def _actualizar_metricas(self):
        for i in range(self.n):
            v = max(self.viento, 1.0)
            self.lambda_[i] = self.omega[i] * R / v
            self.cp[i]      = Cp_curva(self.lambda_[i])
            T_gen           = T_MPPT(self.omega[i])
            if self.roles[i] == 'bat_acelera':
                T_gen *= 0.75
            elif self.roles[i] == 'bat_frena':
                T_gen *= 1.30
            self.P_elec[i]  = max(0.0, T_gen * self.omega[i] / 1e6)

    def step(self):
        self.t += DT
        # Viento turbulento
        v_inst = (self.viento
                  + 1.2 * np.sin(self.t * 0.3)
                  + 0.6 * np.sin(self.t * 1.1)
                  + 0.4 * (np.random.random() - 0.5))
        v_inst = max(v_inst, 1.0)

        for i in range(self.n):
            ta  = T_aerodinamico(v_inst, self.omega[i])
            tm  = T_MPPT(self.omega[i])
            buf = 0.0
            if self.roles[i] == 'bat_acelera':
                buf = -tm * 0.25
            elif self.roles[i] == 'bat_frena':
                buf =  tm * 0.30
            T_gen = max(0.0, tm + buf)
            T_roz = 0.02 * self.omega[i]
            d_omega = (ta - T_gen - T_roz) / J
            self.omega[i] = np.clip(
                self.omega[i] + d_omega * DT, OMEGA_MIN, OMEGA_MAX)

            lam = self.omega[i] * R / v_inst
            self.lambda_[i] = lam
            self.cp[i]      = Cp_curva(lam)
            self.P_elec[i]  = max(0.0, T_gen * self.omega[i] / 1e6)

            E_cin = 0.5 * J * self.omega[i]**2
            E_max = 0.5 * J * OMEGA_MAX**2
            E_min = 0.5 * J * OMEGA_MIN**2
            self.SoE[i] = np.clip((E_cin - E_min) / (E_max - E_min), 0.0, 1.0)

        # Swing equation: frecuencia de red
        P_total  = np.sum(self.P_elec) * 1e6
        P_demand = 8.0e6 + (self.pert_restante * 1e6 if self.pert_activa else 0.0)
        dP       = (P_total - P_demand) / S_TOTAL
        dF       = -dP * F0 / (2.0 * H_SIS)
        self.freq = np.clip(self.freq + dF * DT, 48.5, 51.5)
        self.hist_f = np.roll(self.hist_f, -1)
        self.hist_f[-1] = self.freq

        # Perturbación temporal
        if self.pert_activa:
            self.pert_restante -= DT / 5.0
            if self.pert_restante <= 0:
                self.pert_activa  = False
                self.pert_restante = 0.0
                self.log.append(f"[{self.t:.0f}s] ✓ Perturbación absorbida")
                self.log = self.log[-15:]

        # Rotación automática de roles cada 15 s
        self.rotacion_timer += DT
        if self.rotacion_timer >= 15.0:
            self.rotacion_timer = 0.0
            self._rotar_roles()

    def _rotar_roles(self):
        """Rota el rol bat_acelera a otro índice al azar."""
        idx_acel = [i for i, r in enumerate(self.roles) if r == 'bat_acelera']
        idx_cap  = [i for i, r in enumerate(self.roles) if r == 'captacion']
        if idx_acel and idx_cap:
            old = idx_acel[0]
            new = idx_cap[np.random.randint(len(idx_cap))]
            self.roles[old] = 'captacion'
            self.roles[new] = 'bat_acelera'
            self.log.append(f"[{self.t:.0f}s] ↻ Rol ACEL: T{old+1}→T{new+1}")
            self.log = self.log[-15:]

    def aplicar_reparto(self, pct_cap: float):
        n_cap  = max(1, round(pct_cap / 100.0 * self.n))
        n_acel = max(0, (self.n - n_cap) // 2)
        for i in range(self.n):
            if i < n_cap:
                self.roles[i] = 'captacion'
            elif i < n_cap + n_acel:
                self.roles[i] = 'bat_acelera'
            else:
                self.roles[i] = 'bat_frena'

    def inyectar_perturbacion(self, mag_mw: float):
        self.pert_activa   = True
        self.pert_restante = 1.0
        self.pert_mag      = mag_mw
        self.log.append(f"[{self.t:.0f}s] ⚡ Pert. {mag_mw:.0f} MW")
        self.log = self.log[-15:]

# ═══════════════════════════════════════════════════════════════════
# GENERADOR DE IMAGEN ESTÁTICA (modo sin display)
# ═══════════════════════════════════════════════════════════════════

def barra_h(ax, valor, maximo, color, label, unidad, y_pos, barra_h=0.6):
    """Dibuja una barra horizontal de porcentaje."""
    pct = np.clip(valor / maximo, 0, 1)
    danger = pct > 0.88 or pct < 0.12

    # Fondo de la barra
    ax.barh(y_pos, maximo, left=0, height=barra_h,
            color='#0a1520', edgecolor='#1a3a4a', linewidth=0.5)
    # Barra de valor
    bar_color = '#E74C3C' if danger else color
    ax.barh(y_pos, valor, left=0, height=barra_h,
            color=bar_color, alpha=0.85, edgecolor=bar_color, linewidth=0)

    # Etiqueta izquierda
    ax.text(-maximo * 0.02, y_pos, label,
            va='center', ha='right', fontsize=8, color='#7aabb0',
            fontfamily='monospace')
    # Valor derecha
    ax.text(maximo * 1.02, y_pos, f'{valor:.2f} {unidad}',
            va='center', ha='left', fontsize=8,
            color='#E74C3C' if danger else color,
            fontfamily='monospace', fontweight='bold')
    # Marcas en 25 50 75 %
    for tick in [0.25, 0.5, 0.75]:
        ax.axvline(maximo * tick, color='#2a3a4a', lw=0.5, alpha=0.7,
                   ymin=(y_pos - barra_h/2)/10, ymax=(y_pos + barra_h/2)/10)


def generar_imagen(sis: SistemaTurbinas, viento: float, reparto: float,
                   n_pasos: int, filename: str):
    """Ejecuta la simulación y genera la imagen final."""

    # ─── Simular n_pasos ───────────────────────────────────────────
    sis.viento = viento
    sis.aplicar_reparto(reparto)
    for _ in range(n_pasos):
        sis.step()

    # ─── Layout ────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 13), facecolor='#060a10')
    fig.patch.set_facecolor('#060a10')

    gs = GridSpec(4, 3, figure=fig,
                  left=0.06, right=0.97, top=0.93, bottom=0.05,
                  hspace=0.55, wspace=0.35)

    ax_header = fig.add_axes([0, 0.94, 1, 0.06], facecolor='#0a1520')
    ax_freq   = fig.add_subplot(gs[0, :2])
    ax_glob   = fig.add_subplot(gs[0, 2])
    axs_turb  = [fig.add_subplot(gs[1+i//3, i%3]) for i in range(N_TURB)]
    ax_log    = fig.add_subplot(gs[3, :2])
    ax_info   = fig.add_subplot(gs[3, 2])

    for ax in [ax_freq, ax_glob, ax_log, ax_info] + axs_turb:
        ax.set_facecolor('#090f1a')
        for spine in ax.spines.values():
            spine.set_edgecolor('#1a2a3a')
            spine.set_linewidth(0.7)

    # ─── CABECERA ──────────────────────────────────────────────────
    ax_header.set_xlim(0, 1); ax_header.set_ylim(0, 1)
    ax_header.axis('off')
    ax_header.text(0.5, 0.65, 'ZYPYZAPE  —  MINIGEMELO DIGITAL v1.0',
                   ha='center', va='center', fontsize=14, fontweight='bold',
                   color='#028090', fontfamily='monospace',
                   transform=ax_header.transAxes)
    ax_header.text(0.5, 0.15,
                   f'Víctor Manzanares Alberola · EPSA UPV Alcoy · 2026'
                   f'    |    t = {sis.t:.1f} s    |    N = {N_TURB} turbinas',
                   ha='center', va='center', fontsize=8,
                   color='#4a7a8a', fontfamily='monospace',
                   transform=ax_header.transAxes)

    # ─── FRECUENCIA DE RED ──────────────────────────────────────────
    f_color = ('#E74C3C' if abs(sis.freq - 50) > 0.5
               else '#F4A261' if abs(sis.freq - 50) > 0.2 else '#02C39A')
    ax_freq.set_facecolor('#080c14')
    ax_freq.plot(np.linspace(0, 10, len(sis.hist_f)), sis.hist_f,
                 color='#00ff88', lw=1.5, alpha=0.9)
    ax_freq.axhline(50.0, color='#2a5030', lw=0.8, ls='--', alpha=0.7)
    ax_freq.axhline(49.5, color='#3a2020', lw=0.5, ls=':', alpha=0.5)
    ax_freq.axhline(50.5, color='#3a2020', lw=0.5, ls=':', alpha=0.5)
    ax_freq.fill_between(np.linspace(0, 10, len(sis.hist_f)),
                         sis.hist_f, 50,
                         alpha=0.15, color='#00ff88')
    ax_freq.set_xlim(0, 10); ax_freq.set_ylim(49.0, 51.0)
    ax_freq.set_ylabel('f [Hz]', color='#5a8a9a', fontsize=8,
                       fontfamily='monospace')
    ax_freq.tick_params(colors='#3a5a6a', labelsize=7)
    ax_freq.set_title(
        f'FRECUENCIA DE RED    {sis.freq:.4f} Hz'
        f'   |   ΔF = {abs(sis.freq-50)*1000:.1f} mHz',
        color=f_color, fontsize=9, fontfamily='monospace', pad=4)

    # ─── MÉTRICAS GLOBALES ──────────────────────────────────────────
    P_total = np.sum(sis.P_elec)
    SoE_m   = np.mean(sis.SoE)
    cp_m    = np.mean(sis.cp)
    n_cap   = sum(1 for r in sis.roles if r == 'captacion')
    H_ef    = SoE_m * 0.5 * J * OMEGA_NOM**2 / P_NOM

    ax_glob.set_xlim(0, 2.5); ax_glob.set_ylim(-0.3, 8.5)
    ax_glob.axis('off')
    ax_glob.set_title('MÉTRICAS GLOBALES', color='#4a7a8a',
                      fontsize=8, fontfamily='monospace', pad=4)

    metricas = [
        ('P TOTAL',     f'{P_total:.2f} MW',  '#7EC8E3'),
        ('FREQ',        f'{sis.freq:.4f} Hz',  f_color),
        ('SoE MEDIO',   f'{SoE_m*100:.1f} %',  '#F4A261'),
        ('Cp MEDIO',    f'{cp_m:.4f}',          '#02C39A'),
        ('H EFECTIVO',  f'{H_ef:.2f} s',        '#9B59B6'),
        ('% CAPTACIÓN', f'{n_cap/N_TURB*100:.0f} %', '#028090'),
        ('VIENTO',      f'{viento:.1f} m/s',    '#5ab8c8'),
        ('t SIM',       f'{sis.t:.1f} s',        '#4a7a8a'),
    ]
    for i, (lbl, val, col) in enumerate(metricas):
        y = 7.8 - i * 1.0
        ax_glob.add_patch(FancyBboxPatch((0, y-0.38), 2.5, 0.76,
                          boxstyle="round,pad=0.05",
                          facecolor='#0a1520', edgecolor='#1a3040', lw=0.5))
        ax_glob.text(0.1, y, lbl, va='center', fontsize=7.5,
                     color='#4a7a8a', fontfamily='monospace')
        ax_glob.text(2.4, y, val, va='center', ha='right', fontsize=9,
                     color=col, fontfamily='monospace', fontweight='bold')

    # ─── BARRAS POR TURBINA ─────────────────────────────────────────
    for i, ax in enumerate(axs_turb):
        if i >= N_TURB:
            ax.axis('off')
            continue
        rol_key = sis.roles[i]
        rol     = ROLES[rol_key]
        omega_pct = (sis.omega[i] - OMEGA_MIN) / (OMEGA_MAX - OMEGA_MIN)

        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 4.5)
        ax.axis('off')
        ax.set_title(f'{sis.nombres[i]}   [{rol["label"]}]',
                     color=rol['color'], fontsize=9,
                     fontfamily='monospace', pad=3, fontweight='bold')

        barras = [
            (sis.SoE[i],    1.0,  '#F4A261', 'SoE', '%×100'),
            (omega_pct,     1.0,  '#028090', 'ω',   '%rng'),
            (sis.cp[i],     0.5,  '#02C39A', 'Cp',  ''),
            (sis.P_elec[i], 2.5,  '#7EC8E3', 'P',   'MW'),
        ]
        for j, (val, mx, col, lbl, unit) in enumerate(barras):
            y = 3.6 - j * 1.0
            pct = np.clip(val / mx, 0, 1)
            danger = pct > 0.88 or pct < 0.12
            c = '#E74C3C' if danger else col

            # Fondo
            ax.barh(y, 1.0, left=0, height=0.55,
                    color='#0a1520', edgecolor='#1a3040', lw=0.5)
            # Relleno
            ax.barh(y, pct, left=0, height=0.55,
                    color=c, alpha=0.85)
            # Brillo en borde derecho del relleno
            if pct > 0.02:
                ax.barh(y, min(pct, 0.02), left=max(0, pct-0.02),
                        height=0.55, color='white', alpha=0.3)

            # Etiquetas
            ax.text(-0.02, y, lbl, va='center', ha='right',
                    fontsize=7, color='#5a8a9a', fontfamily='monospace')
            # Valor numérico dentro de la barra
            if lbl == 'SoE':
                vstr = f'{val*100:.0f}%'
            elif lbl == 'ω':
                vstr = f'{sis.omega[i]:.2f}'
            elif lbl == 'Cp':
                vstr = f'{val:.3f}'
            else:
                vstr = f'{val:.2f}'
            ax.text(min(pct, 0.98), y, f' {vstr}', va='center',
                    ha='left' if pct < 0.6 else 'right',
                    fontsize=7, color='white' if pct > 0.15 else c,
                    fontfamily='monospace')

            # Ticks en 25/50/75%
            for tick in [0.25, 0.5, 0.75]:
                ax.axvline(tick, ymin=(y-0.27)/5, ymax=(y+0.27)/5,
                           color='#2a3a4a', lw=0.5, alpha=0.6)

        # Indicador de λ bajo la tarjeta
        ax.text(0.5, -0.35, f'λ = {sis.lambda_[i]:.2f}',
                ha='center', va='center', fontsize=7,
                color='#3a6a7a', fontfamily='monospace')

    # ─── BARRA GLOBAL FRECUENCIA ───────────────────────────────────
    ax_glob2 = fig.add_axes([0.06, 0.03, 0.63, 0.015])
    ax_glob2.set_facecolor('#080c14')
    ax_glob2.set_xlim(48.5, 51.5); ax_glob2.set_ylim(0, 1)
    ax_glob2.axis('off')

    pct_f = (sis.freq - 48.5) / 3.0
    ax_glob2.barh(0.5, 3.0, left=0, height=1.0,
                  color='#0a1520', edgecolor='#1a3040', lw=0.5)
    ax_glob2.barh(0.5, sis.freq - 48.5, left=0, height=1.0,
                  color=f_color, alpha=0.8)
    ax_glob2.axvline(50.0 - 48.5, color='#ffffff', lw=1.5, alpha=0.4)
    ax_glob2.text(0, 0.5, '48.5 Hz', va='center', ha='right',
                  fontsize=6, color='#3a5a6a', fontfamily='monospace')
    ax_glob2.text(3.0, 0.5, '51.5 Hz', va='center', ha='left',
                  fontsize=6, color='#3a5a6a', fontfamily='monospace')
    ax_glob2.text(1.5, 0.5, f'f = {sis.freq:.4f} Hz',
                  va='center', ha='center', fontsize=7,
                  color=f_color, fontfamily='monospace', fontweight='bold')

    # ─── LOG DE EVENTOS ─────────────────────────────────────────────
    ax_log.set_xlim(0, 1); ax_log.set_ylim(0, 1)
    ax_log.axis('off')
    ax_log.set_title('LOG DE EVENTOS', color='#4a7a8a',
                     fontsize=8, fontfamily='monospace', pad=4)
    for i, line in enumerate(sis.log[-8:]):
        ax_log.text(0.02, 0.88 - i * 0.12, line, va='top',
                    fontsize=7, color='#5a9a7a', fontfamily='monospace')

    # ─── INFO / LEYENDA ─────────────────────────────────────────────
    ax_info.set_xlim(0, 1); ax_info.set_ylim(0, 1)
    ax_info.axis('off')
    ax_info.set_title('LEYENDA & FÍSICA', color='#4a7a8a',
                      fontsize=8, fontfamily='monospace', pad=4)
    info = [
        ('CAPT', '#02C39A', 'Captación MPPT: T_gen = K_opt·ω²'),
        ('ACEL', '#F4A261', 'Batería acelera: T_gen = 0.75·T_opt'),
        ('FREN', '#E74C3C', 'Batería frena:  T_gen = 1.30·T_opt'),
    ]
    for i, (lbl, col, desc) in enumerate(info):
        y = 0.82 - i * 0.16
        ax_info.add_patch(FancyBboxPatch((0.01, y-0.05), 0.07, 0.10,
                          boxstyle="round,pad=0.01",
                          facecolor=col, alpha=0.8))
        ax_info.text(0.10, y, lbl, va='center', fontsize=7.5,
                     color=col, fontfamily='monospace', fontweight='bold')
        ax_info.text(0.20, y, desc, va='center', fontsize=6.5,
                     color='#5a8a9a', fontfamily='monospace')
    ax_info.axhline(0.30, color='#1a3040', lw=0.5)
    eqs = [
        'J·dω/dt = T_aero − T_gen − T_roz',
        'λ = ω·R/v        Cp_max = 0.593 (Betz)',
        'RoCoF = −ΔP·f₀ / (2H·S_tot)',
        'SoE = (E_cin−E_min)/(E_max−E_min)',
    ]
    for i, eq in enumerate(eqs):
        ax_info.text(0.02, 0.26 - i * 0.07, eq, va='center',
                     fontsize=6.5, color='#4a7a8a', fontfamily='monospace')

    plt.savefig(filename, dpi=150, bbox_inches='tight',
                facecolor='#060a10', edgecolor='none')
    plt.close(fig)
    print(f"✓ Imagen guardada: {filename}")


# ═══════════════════════════════════════════════════════════════════
# MODO INTERACTIVO CON MATPLOTLIB (animación)
# ═══════════════════════════════════════════════════════════════════

def modo_interactivo():
    """Versión animada con sliders y botones (requiere display gráfico)."""
    import matplotlib
    matplotlib.use('TkAgg')  # o 'Qt5Agg'
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.widgets import Slider, Button

    sis   = SistemaTurbinas()
    pause = [False]

    fig   = plt.figure(figsize=(17, 11), facecolor='#060a10')
    fig.patch.set_facecolor('#060a10')
    plt.rcParams['text.color'] = '#e0e8f0'

    gs = GridSpec(5, 3, figure=fig,
                  left=0.07, right=0.97, top=0.93, bottom=0.14,
                  hspace=0.6, wspace=0.35)

    ax_freq  = fig.add_subplot(gs[0, :2])
    ax_glob  = fig.add_subplot(gs[0, 2])
    axs_turb = [fig.add_subplot(gs[1+i//3, i%3]) for i in range(N_TURB)]
    ax_log   = fig.add_subplot(gs[3, :2])
    ax_info  = fig.add_subplot(gs[3, 2])

    for ax in [ax_freq, ax_glob, ax_log, ax_info] + axs_turb:
        ax.set_facecolor('#090f1a')
        for sp in ax.spines.values():
            sp.set_edgecolor('#1a2a3a'); sp.set_linewidth(0.7)

    # Título
    fig.text(0.5, 0.96, 'ZYPYZAPE — MINIGEMELO DIGITAL v1.0',
             ha='center', fontsize=13, fontweight='bold',
             color='#028090', fontfamily='monospace')
    fig.text(0.5, 0.935,
             'Simulación dinámica · 5 aerogeneradores · Batería cinética sintética',
             ha='center', fontsize=8, color='#4a7a8a', fontfamily='monospace')

    # ── Sliders ──────────────────────────────────────────────────────
    ax_s_viento  = fig.add_axes([0.08, 0.085, 0.25, 0.025], facecolor='#0a1520')
    ax_s_reparto = fig.add_axes([0.40, 0.085, 0.25, 0.025], facecolor='#0a1520')
    sl_viento    = Slider(ax_s_viento,  'VIENTO m/s', 4, 18, valinit=9,
                          color='#028090', initcolor='none')
    sl_reparto   = Slider(ax_s_reparto, 'CAPTACIÓN %', 20, 100, valinit=60,
                          color='#02C39A', initcolor='none')
    for sl in [sl_viento, sl_reparto]:
        sl.label.set_color('#5a9aaa')
        sl.label.set_fontfamily('monospace')
        sl.valtext.set_color('#7ae8f0')
        sl.valtext.set_fontfamily('monospace')

    # ── Botones ───────────────────────────────────────────────────────
    ax_b_pert  = fig.add_axes([0.73, 0.075, 0.11, 0.04], facecolor='#1a0808')
    ax_b_pause = fig.add_axes([0.85, 0.075, 0.07, 0.04], facecolor='#0a1020')
    ax_b_reset = fig.add_axes([0.93, 0.075, 0.05, 0.04], facecolor='#0a1020')
    b_pert     = Button(ax_b_pert,  '⚡ PERT.', color='#1a0808', hovercolor='#3a1010')
    b_pause    = Button(ax_b_pause, '⏸ PAUSE', color='#0a1020', hovercolor='#1a2030')
    b_reset    = Button(ax_b_reset, '↺ RST',   color='#0a1020', hovercolor='#1a2030')
    for btn in [b_pert, b_pause, b_reset]:
        btn.label.set_color('#e0e8f0')
        btn.label.set_fontfamily('monospace')
        btn.label.set_fontsize(8)

    def on_pert(_):
        sis.inyectar_perturbacion(sl_viento.val * 20)

    def on_pause(_):
        pause[0] = not pause[0]

    def on_reset(_):
        nonlocal sis
        sis = SistemaTurbinas()

    b_pert.on_clicked(on_pert)
    b_pause.on_clicked(on_pause)
    b_reset.on_clicked(on_reset)

    def update(_):
        if pause[0]:
            return
        sis.viento = sl_viento.val
        sis.aplicar_reparto(sl_reparto.val)
        for _ in range(3):
            sis.step()

        # Limpiar ejes dinámicos
        for ax in [ax_freq, ax_glob, ax_log] + axs_turb:
            ax.cla()
            ax.set_facecolor('#090f1a')
            for sp in ax.spines.values():
                sp.set_edgecolor('#1a2a3a'); sp.set_linewidth(0.7)

        # ── Frecuencia ────────────────────────────────────────────────
        f_color = ('#E74C3C' if abs(sis.freq-50)>0.5
                   else '#F4A261' if abs(sis.freq-50)>0.2 else '#02C39A')
        ax_freq.set_facecolor('#080c14')
        ax_freq.plot(np.linspace(0,10,len(sis.hist_f)), sis.hist_f,
                     color='#00ff88', lw=1.2)
        ax_freq.axhline(50, color='#2a5030', lw=0.8, ls='--', alpha=0.7)
        ax_freq.fill_between(np.linspace(0,10,len(sis.hist_f)),
                             sis.hist_f, 50, alpha=0.12, color='#00ff88')
        ax_freq.set_ylim(49, 51); ax_freq.set_xlim(0, 10)
        ax_freq.set_ylabel('Hz', color='#5a8a9a', fontsize=8,
                            fontfamily='monospace')
        ax_freq.tick_params(colors='#3a5a6a', labelsize=7)
        ax_freq.set_title(
            f'FRECUENCIA RED   {sis.freq:.4f} Hz   |'
            f'   ΔF = {abs(sis.freq-50)*1000:.1f} mHz   |'
            f'   t = {sis.t:.0f} s',
            color=f_color, fontsize=8, fontfamily='monospace', pad=3)

        # ── Métricas globales ─────────────────────────────────────────
        P_total = np.sum(sis.P_elec)
        SoE_m   = np.mean(sis.SoE)
        cp_m    = np.mean(sis.cp)
        n_cap   = sum(1 for r in sis.roles if r=='captacion')
        H_ef    = SoE_m * 0.5 * J * OMEGA_NOM**2 / P_NOM
        ax_glob.set_xlim(0, 2.5); ax_glob.set_ylim(-0.5, 8.5)
        ax_glob.axis('off')
        ax_glob.set_title('MÉTRICAS GLOBALES', color='#4a7a8a',
                          fontsize=8, fontfamily='monospace', pad=3)
        mets = [
            ('P TOTAL',    f'{P_total:.2f} MW',  '#7EC8E3'),
            ('FREQ',       f'{sis.freq:.4f} Hz',  f_color),
            ('SoE MEDIO',  f'{SoE_m*100:.1f} %',  '#F4A261'),
            ('Cp MEDIO',   f'{cp_m:.4f}',          '#02C39A'),
            ('H EFECTIVO', f'{H_ef:.2f} s',        '#9B59B6'),
            ('%CAPTACIÓN', f'{n_cap/N_TURB*100:.0f} %', '#028090'),
            ('VIENTO',     f'{sis.viento:.1f} m/s','#5ab8c8'),
        ]
        for idx, (lbl, val, col) in enumerate(mets):
            y = 7.8 - idx * 1.1
            ax_glob.add_patch(FancyBboxPatch((0,y-0.4),2.5,0.8,
                              boxstyle="round,pad=0.04",
                              facecolor='#0a1520',edgecolor='#1a3040',lw=0.4))
            ax_glob.text(0.08, y, lbl, va='center', fontsize=7,
                         color='#4a7a8a', fontfamily='monospace')
            ax_glob.text(2.42, y, val, va='center', ha='right',
                         fontsize=8.5, color=col,
                         fontfamily='monospace', fontweight='bold')

        # ── Tarjetas por turbina ──────────────────────────────────────
        for i, ax in enumerate(axs_turb):
            if i >= N_TURB: ax.axis('off'); continue
            rol     = ROLES[sis.roles[i]]
            om_pct  = (sis.omega[i]-OMEGA_MIN)/(OMEGA_MAX-OMEGA_MIN)
            ax.set_xlim(0, 1); ax.set_ylim(-0.5, 4.5); ax.axis('off')
            ax.set_title(f'{sis.nombres[i]}  [{rol["label"]}]',
                         color=rol['color'], fontsize=8.5,
                         fontfamily='monospace', pad=2)
            barras = [
                (sis.SoE[i], 1.0, '#F4A261', 'SoE',
                 f'{sis.SoE[i]*100:.0f}%'),
                (om_pct,     1.0, '#028090', 'ω',
                 f'{sis.omega[i]:.2f}'),
                (sis.cp[i],  0.5, '#02C39A', 'Cp',
                 f'{sis.cp[i]:.3f}'),
                (sis.P_elec[i], 2.5, '#7EC8E3', 'P',
                 f'{sis.P_elec[i]:.2f}MW'),
            ]
            for j, (val, mx, col, lbl, vstr) in enumerate(barras):
                y    = 3.7 - j * 1.0
                pct  = np.clip(val/mx, 0, 1)
                c    = '#E74C3C' if (pct>0.88 or pct<0.12) else col
                ax.barh(y, 1.0, left=0, height=0.60,
                        color='#0a1520', edgecolor='#1a3040', lw=0.4)
                ax.barh(y, pct, left=0, height=0.60, color=c, alpha=0.85)
                ax.text(-0.02, y, lbl, va='center', ha='right',
                        fontsize=6.5, color='#5a8a9a', fontfamily='monospace')
                ax.text(min(pct+0.02, 0.97), y, vstr, va='center',
                        ha='left' if pct<0.55 else 'right',
                        fontsize=6.5, color='white' if pct>0.15 else c,
                        fontfamily='monospace')
                for tick in [0.25, 0.5, 0.75]:
                    ax.axvline(tick, ymin=(y-0.3)/5, ymax=(y+0.3)/5,
                               color='#2a3a4a', lw=0.4, alpha=0.6)
            ax.text(0.5, -0.35, f'λ={sis.lambda_[i]:.2f}',
                    ha='center', fontsize=6.5, color='#3a6a7a',
                    fontfamily='monospace')

        # ── Log ───────────────────────────────────────────────────────
        ax_log.set_xlim(0,1); ax_log.set_ylim(0,1); ax_log.axis('off')
        ax_log.set_title('LOG', color='#4a7a8a',
                         fontsize=8, fontfamily='monospace', pad=2)
        for idx, line in enumerate(sis.log[-7:]):
            ax_log.text(0.01, 0.90 - idx*0.13, line, va='top',
                        fontsize=7, color='#5a9a7a', fontfamily='monospace')

    ani = FuncAnimation(fig, update, interval=120, cache_frame_data=False)
    plt.show()
    return ani


# ═══════════════════════════════════════════════════════════════════
# FUNCIÓN MULTI-ESCENARIO — genera paneles comparativos
# ═══════════════════════════════════════════════════════════════════

def generar_comparativa(filename: str):
    """Genera imagen comparativa de 4 escenarios de reparto."""
    escenarios = [
        (9.0, 90, '90% captación / 10% batería\n(modo producción)'),
        (9.0, 60, '60% captación / 40% batería\n(modo equilibrado)'),
        (12.0, 40, '40% captación / 60% batería\n(viento fuerte + carga batería)'),
        (7.0, 50, '50% captación / 50% batería\n(modo regulación)'),
    ]

    fig, axes = plt.subplots(4, 1, figsize=(16, 14), facecolor='#060a10')
    fig.patch.set_facecolor('#060a10')
    fig.suptitle('ZYPYZAPE — COMPARATIVA DE ESCENARIOS DE REPARTO',
                 color='#028090', fontsize=12, fontfamily='monospace',
                 fontweight='bold', y=0.98)

    for ax, (viento, rep, titulo) in zip(axes, escenarios):
        sis = SistemaTurbinas()
        sis.viento = viento
        sis.aplicar_reparto(rep)
        sis.inyectar_perturbacion(300)
        for _ in range(200):
            sis.step()

        ax.set_facecolor('#080c14')
        for sp in ax.spines.values():
            sp.set_edgecolor('#1a2a3a'); sp.set_linewidth(0.6)

        # Frecuencia
        ax.plot(np.linspace(0, 10, len(sis.hist_f)), sis.hist_f,
                color='#00ff88', lw=1.2, label='f(t)')
        ax.axhline(50, color='#2a5030', lw=0.8, ls='--', alpha=0.7,
                   label='f₀=50 Hz')
        ax.fill_between(np.linspace(0, 10, len(sis.hist_f)),
                        sis.hist_f, 50, alpha=0.15, color='#00ff88')
        ax.set_ylim(49.2, 50.8); ax.set_xlim(0, 10)
        ax.tick_params(colors='#3a5a6a', labelsize=7)
        ax.set_ylabel('Hz', color='#5a8a9a', fontsize=8,
                      fontfamily='monospace')

        # KPIs como texto
        P_t  = np.sum(sis.P_elec)
        SoE  = np.mean(sis.SoE)
        H_ef = SoE * 0.5 * J * OMEGA_NOM**2 / P_NOM
        f_col = ('#E74C3C' if abs(sis.freq-50)>0.3 else '#F4A261'
                 if abs(sis.freq-50)>0.1 else '#02C39A')
        info = (f'P={P_t:.2f} MW  |  SoE={SoE*100:.0f}%  |'
                f'  f={sis.freq:.3f} Hz  |  H_ef={H_ef:.2f} s')
        ax.set_title(f'{titulo}    →    {info}',
                     color=f_col, fontsize=8, fontfamily='monospace', pad=3)

        # Barras SoE por turbina como anotaciones
        for i in range(N_TURB):
            x = 0.3 + i * 1.8
            ax.bar(x, sis.SoE[i] * 0.3, bottom=49.2, width=0.8,
                   color=ROLES[sis.roles[i]]['color'], alpha=0.6)
            ax.text(x, 49.18, sis.nombres[i], ha='center', fontsize=5.5,
                    color='#5a8a9a', fontfamily='monospace')

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(filename, dpi=130, bbox_inches='tight',
                facecolor='#060a10', edgecolor='none')
    plt.close(fig)
    print(f"✓ Comparativa guardada: {filename}")


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import sys
    import os

    print("╔══════════════════════════════════════════════════════╗")
    print("║   ZYPYZAPE MiniGemelo Digital v1.0                  ║")
    print("║   Víctor Manzanares Alberola · UPV-EPSA · 2026      ║")
    print("╚══════════════════════════════════════════════════════╝")
    print()

    # Intentar modo interactivo primero
    try:
        import matplotlib
        pass
        import tkinter  # sólo para comprobar que hay display
        print("► Modo interactivo detectado — abriendo ventana animada...")
        ani = modo_interactivo()

    except Exception as e:
        print(f"► Sin display gráfico ({e})")
        print("► Modo estático: generando imágenes PNG...")
        print()

        out_dir = os.path.join(os.path.dirname(__file__), 'output_gemelo')
        os.makedirs(out_dir, exist_ok=True)

        # Imagen 1: Estado inicial (t=10s)
        sis = SistemaTurbinas()
        generar_imagen(sis, viento=9.0, reparto=80,
                       n_pasos=200,
                       filename=os.path.join(out_dir, '01_estado_inicial.png'))

        # Imagen 2: Tras perturbación
        sis2 = SistemaTurbinas()
        sis2.inyectar_perturbacion(300)
        generar_imagen(sis2, viento=9.0, reparto=60,
                       n_pasos=300,
                       filename=os.path.join(out_dir, '02_post_perturbacion.png'))

        # Imagen 3: Viento fuerte + máxima batería
        sis3 = SistemaTurbinas()
        generar_imagen(sis3, viento=14.0, reparto=40,
                       n_pasos=300,
                       filename=os.path.join(out_dir, '03_viento_fuerte_bat.png'))

        # Comparativa de escenarios
        generar_comparativa(
            os.path.join(out_dir, '04_comparativa_escenarios.png'))

        print()
        print(f"✓ 4 imágenes generadas en: {out_dir}/")
        print()
        print("Para modo interactivo (requiere display):")
        print("  pip install matplotlib numpy")
        print("  python3 zypyzape_minigemelo.py")
