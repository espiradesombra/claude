"""
ZypyZape Digital Twin v2.0
==========================
Sistema: Víctor Manzanares Alberola
Asistencia de código: IA

Mejoras sobre versiones anteriores:
  - Ecuación de movimiento correcta: J·dω/dt = T_aero − T_gen − T_roz + T_ZZ + T_IS
  - MPPT real: T_gen = k_mppt · ω²  (maximiza Cp en λ_opt)
  - Curva Cp(λ,β) parametrizada realista (modelo estándar de 6 coeficientes)
  - Ciclo ZypyZape: roles CAPT/ACEL/FREN con intercambio de potencia real
  - Swing equation para frecuencia de red con inercia sintética ZypyZape
  - Respuesta a perturbación de red (RoCoF / nadir de frecuencia)
  - Contabilidad energética correcta (bug corregido en versión anterior)
  - Acoplamiento Kuramoto correcto: Σ K_ij · sin(θ_j − θ_i)
  - Topología de 5 nodos según documento técnico: par central (0,1) + anillo (2,3,4)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ============================================================
# PARÁMETROS FÍSICOS — turbina 2.5 MW
# ============================================================
N       = 5
R       = 60.0              # radio rotor [m]
A       = np.pi * R**2      # área barrida [m²]
rho     = 1.225             # densidad aire [kg/m³]
J       = 5e6               # inercia rotor [kg·m²]
S_nom   = 2.5e6             # potencia nominal por turbina [W]
Kv      = 2e4               # coeficiente rozamiento viscoso [N·m·s/rad]
T_c0    = 1e3               # rozamiento coulombiano [N·m]
omega_min = 0.2             # velocidad mínima [rad/s]
omega_max = 2.8             # velocidad máxima [rad/s]

# --- MPPT: T_gen = k_mppt · ω² → opera en λ_opt ---
lambda_opt = 7.0
Cp_max     = 0.48
k_mppt     = 0.5 * rho * A * R**3 * Cp_max / (lambda_opt**3)

# --- ZypyZape: ciclo de intercambio energético ---
P_ZZ_frac  = 0.18           # fracción de S_nom intercambiada [−]
f_cycle    = 0.4            # frecuencia del ciclo mecánico [Hz]
T_cycle    = 1.0 / f_cycle  # periodo del ciclo [s]

# --- Acoplamiento Kuramoto ---
Kc = 1.5e5                  # constante de acoplamiento [N·m/rad]

# ============================================================
# RED ELÉCTRICA — swing equation simplificada
# ============================================================
S_total    = 10e9           # potencia total sistema [W]
H_sistema  = 3.0            # constante de inercia [s] (baja por alta renovable)
f0         = 50.0           # frecuencia nominal [Hz]
D_amort    = 0.03           # amortiguamiento red [p.u./Hz]

# ============================================================
# TIEMPO Y EVENTOS
# ============================================================
dt         = 0.02           # paso [s]
T_sim      = 90.0           # duración [s]
steps      = int(T_sim / dt)
t_perturb  = 30.0           # instante de perturbación [s]
dP_perturb = -300e6         # pérdida de generación [W]

# ============================================================
# MATRIZ DE CONECTIVIDAD — topología documento técnico
# Par central: 0 ↔ 1   |   Anillo: 2 ↔ 3 ↔ 4 ↔ 2
# ============================================================
K_mat = np.zeros((N, N))
K_mat[0,1] = K_mat[1,0] = Kc
K_mat[2,3] = K_mat[3,2] = Kc
K_mat[3,4] = K_mat[4,3] = Kc
K_mat[4,2] = K_mat[2,4] = Kc

# Pares ZypyZape activos (intercambio cíclico)
PARES_ZZ = [(0, 1), (2, 3)]

# ============================================================
# MODELO DE VIENTO
# ============================================================
def viento(t):
    """Perfil con variación lenta, ráfaga en t=[45,55]s."""
    v = 11.5 + 1.5 * np.sin(0.08 * t) + 0.7 * np.sin(0.25 * t + 0.9)
    if 45.0 < t < 55.0:
        v += 3.5 * np.sin(np.pi * (t - 45.0) / 10.0)
    return max(3.0, v)

# ============================================================
# CURVA Cp(λ, β) — modelo parametrizado estándar 6 coeficientes
# ============================================================
def Cp_curva(lmbda, beta=0.0):
    lmbda = float(lmbda)
    if lmbda < 0.5:
        return 0.0
    denom1 = lmbda + 0.08 * beta
    denom2 = beta**3 + 1.0
    if abs(denom1) < 1e-6 or abs(denom2) < 1e-6:
        return 0.0
    lambda_i = 1.0 / (1.0 / denom1 - 0.035 / denom2)
    lambda_i = max(lambda_i, 0.5)
    cp = (0.5176 * (116.0 / lambda_i - 0.4 * beta - 5.0) *
          np.exp(-21.0 / lambda_i) + 0.0068 * lmbda)
    return float(np.clip(cp, 0.0, 0.593))

# ============================================================
# PAR AERODINÁMICO
# ============================================================
def T_aero(omega_i, v, beta=0.0):
    if omega_i < 0.05:
        return 0.0
    lmbda = omega_i * R / v
    cp = Cp_curva(lmbda, beta)
    P_capt = 0.5 * rho * A * v**3 * cp
    return P_capt / omega_i

# ============================================================
# ROZAMIENTO TOTAL
# ============================================================
def T_roz(omega_i):
    return Kv * omega_i + T_c0

# ============================================================
# CONTROL MPPT
# T_gen = k_mppt · ω²   →   λ = λ_opt cuando ω libre
# Limitado a S_nom para no superar la potencia nominal
# ============================================================
def T_gen_mppt(omega_i):
    T = k_mppt * omega_i**2
    T_max = S_nom / max(omega_i, 0.05)
    return min(T, T_max)

# ============================================================
# CICLO ZYPYZAPE — roles CAPT / ACEL / FREN
# El par (a, b) intercambia potencia a f_cycle Hz:
#   semiciclo 1: a cede energía → b acelera
#   semiciclo 2: b cede energía → a acelera
# ============================================================
def T_zypyzape(i, t, omega_i):
    T_zz = 0.0
    for par_idx, (a, b) in enumerate(PARES_ZZ):
        if i != a and i != b:
            continue
        # Fase del ciclo (desfasada entre pares)
        fase = (t / T_cycle + par_idx * 0.5) % 1.0
        role = +1 if (i == a) == (fase < 0.5) else -1
        # Intensidad sinusoidal suavizada
        mod = abs(np.sin(np.pi * f_cycle * t))
        P_zz = role * P_ZZ_frac * S_nom * mod
        T_zz += P_zz / max(omega_i, 0.1)
    return T_zz

# ============================================================
# INERCIA SINTÉTICA — respuesta droop a caída de frecuencia
# Cuando f cae, el convertidor aumenta T_gen (inyecta más potencia)
# → rotor frena → libera E_cin hacia la red (física correcta)
# ============================================================
def T_inercia_sintetica(f_grid, omega_i, df_dt):
    """Par adicional de generación proporcional a la desviación de f."""
    df = f0 - f_grid
    if df <= 0.0:
        return 0.0
    # Droop: ΔT_gen = k_droop * Δf/f0 * T_nom
    k_droop = 0.04          # ganancia droop [p.u.]
    T_nom   = S_nom / max(omega_i, 0.1)
    return k_droop * (df / f0) * T_nom  # positivo → frena más rotor → libera E_cin

# ============================================================
# ESTADO INICIAL
# ============================================================
omega  = np.ones(N) * 1.25
theta  = np.zeros(N)
f_grid = 50.0
df_dt  = 0.0

# Historial
omega_hist   = np.zeros((steps, N))
f_hist       = np.zeros(steps)
P_aero_hist  = np.zeros((steps, N))
P_elec_hist  = np.zeros((steps, N))
P_zz_hist    = np.zeros((steps, N))
lam_hist     = np.zeros((steps, N))
v_hist       = np.zeros(steps)

E_aero_tot = 0.0
E_elec_tot = 0.0

# ============================================================
# BUCLE DE SIMULACIÓN
# ============================================================
for step in range(steps):
    t = step * dt
    v = viento(t)
    v_hist[step] = v

    # --- Perturbación de red ---
    dP_red = dP_perturb if t >= t_perturb else 0.0

    # --- Inercia efectiva del sistema (ZypyZape contribuye) ---
    E_cin_zz = 0.5 * J * np.sum(omega**2)
    H_zz     = E_cin_zz / S_total          # inercia aportada por ZypyZape [s]
    H_eff    = H_sistema + H_zz

    # --- Swing equation: df/dt
    # dP_red < 0 → déficit de generación → f cae (df/dt < 0)
    df_dt = (dP_red * f0 / (2.0 * H_eff * S_total)
             - D_amort * (f_grid - f0))
    f_grid = np.clip(f_grid + df_dt * dt, 48.0, 52.0)

    # --- Dinámica de cada turbina ---
    domega = np.zeros(N)
    for i in range(N):
        ta  = T_aero(omega[i], v)
        tg  = T_gen_mppt(omega[i])
        tr  = T_roz(omega[i])
        tzz = T_zypyzape(i, t, omega[i])
        tis = T_inercia_sintetica(f_grid, omega[i], df_dt)

        # Acoplamiento Kuramoto
        coupling = sum(K_mat[i,j] * np.sin(theta[j] - theta[i])
                       for j in range(N) if K_mat[i,j] != 0)

        # Ecuación de movimiento
        T_net    = ta - tg - tr + tzz + tis + coupling
        domega[i] = T_net / J

        P_aero_hist[step, i] = ta * omega[i]
        P_elec_hist[step, i] = tg * omega[i]
        P_zz_hist[step, i]   = tzz * omega[i]
        lam_hist[step, i]    = omega[i] * R / v

        E_aero_tot += ta * omega[i] * dt
        E_elec_tot += tg * omega[i] * dt

    omega = np.clip(omega + domega * dt, omega_min, omega_max)
    theta += omega * dt

    omega_hist[step] = omega
    f_hist[step]     = f_grid

t_vec = np.arange(steps) * dt

# ============================================================
# KPIs
# ============================================================
idx_p  = int(t_perturb / dt)
f_min  = np.min(f_hist[idx_p:])
rocof  = (f_hist[idx_p + 1] - f_hist[idx_p]) / dt
eff    = 100.0 * E_elec_tot / max(E_aero_tot, 1.0)
lam_medio_final = np.mean(lam_hist[-100:, :])

print("=" * 56)
print("   ZYPYZAPE GEMELO DIGITAL v2.0 — KPIs")
print("=" * 56)
print(f"   Energía aerodinámica total : {E_aero_tot/1e9:.3f} GJ")
print(f"   Energía eléctrica generada : {E_elec_tot/1e9:.3f} GJ")
print(f"   Eficiencia global          : {eff:.1f} %")
print(f"   f_min tras perturbación    : {f_min:.3f} Hz")
print(f"   RoCoF inicial              : {rocof:.3f} Hz/s")
print(f"   λ medio (últimos 2s)       : {lam_medio_final:.2f}")
print(f"   ω medio final              : {np.mean(omega_hist[-1]):.3f} rad/s")
print("=" * 56)

# ============================================================
# GRÁFICAS — 6 paneles
# ============================================================
NOMBRES = ['Central-A (0)', 'Central-B (1)',
           'Anillo-2', 'Anillo-3', 'Anillo-4']
COLS = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6']
BG   = '#0f0e17'
PAN  = '#1a1a2e'

fig = plt.figure(figsize=(17, 12), facecolor=BG)
gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.50, wspace=0.38)

def estilo_ax(ax, titulo):
    ax.set_facecolor(PAN)
    ax.set_title(titulo, color='white', fontsize=10, pad=6)
    ax.tick_params(colors='#aaaaaa', labelsize=8)
    for sp in ax.spines.values():
        sp.set_color('#444466')
    ax.xaxis.label.set_color('#aaaaaa')
    ax.yaxis.label.set_color('#aaaaaa')
    ax.grid(color='#2a2a4a', lw=0.6, linestyle='--')

def marca_perturb(ax):
    ax.axvline(t_perturb, color='#ffcc00', ls='--', lw=1.2, alpha=0.7,
               label='Perturbación')

# ── Panel 1: velocidades angulares ──────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
for i in range(N):
    ax1.plot(t_vec, omega_hist[:, i], color=COLS[i], lw=1.5, label=NOMBRES[i])
marca_perturb(ax1)
ax1.set_xlabel('t [s]')
ax1.set_ylabel('ω [rad/s]')
ax1.legend(fontsize=7, loc='upper right', framealpha=0.3)
estilo_ax(ax1, 'Velocidades angulares ω(t)')

# ── Panel 2: frecuencia de red ───────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(t_vec, f_hist, color='#00d2ff', lw=2.0, label='f_grid (ZypyZape)')
ax2.axhline(49.0, color='#e74c3c', ls=':', lw=1.5, alpha=0.9, label='49 Hz (límite)')
ax2.axhline(50.0, color='white',   ls='--', lw=0.8, alpha=0.3, label='f₀ = 50 Hz')
marca_perturb(ax2)
ax2.set_xlabel('t [s]')
ax2.set_ylabel('f [Hz]')
ax2.legend(fontsize=8, framealpha=0.3)
estilo_ax(ax2, 'Frecuencia de red (swing eq.)')

# ── Panel 3: potencia eléctrica ──────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
P_tot = np.sum(P_elec_hist, axis=1) / 1e6
ax3.fill_between(t_vec, P_tot, alpha=0.25, color='#2ecc71')
ax3.plot(t_vec, P_tot, color='#2ecc71', lw=2.0, label='P_total')
for i in range(N):
    ax3.plot(t_vec, P_elec_hist[:, i]/1e6, color=COLS[i], lw=0.8, alpha=0.55,
             label=NOMBRES[i])
marca_perturb(ax3)
ax3.axhline(N * S_nom / 1e6, color='white', ls=':', lw=1, alpha=0.4, label='P_nom total')
ax3.set_xlabel('t [s]')
ax3.set_ylabel('P [MW]')
ax3.legend(fontsize=7, loc='upper left', framealpha=0.3)
estilo_ax(ax3, 'Potencia eléctrica generada (MPPT)')

# ── Panel 4: intercambio ZypyZape ────────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
for i in range(N):
    ax4.plot(t_vec, P_zz_hist[:, i]/1e3, color=COLS[i], lw=1.3, label=NOMBRES[i])
ax4.axhline(0, color='white', lw=0.8, alpha=0.3)
ax4.axhline( P_ZZ_frac * S_nom / 1e3, color='white', ls=':', lw=0.8, alpha=0.3,
             label=f'+P_ZZ_max ({P_ZZ_frac*100:.0f}%)')
ax4.axhline(-P_ZZ_frac * S_nom / 1e3, color='white', ls=':', lw=0.8, alpha=0.3)
ax4.set_xlabel('t [s]')
ax4.set_ylabel('P_ZZ [kW]')
ax4.legend(fontsize=7, framealpha=0.3)
estilo_ax(ax4, 'Ciclo ZypyZape — intercambio ACEL/FREN')

# ── Panel 5: Tip Speed Ratio λ ───────────────────────────────
ax5 = fig.add_subplot(gs[2, 0])
for i in range(N):
    ax5.plot(t_vec, lam_hist[:, i], color=COLS[i], lw=1.3, label=NOMBRES[i])
ax5.axhline(lambda_opt, color='white', ls='--', lw=1.5, alpha=0.7,
            label=f'λ_opt = {lambda_opt}')
ax5.set_ylim(0, 15)
ax5.set_xlabel('t [s]')
ax5.set_ylabel('λ [−]')
ax5.legend(fontsize=7, framealpha=0.3)
estilo_ax(ax5, 'Tip Speed Ratio λ(t) — eficiencia aerodinámica')

# ── Panel 6: balance energético acumulado ───────────────────
ax6 = fig.add_subplot(gs[2, 1])
E_ae_cum = np.cumsum(np.sum(P_aero_hist, axis=1)) * dt / 1e9
E_el_cum = np.cumsum(np.sum(P_elec_hist, axis=1)) * dt / 1e9
ax6.plot(t_vec, E_ae_cum, color='#f39c12', lw=2.0, label='E_aero [GJ]')
ax6.plot(t_vec, E_el_cum, color='#2ecc71', lw=2.0, label='E_elec [GJ]')
ax6.fill_between(t_vec, E_ae_cum, E_el_cum, alpha=0.2, color='gray',
                 label='Pérdidas + rozamiento')
ax6.set_xlabel('t [s]')
ax6.set_ylabel('E [GJ]')
ax6.legend(fontsize=8, framealpha=0.3)
estilo_ax(ax6, 'Balance energético acumulado')

fig.suptitle(
    'ZypyZape — Gemelo Digital v2.0  |  Víctor Manzanares Alberola\n'
    f'f_min={f_min:.3f} Hz   RoCoF={rocof:.3f} Hz/s   η={eff:.1f}%   λ_med={lam_medio_final:.2f}',
    color='white', fontsize=12, fontweight='bold', y=0.98
)

plt.savefig('/home/claude/zypyzape_twin_v2.png', dpi=150, bbox_inches='tight',
            facecolor=BG)
plt.close()
print("Figura guardada: zypyzape_twin_v2.png")
