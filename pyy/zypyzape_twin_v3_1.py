"""
ZypyZape Digital Twin v3.0
==========================
Sistema: Víctor Manzanares Alberola
Asistencia de código: IA

Novedades v3 sobre v2:
  - Gobernador primario (droop R=5%, T_gov=5s): equilibrio de frecuencia real
  - Sistema escalado a 2 GW: ZypyZape (12.5 MW) es ~0.6%, contribución visible
  - Nadir sin clip artificial (min 47 Hz, rango fisico real)
  - H_eff(t): tracking en tiempo real de la inercia aportada
  - Analisis de sensibilidad: f_min y RoCoF vs k_droop
  - Panel Cp(lambda): verificacion curva aerodinamica
  - Exportacion CSV completa
  - Animacion GIF: barras omega + osciloscopio f(t) + flujo P_ZZ

Modelo de red (swing eq. + gobernador primario):
  (2H/w0)*d2delta/dt2 = dP_mec - dP_elec          [swing]
  T_gov * dP_gov/dt   = -(f-f0)/(f0*R) - P_gov    [gobernador]
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import csv, os

# ============================================================
# PARAMETROS FISICOS -- turbina 2.5 MW HAWT
# ============================================================
N         = 5
R         = 60.0
A         = np.pi * R**2
rho       = 1.225
J         = 5e6
S_nom     = 2.5e6
Kv        = 2e4
T_c0      = 1e3
omega_min = 0.15
omega_max = 2.8

lambda_opt  = 7.0
Cp_max      = 0.48
k_mppt      = 0.5 * rho * A * R**3 * Cp_max / (lambda_opt**3)
v_rated     = (S_nom / (0.5 * rho * A * Cp_max)) ** (1.0/3.0)
omega_rated = lambda_opt * v_rated / R
T_rated     = S_nom / omega_rated
K_pitch     = 15e6

P_ZZ_frac  = 0.18
f_cycle    = 0.4
T_cycle    = 1.0 / f_cycle

Kc         = 1.5e5

# ============================================================
# RED ELECTRICA -- sistema 2 GW
# ============================================================
S_total   = 2e9
H_sistema = 4.0
f0        = 50.0
D_amort   = 0.05

T_gov     = 5.0
R_gov     = 0.05
P_gov_max = 0.15 * S_total

# ============================================================
# TIEMPO Y EVENTOS
# ============================================================
dt         = 0.02
T_sim      = 90.0
steps      = int(T_sim / dt)
t_perturb  = 30.0
dP_perturb = -100e6

# ============================================================
# TOPOLOGIA
# ============================================================
K_mat = np.zeros((N, N))
K_mat[0,1] = K_mat[1,0] = Kc
K_mat[2,3] = K_mat[3,2] = Kc
K_mat[3,4] = K_mat[4,3] = Kc
K_mat[4,2] = K_mat[2,4] = Kc
PARES_ZZ = [(0, 1), (2, 3)]

# ============================================================
# FISICA
# ============================================================
def viento(t):
    v = 11.5 + 1.5*np.sin(0.08*t) + 0.7*np.sin(0.25*t + 0.9)
    if 45.0 < t < 55.0:
        v += 3.5 * np.sin(np.pi*(t - 45.0)/10.0)
    return max(3.0, float(v))

def Cp_curva(lmbda):
    lmbda = float(lmbda)
    if lmbda < 1.0:
        return 0.0
    s = 3.2 if lmbda <= lambda_opt else 1.8
    return float(Cp_max * np.exp(-0.5 * ((lmbda - lambda_opt) / s)**2))

def T_aero(omega_i, v):
    if omega_i < 0.05:
        return 0.0
    return 0.5 * rho * A * v**3 * Cp_curva(omega_i * R / v) / omega_i

def T_roz(omega_i):
    return Kv * omega_i + T_c0

def T_gen_mppt(omega_i):
    if omega_i <= omega_rated:
        return k_mppt * omega_i**2
    return min(T_rated + K_pitch * (omega_i - omega_rated), T_rated * 3.5)

def T_zypyzape(i, t, omega_i):
    T_zz = 0.0
    for idx, (a, b) in enumerate(PARES_ZZ):
        if i != a and i != b:
            continue
        fase = (t / T_cycle + idx * 0.5) % 1.0
        role = +1 if (i == a) == (fase < 0.5) else -1
        T_zz += role * P_ZZ_frac * S_nom * abs(np.sin(np.pi * f_cycle * t)) / max(omega_i, 0.1)
    return T_zz

def T_inercia_sintetica(f_grid, omega_i, kd):
    df = f0 - f_grid
    if df <= 0.0:
        return 0.0
    return -kd * (df / f0) * S_nom / max(omega_i, 0.1)

# ============================================================
# SIMULACION REUTILIZABLE
# ============================================================
def simular(mode='ZZ', kd=0.06):
    omega  = np.ones(N) * 1.25
    theta  = np.zeros(N)
    f_grid = 50.0
    P_gov  = 0.0

    hs = {k: np.zeros((steps, N)) if k in ('omega','P_aero','P_elec','P_zz','lam')
          else np.zeros(steps)
          for k in ('omega','P_aero','P_elec','P_zz','lam','f','v','H_eff','P_gov')}
    E_aero = E_elec = 0.0

    for step in range(steps):
        t  = step * dt
        v  = viento(t)
        dP = dP_perturb if t >= t_perturb else 0.0

        # Gobernador primario
        dPg  = (-(f_grid - f0) / (f0 * R_gov) - P_gov) / T_gov
        P_gov = np.clip(P_gov + dPg * dt, 0.0, P_gov_max)

        # Inercia efectiva
        H_zz  = (0.5 * J * float(np.sum(omega**2)) / S_total) if mode == 'ZZ' else 0.0
        H_eff = H_sistema + H_zz

        # Swing equation
        df_dt  = ((dP + P_gov) * f0 / (2.0 * H_eff * S_total)
                  - D_amort * (f_grid - f0))
        f_grid = np.clip(f_grid + df_dt * dt, 47.0, 52.0)

        domega = np.zeros(N)
        for i in range(N):
            ta = T_aero(omega[i], v)
            tg = T_gen_mppt(omega[i])
            tr = T_roz(omega[i])
            tzz = T_zypyzape(i, t, omega[i]) if mode == 'ZZ' else 0.0
            tis = T_inercia_sintetica(f_grid, omega[i], kd) if mode == 'ZZ' else 0.0
            coup = sum(K_mat[i,j]*np.sin(theta[j]-theta[i])
                       for j in range(N) if K_mat[i,j] != 0)
            domega[i] = (ta - tg - tr + tzz + tis + coup) / J

            hs['P_aero'][step,i] = ta * omega[i]
            hs['P_elec'][step,i] = tg * omega[i]
            hs['P_zz'][step,i]   = tzz * omega[i]
            hs['lam'][step,i]    = omega[i] * R / v
            E_aero += ta * omega[i] * dt
            E_elec += tg * omega[i] * dt

        omega = np.clip(omega + domega * dt, omega_min, omega_max)
        theta += omega * dt
        hs['omega'][step] = omega
        hs['f'][step]     = f_grid
        hs['v'][step]     = v
        hs['H_eff'][step] = H_eff
        hs['P_gov'][step] = P_gov

    hs['E_aero'] = E_aero
    hs['E_elec'] = E_elec
    return hs

# ============================================================
# SIMULACIONES PRINCIPALES
# ============================================================
print("Simulando CON ZypyZape...")
hz = simular('ZZ', kd=0.06)
print("Simulando SIN ZypyZape (referencia)...")
hr = simular('REF', kd=0.0)

t_vec = np.arange(steps) * dt
idx_p = int(t_perturb / dt)

f_min_zz  = float(np.min(hz['f'][idx_p:]))
f_min_ref = float(np.min(hr['f'][idx_p:]))
rocof_zz  = float((hz['f'][idx_p+1] - hz['f'][idx_p]) / dt)
rocof_ref = float((hr['f'][idx_p+1] - hr['f'][idx_p]) / dt)
eff_zz    = 100.0 * hz['E_elec'] / max(hz['E_aero'], 1.0)
lam_final = float(np.mean(hz['lam'][-100:]))

print()
print("=" * 62)
print("   ZYPYZAPE GEMELO DIGITAL v3.0")
print("=" * 62)
print(f"   E_aero          : {hz['E_aero']/1e9:.3f} GJ")
print(f"   E_elec          : {hz['E_elec']/1e9:.3f} GJ")
print(f"   Eficiencia      : {eff_zz:.1f} %")
print(f"   lambda_medio    : {lam_final:.2f}  (opt={lambda_opt})")
print(f"   omega_medio fin : {float(np.mean(hz['omega'][-1])):.3f} rad/s")
print(f"   H_eff medio     : {float(np.mean(hz['H_eff'])):.5f} s")
print(f"   P_gov final     : {float(hz['P_gov'][-1])/1e6:.3f} MW")
print()
print(f"   {'Metrica':<24} {'SIN ZypyZape':>15} {'CON ZypyZape':>15}")
print(f"   {'-'*54}")
print(f"   {'f_min (nadir)':<24} {f_min_ref:>14.5f} Hz {f_min_zz:>13.5f} Hz")
print(f"   {'RoCoF inicial':<24} {rocof_ref:>14.5f} Hz/s {rocof_zz:>10.5f} Hz/s")
if abs(rocof_ref) > 1e-9:
    print(f"   {'Mejora RoCoF':<24} {(abs(rocof_ref)-abs(rocof_zz))/abs(rocof_ref)*100:>14.2f} %")
print(f"   {'Mejora nadir':<24} {f_min_zz - f_min_ref:>+14.5f} Hz")
print("=" * 62)

# ============================================================
# SENSIBILIDAD k_droop
# ============================================================
print("\nSensibilidad k_droop (18 puntos)...")
k_vals = np.linspace(0.0, 0.18, 18)
fmin_s, rocof_s = [], []
for kd in k_vals:
    hs_s = simular('ZZ', kd=kd)
    fmin_s.append(float(np.min(hs_s['f'][idx_p:])))
    rocof_s.append(float(abs((hs_s['f'][idx_p+1] - hs_s['f'][idx_p]) / dt)))

# ============================================================
# CSV
# ============================================================
csv_path = '/home/claude/zypyzape_v3_datos.csv'
with open(csv_path, 'w', newline='') as cf:
    w = csv.writer(cf)
    hdr = ['t_s','v_ms','f_ZZ_Hz','f_REF_Hz','H_eff_s','P_gov_MW']
    for i in range(N):
        hdr += [f'omega{i}_rads',f'Pelec{i}_MW',f'Pzz{i}_kW',f'lam{i}']
    w.writerow(hdr)
    for s in range(steps):
        row = [round(s*dt,3), round(float(hz['v'][s]),3),
               round(float(hz['f'][s]),5), round(float(hr['f'][s]),5),
               round(float(hz['H_eff'][s]),6), round(float(hz['P_gov'][s])/1e6,4)]
        for i in range(N):
            row += [round(float(hz['omega'][s,i]),4),
                    round(float(hz['P_elec'][s,i])/1e6,4),
                    round(float(hz['P_zz'][s,i])/1e3,3),
                    round(float(hz['lam'][s,i]),3)]
        w.writerow(row)
print(f"CSV: {csv_path}  ({steps} filas)")

# ============================================================
# FIGURA ESTATICA 8 PANELES
# ============================================================
NOMBRES = ['Central-A','Central-B','Anillo-2','Anillo-3','Anillo-4']
COLS    = ['#e74c3c','#3498db','#2ecc71','#f39c12','#9b59b6']
BG, PAN = '#0d0d1a', '#13132b'

def estilar(ax, titulo, xlabel='t [s]'):
    ax.set_facecolor(PAN)
    ax.set_title(titulo, color='white', fontsize=9.5, pad=5)
    ax.tick_params(colors='#aaaaaa', labelsize=7.5)
    for sp in ax.spines.values(): sp.set_color('#333355')
    ax.set_xlabel(xlabel, color='#aaaaaa', fontsize=8)
    ax.yaxis.label.set_color('#aaaaaa')
    ax.grid(color='#1e1e40', lw=0.6, ls='--')

def vp(ax): ax.axvline(t_perturb, color='#ffcc00', ls='--', lw=1.1, alpha=0.65)

fig = plt.figure(figsize=(18, 14), facecolor=BG)
gs  = gridspec.GridSpec(4, 2, figure=fig, hspace=0.58, wspace=0.38)

# P1 velocidades angulares
ax1 = fig.add_subplot(gs[0,0])
for i in range(N):
    ax1.plot(t_vec, hz['omega'][:,i], color=COLS[i], lw=1.6, label=NOMBRES[i])
ax1.axhline(omega_rated, color='white', ls=':', lw=0.9, alpha=0.4,
            label=f'omega_rat={omega_rated:.2f}')
vp(ax1)
ax1.set_ylabel('omega [rad/s]')
ax1.legend(fontsize=7, framealpha=0.25, ncol=2)
estilar(ax1, 'Velocidades angulares omega(t)')

# P2 frecuencia comparativa
ax2 = fig.add_subplot(gs[0,1])
ax2.plot(t_vec, hr['f'], color='#e74c3c', lw=1.8, ls='--', alpha=0.85,
         label=f'Sin ZZ  nadir={f_min_ref:.4f} Hz')
ax2.plot(t_vec, hz['f'], color='#00d2ff', lw=2.3,
         label=f'Con ZZ  nadir={f_min_zz:.4f} Hz')
ax2.axhline(49.0, color='#ffcc00', ls=':', lw=1.2, alpha=0.7, label='49 Hz ENTSO-E')
ax2.axhline(50.0, color='white', ls='--', lw=0.7, alpha=0.2)
vp(ax2)
if f_min_zz > f_min_ref:
    ax2.annotate(f'+{f_min_zz-f_min_ref:.4f} Hz',
                 xy=(t_perturb+5, (f_min_zz+f_min_ref)/2),
                 fontsize=9, color='#2ecc71',
                 arrowprops=dict(arrowstyle='->', color='#2ecc71', lw=1.2))
ax2.set_ylabel('f [Hz]')
ax2.legend(fontsize=7.5, framealpha=0.25)
estilar(ax2, 'Frecuencia de red -- sin vs. con ZypyZape')

# P3 potencia electrica + P_gov
ax3 = fig.add_subplot(gs[1,0])
Pt = np.sum(hz['P_elec'], axis=1)/1e6
ax3.fill_between(t_vec, Pt, alpha=0.15, color='#2ecc71')
ax3.plot(t_vec, np.sum(hr['P_elec'],axis=1)/1e6, color='#e74c3c', lw=1.3,
         ls='--', alpha=0.7, label='Ref')
ax3.plot(t_vec, Pt, color='#2ecc71', lw=2.0, label='ZypyZape total')
ax3.plot(t_vec, hz['P_gov']/1e6, color='#f39c12', lw=1.5, ls='-.', label='P_gov')
for i in range(N):
    ax3.plot(t_vec, hz['P_elec'][:,i]/1e6, color=COLS[i], lw=0.7, alpha=0.3)
ax3.axhline(N*S_nom/1e6, color='white', ls=':', lw=0.8, alpha=0.3, label='P_nom')
vp(ax3)
ax3.set_ylabel('P [MW]')
ax3.legend(fontsize=7.5, framealpha=0.25)
estilar(ax3, 'Potencia electrica + respuesta gobernador')

# P4 ciclo ZypyZape
ax4 = fig.add_subplot(gs[1,1])
for i in range(N):
    ax4.plot(t_vec, hz['P_zz'][:,i]/1e3, color=COLS[i], lw=1.3, label=NOMBRES[i])
ax4.axhline(0, color='white', lw=0.7, alpha=0.25)
Pm = P_ZZ_frac*S_nom/1e3
ax4.axhline(Pm,  color='#aaaaff', ls=':', lw=0.8, alpha=0.5)
ax4.axhline(-Pm, color='#aaaaff', ls=':', lw=0.8, alpha=0.5)
ax4.set_ylabel('P_ZZ [kW]')
ax4.legend(fontsize=7, framealpha=0.25, ncol=2)
estilar(ax4, 'Ciclo ZypyZape -- intercambio ACEL/FREN')

# P5 H_eff(t)
ax5 = fig.add_subplot(gs[2,0])
ax5.plot(t_vec, hz['H_eff'], color='#f39c12', lw=1.8, label='H_eff(t)')
ax5.axhline(H_sistema, color='white', ls='--', lw=1.0, alpha=0.5,
            label=f'H_sistema={H_sistema} s')
ax5.fill_between(t_vec, H_sistema, hz['H_eff'], alpha=0.3, color='#f39c12',
                 label='DeltaH ZypyZape')
vp(ax5)
ax5.set_ylabel('H [s]')
ax5.legend(fontsize=7.5, framealpha=0.25)
estilar(ax5, 'Inercia efectiva H_eff(t) aportada al sistema')

# P6 Cp(lambda) + puntos operativos
ax6 = fig.add_subplot(gs[2,1])
lr = np.linspace(0.5, 14, 300)
ax6.plot(lr, [Cp_curva(l) for l in lr], color='#00d2ff', lw=2.2, label='Cp(lambda)')
ax6.axvline(lambda_opt, color='#f39c12', ls='--', lw=1.2, label=f'lambda_opt={lambda_opt}')
ax6.axhline(Cp_max, color='#2ecc71', ls=':', lw=1.0, alpha=0.7, label=f'Cp_max={Cp_max}')
ax6.axhline(0.593, color='white', ls=':', lw=0.8, alpha=0.3, label='Betz 0.593')
for i in range(N):
    lm = float(np.mean(hz['lam'][-100:,i]))
    ax6.scatter(lm, Cp_curva(lm), color=COLS[i], s=60, zorder=5,
                label=f'{NOMBRES[i]} lam={lm:.1f}')
ax6.set_xlim(0,14); ax6.set_ylim(0, 0.62)
ax6.set_xlabel('lambda [-]', color='#aaaaaa', fontsize=8)
ax6.set_ylabel('Cp [-]')
ax6.legend(fontsize=6.5, framealpha=0.25, ncol=2)
estilar(ax6, 'Curva Cp(lambda) -- puntos = estado final turbinas', xlabel='lambda [-]')

# P7 balance energetico
ax7 = fig.add_subplot(gs[3,0])
Eae = np.cumsum(np.sum(hz['P_aero'],axis=1))*dt/1e9
Eel = np.cumsum(np.sum(hz['P_elec'],axis=1))*dt/1e9
ax7.plot(t_vec, Eae, color='#f39c12', lw=2.0, label='E_aero [GJ]')
ax7.plot(t_vec, Eel, color='#2ecc71', lw=2.0, label='E_elec [GJ]')
ax7.fill_between(t_vec, Eae, Eel, alpha=0.2, color='gray', label='Perdidas')
ax7.set_ylabel('E [GJ]')
ax7.legend(fontsize=8, framealpha=0.25)
estilar(ax7, 'Balance energetico acumulado')

# P8 sensibilidad f_min y RoCoF vs k_droop
ax8 = fig.add_subplot(gs[3,1])
ax8b = ax8.twinx()
ax8.plot(k_vals*100, fmin_s, color='#9b59b6', lw=2.2, marker='o', ms=4, label='f_min [Hz]')
ax8b.plot(k_vals*100, rocof_s, color='#e74c3c', lw=1.8, marker='s', ms=3,
          ls='--', label='|RoCoF| [Hz/s]')
ax8.axhline(f_min_ref, color='#e74c3c', ls=':', lw=1.2, alpha=0.7,
            label=f'f_min sin ZZ ({f_min_ref:.3f} Hz)')
ax8.axhline(49.0, color='#ffcc00', ls=':', lw=1.0, alpha=0.6, label='49 Hz')
ax8.axvline(0.06*100, color='white', ls=':', lw=0.9, alpha=0.4, label='k nominal')
ax8.set_xlabel('k_droop [%]', color='#aaaaaa', fontsize=8)
ax8.set_ylabel('f_min [Hz]', color='#9b59b6')
ax8b.set_ylabel('|RoCoF| [Hz/s]', color='#e74c3c')
ax8b.tick_params(colors='#e74c3c', labelsize=7.5)
ax8b.set_facecolor(PAN)
l1,la1 = ax8.get_legend_handles_labels()
l2,la2 = ax8b.get_legend_handles_labels()
ax8.legend(l1+l2, la1+la2, fontsize=7, framealpha=0.25)
ax8.set_facecolor(PAN)
ax8.tick_params(colors='#aaaaaa', labelsize=7.5)
for sp in ax8.spines.values(): sp.set_color('#333355')
ax8.grid(color='#1e1e40', lw=0.6, ls='--')
ax8.set_title('Sensibilidad: f_min y |RoCoF| vs. k_droop', color='white', fontsize=9.5, pad=5)

fig.suptitle(
    'ZypyZape -- Gemelo Digital v3.0  |  Victor Manzanares Alberola\n'
    f'nadir: {f_min_ref:.4f} -> {f_min_zz:.4f} Hz  (+{f_min_zz-f_min_ref:.4f} Hz)  |  '
    f'RoCoF: {rocof_ref:.4f} -> {rocof_zz:.4f} Hz/s  |  '
    f'eta={eff_zz:.1f}%  |  S=2GW  |  dP=-100MW',
    color='white', fontsize=11, fontweight='bold', y=0.997
)
png_path = '/home/claude/zypyzape_v3_estatica.png'
plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor=BG)
plt.close()
print(f"Figura estatica: {png_path}")

# ============================================================
# ANIMACION GIF
# ============================================================
print("Generando animacion GIF...")

NFRAMES        = 270
SPF            = steps // NFRAMES
win_steps_anim = int(15.0 / dt)

fig_a = plt.figure(figsize=(15, 7), facecolor=BG)
gsa   = gridspec.GridSpec(2, 3, figure=fig_a, hspace=0.52, wspace=0.42)

def init_ax(ax, titulo):
    ax.set_facecolor(PAN)
    ax.set_title(titulo, color='white', fontsize=9)
    ax.tick_params(colors='#aaaaaa', labelsize=7)
    for sp in ax.spines.values(): sp.set_color('#333355')
    ax.grid(color='#1e1e40', lw=0.6, ls='--')

ax_bar = fig_a.add_subplot(gsa[0,0])
init_ax(ax_bar, 'omega rotor [rad/s]')
ax_bar.set_xlim(-0.5, N-0.5); ax_bar.set_ylim(0, omega_max*1.05)
ax_bar.axhline(omega_rated, color='#ffcc00', ls='--', lw=0.9, alpha=0.7)
ax_bar.set_xticks(range(N))
ax_bar.set_xticklabels(['CA','CB','A2','A3','A4'], color='#aaaaaa', fontsize=8)
bars_w = ax_bar.bar(range(N), hz['omega'][0], color=COLS, width=0.65, alpha=0.85)

ax_osc = fig_a.add_subplot(gsa[0,1:])
init_ax(ax_osc, 'Frecuencia de red f(t) [Hz]')
ax_osc.set_xlim(0, 15); ax_osc.set_ylim(49.2, 50.3)
ax_osc.axhline(50.0, color='white', lw=0.7, alpha=0.2)
ax_osc.axhline(49.5, color='#ffcc00', ls=':', lw=1.0, alpha=0.6)
ax_osc.set_xlabel('Dt [s]', color='#aaaaaa', fontsize=8)
ax_osc.set_ylabel('f [Hz]', color='#aaaaaa', fontsize=8)
ln_zz,  = ax_osc.plot([], [], color='#00d2ff', lw=2.1, label='Con ZypyZape')
ln_ref, = ax_osc.plot([], [], color='#e74c3c', lw=1.5, ls='--', alpha=0.75,
                       label='Sin ZypyZape')
ax_osc.legend(fontsize=7.5, framealpha=0.3, loc='lower right')

ax_pzz = fig_a.add_subplot(gsa[1,0])
init_ax(ax_pzz, 'Ciclo ZypyZape [kW]')
ax_pzz.set_xlim(-0.5, N-0.5)
Pm_bar = P_ZZ_frac*S_nom/1e3*1.15
ax_pzz.set_ylim(-Pm_bar, Pm_bar)
ax_pzz.axhline(0, color='white', lw=0.7, alpha=0.25)
ax_pzz.set_xticks(range(N))
ax_pzz.set_xticklabels(['CA','CB','A2','A3','A4'], color='#aaaaaa', fontsize=8)
bars_zz = ax_pzz.bar(range(N), [0]*N, color=COLS, width=0.65, alpha=0.85)

ax_pel = fig_a.add_subplot(gsa[1,1:])
init_ax(ax_pel, 'Potencia electrica por turbina [MW]')
ax_pel.set_xlim(0, 15); ax_pel.set_ylim(0, S_nom/1e6*1.25)
ax_pel.set_xlabel('Dt [s]', color='#aaaaaa', fontsize=8)
ax_pel.set_ylabel('P [MW]', color='#aaaaaa', fontsize=8)
lns_p = [ax_pel.plot([], [], color=COLS[i], lw=1.5, label=NOMBRES[i])[0]
         for i in range(N)]
ax_pel.legend(fontsize=7, framealpha=0.25, ncol=3)

txt_top = fig_a.text(0.5, 0.97, '', ha='center', va='top',
                     color='white', fontsize=10, fontweight='bold')

def init_an():
    ln_zz.set_data([],[]); ln_ref.set_data([],[])
    for ln in lns_p: ln.set_data([],[])
    return [ln_zz,ln_ref]+lns_p+list(bars_w)+list(bars_zz)

def update_an(frame):
    s   = min(frame*SPF, steps-1)
    t   = s*dt
    i0  = max(0, s-win_steps_anim)
    tw  = t_vec[i0:s+1] - t_vec[i0]
    for i,b in enumerate(bars_w):
        b.set_height(float(hz['omega'][s,i]))
    for i,b in enumerate(bars_zz):
        val = float(hz['P_zz'][s,i])/1e3
        b.set_height(abs(val)); b.set_y(min(0.0, val))
    ln_zz.set_data(tw,  hz['f'][i0:s+1])
    ln_ref.set_data(tw, hr['f'][i0:s+1])
    for i,ln in enumerate(lns_p):
        ln.set_data(tw, hz['P_elec'][i0:s+1,i]/1e6)
    perturb = '  !! dP=-100MW' if t >= t_perturb else ''
    txt_top.set_text(
        f't={t:.1f}s  v={float(hz["v"][s]):.1f}m/s  '
        f'f={float(hz["f"][s]):.4f}Hz  '
        f'H_eff={float(hz["H_eff"][s]):.5f}s  '
        f'P_gov={float(hz["P_gov"][s])/1e6:.1f}MW{perturb}'
    )
    return [ln_zz,ln_ref]+lns_p+list(bars_w)+list(bars_zz)+[txt_top]

anim_obj = animation.FuncAnimation(
    fig_a, update_an, init_func=init_an,
    frames=NFRAMES, interval=40, blit=False
)
gif_path = '/home/claude/zypyzape_v3_animacion.gif'
anim_obj.save(gif_path, writer=animation.PillowWriter(fps=25), dpi=90)
plt.close()
print(f"Animacion GIF: {gif_path}")

print()
print("Archivos generados:")
for p in [png_path, gif_path, csv_path]:
    print(f"  {p}   ({os.path.getsize(p)//1024} KB)")
