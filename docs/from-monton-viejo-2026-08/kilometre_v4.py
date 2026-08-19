"""
KILÒMETRE — Simulació v4
Reescriptura completa de la comptabilitat energètica.

REGLA ÚNICA: tota l'energia ve de dW = tau_i * omega * dt
  - dW_grav  = tau_g   * omega * dt   (pot ser + o -)
  - dW_fric  = tau_fr  * omega * dt   (sempre ≤ 0, pèrdua)
  - dW_ctrl  = tau_c   * omega * dt   (+ si motor aporta, - si gen extreu)

  Si mode == MOTOR:
    dW_ctrl > 0  → energia mecànica injectada
    energia elèctrica consumida = dW_ctrl / eta_mot

  Si mode == GEN:
    dW_ctrl < 0  → energia mecànica extreta del sistema
    energia elèctrica produïda  = |dW_ctrl| * eta_gen

Primera Llei (verificació):
  ΔKE = dW_grav + dW_ctrl + dW_fric
  → residu = |ΔKE_calculada - ΔKE_real| < 0.1%

Política de control:
  ph_mod ∈ [0°, 90°)    → MOTOR   (zona desfavorable, passar ràpid)
  ph_mod ∈ [90°, 180°)  → LLIURE  (inercia, zero control)
  ph_mod ∈ [180°, 360°) → GEN     (zona favorable, embriagar sinfín)
    Intensitat de frenada ∝ v_obj (velocitat axial de l'objecte)
    Llindar mínim: no frena si omega < omega_min
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ════════════════════════════════════════════════════════════
#  PARÀMETRES  — ajusta aquí
# ════════════════════════════════════════════════════════════
m         = 10.0    # kg    massa objecte
R         = 0.15    # m     radi objecte-eix
g         = 9.81    # m/s²
pas       = 0.04    # m/rad avanç sinfín per radian

m_tub     = 2.0     # kg
R_tub     = 0.18    # m
I         = 0.5*m_tub*R_tub**2 + m*R**2   # [kg·m²]

beta      = 0.02    # N·m·s/rad   fricció viscosa

omega_0   = 3.0     # rad/s  velocitat inicial
omega_min = 0.5     # rad/s  llindar mínim per activar generador

# Control motor (zona 0°→90°):
#   tau_motor = alfa_m * tau_grav_max
alfa_motor = 1.10   # sobrecompensació 10% sobre el màxim gravitatori

# Control generador (zona 180°→360°):
#   tau_gen = -alfa_g * tau_grav   (frena una fracció del par gravitatori)
alfa_gen  = 0.70    # fracció de frenada

eta_gen   = 0.90
eta_mot   = 0.92

# Simulació
dt        = 0.0001  # s
t_max     = 10.0    # s
N         = int(t_max / dt)

# ════════════════════════════════════════════════════════════
#  BUCLE PRINCIPAL
# ════════════════════════════════════════════════════════════
tau_grav_max = m * g * R

phi    = np.empty(N); phi[0]   = 0.01
omega  = np.empty(N); omega[0] = omega_0

# Acumuladors energètics (pas a pas, s'agreguen al final)
dW_grav_v  = np.zeros(N)
dW_fric_v  = np.zeros(N)
dW_ctrl_v  = np.zeros(N)   # negatiu si gen, positiu si motor
dW_elec_in = np.zeros(N)   # energia elèctrica consumida (motor)
dW_elec_out= np.zeros(N)   # energia elèctrica produïda  (gen)

mode_v     = np.zeros(N, dtype=int)
tau_g_v    = np.zeros(N)
tau_c_v    = np.zeros(N)

for i in range(1, N):
    ph = phi[i-1];  om = omega[i-1]
    ph_mod = ph % (2*np.pi)

    # Parell gravitatori
    tau_g = -m * g * R * np.cos(ph)

    # Fricció
    tau_fr = -beta * om

    # Control per fase
    if ph_mod < np.pi/2:
        # MOTOR
        tau_c = alfa_motor * tau_grav_max
        mode  = 1
    elif ph_mod < np.pi:
        # LLIURE
        tau_c = 0.0
        mode  = 0
    else:
        # GENERADOR (si va prou de pressa)
        if abs(om) > omega_min:
            # Frena una fracció del parell gravitatori disponible
            # En aquesta zona cos(ph) < 0, tau_g > 0 (ajuda el gir)
            # El generador frena → tau_c negatiu
            tau_c = -alfa_gen * abs(tau_g)
            mode  = 2
        else:
            tau_c = 0.0
            mode  = 0

    # Dinàmica
    alpha    = (tau_g + tau_c + tau_fr) / I
    om_new   = om + alpha * dt
    ph_new   = ph + om * dt

    phi[i]   = ph_new
    omega[i] = om_new

    # ── Comptabilitat (regla única: tau * omega * dt) ────
    dWg  = tau_g  * om * dt
    dWfr = tau_fr * om * dt      # ≤ 0
    dWc  = tau_c  * om * dt      # + si motor, - si gen

    dW_grav_v[i] = dWg
    dW_fric_v[i] = dWfr
    dW_ctrl_v[i] = dWc

    if mode == 1 and dWc > 0:
        # Motor aporta dWc de mecànica → consumeix dWc/eta_mot d'elèctrica
        dW_elec_in[i]  = dWc / eta_mot
        dW_elec_out[i] = 0.0
    elif mode == 2 and dWc < 0:
        # Gen extreu |dWc| de mecànica → produeix |dWc|*eta_gen d'elèctrica
        dW_elec_in[i]  = 0.0
        dW_elec_out[i] = abs(dWc) * eta_gen
    else:
        dW_elec_in[i]  = 0.0
        dW_elec_out[i] = 0.0

    mode_v[i]  = mode
    tau_g_v[i] = tau_g
    tau_c_v[i] = tau_c

# ════════════════════════════════════════════════════════════
#  VERIFICACIÓ PRIMERA LLEI (pas a pas)
#  ΔKE[i] = dW_grav + dW_ctrl + dW_fric
# ════════════════════════════════════════════════════════════
KE        = 0.5 * I * omega**2
dKE_real  = np.diff(KE)
dKE_calc  = dW_grav_v[1:] + dW_ctrl_v[1:] + dW_fric_v[1:]
residu_v  = np.abs(dKE_real - dKE_calc)
residu_rel= residu_v / (np.abs(dKE_real) + 1e-12)

max_residu = np.max(residu_rel[np.abs(dKE_real) > 1e-8])
mean_residu= np.mean(residu_rel[np.abs(dKE_real) > 1e-8])

# Totals
W_grav    = np.sum(dW_grav_v)
W_fric    = np.sum(np.abs(dW_fric_v))
W_ctrl    = np.sum(dW_ctrl_v)
W_el_in   = np.sum(dW_elec_in)
W_el_out  = np.sum(dW_elec_out)
dKE_total = KE[-1] - KE[0]
W_net_el  = W_el_out - W_el_in

n_voltes  = phi[-1] / (2*np.pi)

# Residu global
residu_global     = abs(dKE_total - (W_grav + W_ctrl + (-W_fric)))
escala_global     = max(abs(dKE_total), abs(W_grav), abs(W_el_in), 0.001)

print("=" * 62)
print("KILÒMETRE v4 — Balanç Energètic (comptabilitat corregida)")
print("=" * 62)
print(f"I = {I:.4f} kg·m²  |  tau_grav_max = {tau_grav_max:.3f} N·m")
print(f"Temps: {t_max}s  |  Voltes: {n_voltes:.2f}")
print(f"omega: {omega[0]:.3f} → {omega[-1]:.3f} rad/s  |  ΔKE = {dKE_total:+.4f} J")
print("-" * 62)
print(f"W_grav   (net gravetat):      {W_grav:+10.4f} J  ← ~0 cicle tancat")
print(f"W_ctrl   (net control mec.):  {W_ctrl:+10.4f} J  (+ motor, - gen)")
print(f"W_fric   (pèrdues mec.):     {-W_fric:+10.4f} J")
print(f"  Check: W_grav+W_ctrl-W_fric = {W_grav+W_ctrl-W_fric:+.4f} J  vs ΔKE = {dKE_total:+.4f} J")
print("-" * 62)
print(f"W_el_in  (elèctric consumit): {W_el_in:+10.4f} J  (motor, paret)")
print(f"W_el_out (elèctric produït):  {W_el_out:+10.4f} J  (generador)")
print(f"BALANÇ NET ELÈCTRIC:          {W_net_el:+10.4f} J")
print("-" * 62)
print(f"Residu 1a Llei (màx pas):     {max_residu*100:.4f}%")
print(f"Residu 1a Llei (mitjà pas):   {mean_residu*100:.6f}%")
print(f"Residu global:                {residu_global:.4f} J  ({100*residu_global/escala_global:.3f}%)")
print("-" * 62)
if W_net_el > 0.01:
    print("⚠  BALANÇ ELÈCTRIC POSITIU — comprova física")
else:
    print(f"✓  Pèrdua neta: {abs(W_net_el):.4f} J  ({abs(W_net_el)/t_max*1000:.2f} mW)")
    print(f"   Explicació: el motor paga les pèrdues de fricció + conversió")
print("=" * 62)

# Anàlisi per volta
print("\nAnàlisi per volta (mecànica):")
print(f"{'V':>3} | {'W_grav':>8} | {'W_ctrl':>8} | {'W_fric':>8} | "
      f"{'ΔKE_real':>9} | {'residu%':>8}")
print("-" * 58)
prev = 0
vc   = 0
for i in range(1, N):
    if int(phi[i] / (2*np.pi)) > int(phi[i-1] / (2*np.pi)):
        s   = slice(prev, i)
        wg  = np.sum(dW_grav_v[s])
        wc  = np.sum(dW_ctrl_v[s])
        wf  = np.sum(np.abs(dW_fric_v[s]))
        dke = KE[i] - KE[prev]
        res = abs(wg + wc - wf - dke) / (abs(dke)+abs(wg)+1e-6) * 100
        print(f"{vc+1:>3} | {wg:>8.4f} | {wc:>8.4f} | {wf:>8.4f} | "
              f"{dke:>9.4f} | {res:>8.4f}%")
        prev = i;  vc += 1
        if vc >= 8: break

# ════════════════════════════════════════════════════════════
#  GRÀFICS
# ════════════════════════════════════════════════════════════
t = np.linspace(0, t_max, N)
fig, axes = plt.subplots(4, 1, figsize=(12, 14))
fig.suptitle("KILÒMETRE v4 — Control angle+velocitat objecte\n"
             "Comptabilitat energètica verificada (1a Llei)",
             fontsize=12, fontweight='bold')

# Ombres de mode (mostra primer 4s per claredat)
def ombres(ax, t, mode_v, fin=4.0):
    mask_m = (mode_v == 1) & (t <= fin)
    mask_g = (mode_v == 2) & (t <= fin)
    for j in np.where(np.diff(mask_m.astype(int)) != 0)[0]:
        pass
    # Ombra per blocs
    in_m = False; in_g = False; t0 = 0
    for j in range(len(t)):
        if t[j] > fin: break
        if mode_v[j] == 1 and not in_m:
            t0 = t[j]; in_m = True
        elif mode_v[j] != 1 and in_m:
            ax.axvspan(t0, t[j], alpha=0.12, color='red', lw=0)
            in_m = False
        if mode_v[j] == 2 and not in_g:
            t0 = t[j]; in_g = True
        elif mode_v[j] != 2 and in_g:
            ax.axvspan(t0, t[j], alpha=0.12, color='green', lw=0)
            in_g = False

# 1. Velocitat angular
ax = axes[0]
ax.plot(t, omega, color='steelblue', lw=1.2)
ax.axhline(omega_0, color='red', lw=1, ls='--', label=f'ω₀={omega_0} rad/s')
ax.axhline(omega_min, color='gray', lw=1, ls=':', label=f'ω_min={omega_min} rad/s')
ombres(ax, t, mode_v)
ax.set_ylabel('ω [rad/s]')
ax.set_title('Velocitat angular  (vermell=motor, verd=generador)')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# 2. Parells
ax = axes[1]
ax.plot(t, tau_g_v, color='firebrick',   lw=1,   label='τ_grav')
ax.plot(t, tau_c_v, color='forestgreen', lw=1,   label='τ_ctrl', alpha=0.8)
ax.axhline(0, color='gray', lw=0.5)
ombres(ax, t, mode_v)
ax.set_ylabel('Parell [N·m]')
ax.set_title('Parell gravitatori vs control')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# 3. Energies acumulades
Wg_c  = np.cumsum(dW_grav_v)
Wfr_c = np.cumsum(np.abs(dW_fric_v))
Wei_c = np.cumsum(dW_elec_in)
Weo_c = np.cumsum(dW_elec_out)

ax = axes[2]
ax.plot(t, Wei_c,  color='red',    lw=1.5, label='W_el_in  (consumit)')
ax.plot(t, Weo_c,  color='green',  lw=1.5, label='W_el_out (produït)')
ax.plot(t, Wfr_c,  color='orange', lw=1,   ls='--', label='W_fric')
ax.plot(t, Wg_c,   color='blue',   lw=1,   ls=':',  label='W_grav (net)')
ax.set_ylabel('Energia [J]')
ax.set_title('Energies acumulades')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# 4. Balanç net elèctric
Wnet_c = Weo_c - Wei_c
ax = axes[3]
ax.plot(t, Wnet_c, color='purple', lw=2, label='Balanç net elèctric')
ax.axhline(0, color='black', lw=1)
ax.fill_between(t, Wnet_c, 0, where=(Wnet_c>0), color='green', alpha=0.25, label='Guany')
ax.fill_between(t, Wnet_c, 0, where=(Wnet_c<0), color='red',   alpha=0.25, label='Cost')
ax.set_ylabel('Energia neta [J]')
ax.set_xlabel('Temps [s]')
ax.set_title('Balanç net elèctric acumulat (el_out − el_in)')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/kilometre_v4.png', dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v4.py  |  kilometre_v4.png")
