"""
KILÒMETRE — Simulació v6
Sistema clarificat:
  - Motor extern manté omega constant (omega_ref)
  - El motor aporta exactament el que cal per compensar
    tau_grav + tau_gen + tau_fric en cada instant
  - Generador embriagat al sinfín NOMÉS quan tau_grav > 0
    (fase favorable: 180°→360°, objecte a la part baixa)
  - Comptabilitat neta: W_gen_elec vs W_motor_elec

Amb omega constant, la física és neta:
  tau_motor = -(tau_grav + tau_gen + tau_fric)
  en cada instant (no hi ha acceleració → ΔKE = 0)

Primera Llei verificació:
  W_grav_net + W_motor_mec + W_gen_mec - W_fric = ΔKE = 0
  → W_grav_net ha de ser ~0 en cicle complet
  → W_motor_mec + W_gen_mec = W_fric (tot va a fricció si és neutre)

Clau: el generador extreu energia MECÀNICA del sistema via sinfín.
Aquesta energia ve de: la gravetat (fase favorable) + el motor (si cal).
Si en fase favorable tau_grav > |tau_gen|, el motor fins i tot pot
rebre parell (funcionar com a generador secundari).
"""

import numpy as np
import matplotlib.pyplot as plt

# ════════════════════════════════════════════════════════════
#  PARÀMETRES
# ════════════════════════════════════════════════════════════
m        = 10.0    # kg
R        = 0.15    # m
g        = 9.81    # m/s²
pas      = 0.04    # m/rad

m_tub    = 2.0;  R_tub = 0.18
I        = 0.5*m_tub*R_tub**2 + m*R**2

beta     = 0.02    # N·m·s/rad  fricció viscosa
tau_gmax = m * g * R   # 14.715 N·m

omega_ref = 3.0    # rad/s  (constant, mantingut pel motor)

# Generador: fracció del parell gravitatori que s'embriaga
# quan tau_grav > 0 (fase 180°→360°)
alfa_gen  = 0.60   # 0=res, 1=tot el parell gravitatori

eta_gen   = 0.90
eta_mot   = 0.92

# Simulació — omega és CONSTANT, una sola volta per analitzar
dt        = 0.0001   # s
# Una volta completa = 2π / omega_ref
T_volta   = 2 * np.pi / omega_ref
N_volta   = int(T_volta / dt)
N_cicles  = 20       # simula 20 voltes per veure règim estacionari
N         = N_volta * N_cicles

t_arr = np.linspace(0, N * dt, N)

# ════════════════════════════════════════════════════════════
#  BUCLE (omega = omega_ref constant)
# ════════════════════════════════════════════════════════════
om = omega_ref

dW_grav_v   = np.zeros(N)
dW_fric_v   = np.zeros(N)
dW_gen_mec  = np.zeros(N)   # energia mecànica extreta pel gen (>0)
dW_mot_mec  = np.zeros(N)   # energia mecànica aportada pel motor (>0)
dW_elec_in  = np.zeros(N)   # elèctric consumit pel motor
dW_elec_out = np.zeros(N)   # elèctric produït pel gen
tau_g_v     = np.zeros(N)
tau_gen_v   = np.zeros(N)
tau_mot_v   = np.zeros(N)
mode_v      = np.zeros(N, dtype=int)  # 0=lliure, 2=gen actiu

for i in range(N):
    ph     = omega_ref * i * dt
    ph_mod = ph % (2 * np.pi)

    # Parell gravitatori
    tau_g  = -m * g * R * np.cos(ph)

    # Fricció
    tau_fr = -beta * om

    # Generador (sinfín embriagat) NOMÉS en fase favorable
    # Fase favorable: tau_g > 0 → cos(ph) < 0 → ph_mod ∈ (π/2, 3π/2)
    # Però per la idea del sistema: 180°→360°
    if ph_mod > np.pi and om > 0:
        # tau_g pot ser positiu o negatiu en 180°→360°:
        #   180°→270°: cos<0 → tau_g>0 (ajuda)
        #   270°→360°: cos>0 → tau_g<0 (frena)
        # Embriguem SEMPRE en 180°→360° per capturar la inèrcia
        tau_gen = -alfa_gen * abs(tau_g + beta*om)
        mode    = 2
    else:
        tau_gen = 0.0
        mode    = 0

    # Motor: aporta exactament el que cal per mantindre omega constant
    # tau_motor + tau_g + tau_gen + tau_fr = 0  (ΔKE = 0)
    tau_mot = -(tau_g + tau_gen + tau_fr)

    # Comptabilitat energètica (tot en mecànica primer)
    dWg   = tau_g   * om * dt   # pot ser + o -
    dWfr  = tau_fr  * om * dt   # ≤ 0
    dWgen = tau_gen * om * dt   # ≤ 0 (extreu mecànica)
    dWmot = tau_mot * om * dt   # pot ser + o - !

    dW_grav_v[i] = dWg
    dW_fric_v[i] = dWfr

    # Generador: extreu mecànica → produeix elèctrica
    if dWgen < 0:
        dW_gen_mec[i]  = abs(dWgen)
        dW_elec_out[i] = abs(dWgen) * eta_gen

    # Motor: pot estar aportant O recuperant energia
    if dWmot > 0:
        # Motor consumix elèctrica per aportar mecànica
        dW_mot_mec[i]  = dWmot
        dW_elec_in[i]  = dWmot / eta_mot
    elif dWmot < 0:
        # El tub "empeny" el motor → motor actua com a gen secundari
        dW_mot_mec[i]  = 0.0
        dW_elec_in[i]  = 0.0
        # Energia recuperada pel motor en mode gen
        dW_elec_out[i] += abs(dWmot) * eta_gen  # suma al gen principal

    tau_g_v[i]   = tau_g
    tau_gen_v[i] = tau_gen
    tau_mot_v[i] = tau_mot
    mode_v[i]    = mode

# ════════════════════════════════════════════════════════════
#  BALANÇ
# ════════════════════════════════════════════════════════════
W_grav    = np.sum(dW_grav_v)
W_fric    = np.sum(np.abs(dW_fric_v))
W_gen_mec = np.sum(dW_gen_mec)
W_mot_mec = np.sum(dW_mot_mec)
W_el_in   = np.sum(dW_elec_in)
W_el_out  = np.sum(dW_elec_out)
W_net_el  = W_el_out - W_el_in
n_voltes  = (omega_ref * N * dt) / (2*np.pi)

# Verificació 1a Llei (omega constant → ΔKE = 0)
# W_grav + W_motor_mec + W_gen_mec_net - W_fric = 0
check_1llei = W_grav + W_mot_mec + (-W_gen_mec) - W_fric
escala      = max(W_fric, W_gen_mec, W_mot_mec, 0.001)
residu_rel  = abs(check_1llei) / escala * 100

print("=" * 66)
print("KILÒMETRE v6 — Motor manté omega constant, gen via sinfín")
print("=" * 66)
print(f"omega_ref={omega_ref} rad/s ({omega_ref*60/(2*np.pi):.1f} rpm)  "
      f"|  alfa_gen={alfa_gen}  |  beta={beta}")
print(f"Voltes simulades: {n_voltes:.1f}  |  T_volta={T_volta:.3f}s")
print("-" * 66)
print(f"W_grav  (net {n_voltes:.0f} voltes): {W_grav:+10.4f} J  ← ha de ser ~0")
print(f"W_mot   (mec. aportada):   {W_mot_mec:+10.4f} J")
print(f"W_gen   (mec. extreta):   {-W_gen_mec:+10.4f} J")
print(f"W_fric  (pèrdues):        {-W_fric:+10.4f} J")
print(f"  Check 1a Llei (ha ser 0): {check_1llei:+.6f} J  ({residu_rel:.4f}%)")
print("-" * 66)
print(f"W_el_in  (motor, xarxa):   {W_el_in:+10.4f} J")
print(f"W_el_out (gen elèctric):   {W_el_out:+10.4f} J")
print(f"BALANÇ NET ELÈCTRIC:       {W_net_el:+10.4f} J")
print(f"Potència neta:             {W_net_el/(N*dt)*1000:+10.2f} mW")
print("-" * 66)

# Anàlisi per volta
print("\nPer volta (règim estacionari = totes iguals si cicle tancat):")
print(f"{'V':>3} | {'W_grav':>8} | {'W_mot':>8} | {'W_gen':>8} | "
      f"{'W_fric':>8} | {'el_in':>8} | {'el_out':>8} | {'net_el':>8}")
print("-" * 80)
for c in range(N_cicles):
    s   = slice(c*N_volta, (c+1)*N_volta)
    wg  = np.sum(dW_grav_v[s])
    wm  = np.sum(dW_mot_mec[s])
    wgn = np.sum(dW_gen_mec[s])
    wf  = np.sum(np.abs(dW_fric_v[s]))
    ei  = np.sum(dW_elec_in[s])
    eo  = np.sum(dW_elec_out[s])
    print(f"{c+1:>3} | {wg:>8.5f} | {wm:>8.5f} | {-wgn:>8.5f} | "
          f"{-wf:>8.5f} | {ei:>8.5f} | {eo:>8.5f} | {eo-ei:>8.5f}")

# Resum per una volta (promig règim estacionari)
s_re = slice(5*N_volta, N)   # descarta les primeres 5 voltes (transitori)
wg_re  = np.sum(dW_grav_v[s_re])   / (N_cicles-5)
wm_re  = np.sum(dW_mot_mec[s_re])  / (N_cicles-5)
wgn_re = np.sum(dW_gen_mec[s_re])  / (N_cicles-5)
wf_re  = np.sum(np.abs(dW_fric_v[s_re])) / (N_cicles-5)
ei_re  = np.sum(dW_elec_in[s_re])  / (N_cicles-5)
eo_re  = np.sum(dW_elec_out[s_re]) / (N_cicles-5)
net_re = eo_re - ei_re

print("-" * 80)
print(f"{'PROMIG':>3} | {wg_re:>8.5f} | {wm_re:>8.5f} | {-wgn_re:>8.5f} | "
      f"{-wf_re:>8.5f} | {ei_re:>8.5f} | {eo_re:>8.5f} | {net_re:>8.5f}")
print()
if residu_rel < 1.0:
    print(f"✓ Primera Llei tancada al {residu_rel:.4f}%")
else:
    print(f"⚠ Residu 1a Llei: {residu_rel:.3f}% — revisar")

if net_re > 0.001:
    print(f"\n⚠  Net elèctric positiu per volta: +{net_re:.5f} J")
    print(f"   Però: W_grav net = {wg_re:+.5f} J (ha de ser 0 per ser vàlid)")
    print(f"   Si W_grav≠0, estem consumint energia potencial, no hurtant-la.")
else:
    print(f"\n✓  Net elèctric negatiu per volta: {net_re:.5f} J")
    print(f"   El motor paga: fricció + pèrdues conversió")
    print(f"   El gen recupera una part, però no tota.")

# ════════════════════════════════════════════════════════════
#  ANÀLISI D'UNA VOLTA — tau per angle
# ════════════════════════════════════════════════════════════
print("\n--- Anàlisi d'una volta (angles 0°→360°) ---")
angles  = np.linspace(0, 360, 361)
ph_rad  = np.radians(angles)
tau_g_a = -m * g * R * np.cos(ph_rad)
tau_fr_a= -beta * omega_ref * np.ones_like(ph_rad)
tau_gen_a = np.where(ph_rad % (2*np.pi) > np.pi,
                     -alfa_gen * np.abs(tau_g_a + beta*omega_ref), 0.0)
tau_mot_a = -(tau_g_a + tau_gen_a + tau_fr_a)

print(f"Zona motor     (tau_mot>0): "
      f"{np.sum(tau_mot_a>0)} graus  →  "
      f"W_mot = {np.trapezoid(np.maximum(tau_mot_a,0), ph_rad)*omega_ref*dt*N_volta/(2*np.pi):.4f} J/volta est.")
print(f"Zona gen       (tau_gen<0): "
      f"{np.sum(tau_gen_a<0)} graus  →  "
      f"W_gen = {np.trapezoid(np.abs(np.minimum(tau_gen_a,0)), ph_rad)*omega_ref*dt*N_volta/(2*np.pi):.4f} J/volta est.")
print(f"Integral tau_g sobre 2π:   {np.trapezoid(tau_g_a, ph_rad):.6f} N·m·rad  (ha ser 0)")

# ════════════════════════════════════════════════════════════
#  GRÀFICS
# ════════════════════════════════════════════════════════════
t  = t_arr
ph = omega_ref * np.arange(N) * dt

fig, axes = plt.subplots(4, 1, figsize=(13, 15))
fig.suptitle(f"KILÒMETRE v6 — omega constant={omega_ref} rad/s  |  "
             f"alfa_gen={alfa_gen}  |  beta={beta}\n"
             f"Balanç net elèctric per volta: {net_re:+.5f} J  "
             f"({'GUANY' if net_re>0 else 'COST'})",
             fontsize=11, fontweight='bold')

# 1. Parells vs angle (una volta)
ax = axes[0]
ang = np.degrees(ph[:N_volta] % (2*np.pi))
ax.plot(ang, tau_g_v[:N_volta],   'firebrick',   lw=1.5, label='τ_grav')
ax.plot(ang, tau_gen_v[:N_volta], 'forestgreen', lw=1.5, label='τ_gen (sinfín)')
ax.plot(ang, tau_mot_v[:N_volta], 'steelblue',   lw=1.5, label='τ_motor')
ax.axhline(0, color='gray', lw=0.5)
ax.axvline(180, color='purple', lw=1, ls='--', label='180° (inici gen)')
ax.fill_between(ang, tau_gen_v[:N_volta], 0,
                where=(tau_gen_v[:N_volta]<0), color='green', alpha=0.2)
ax.fill_between(ang, tau_mot_v[:N_volta], 0,
                where=(tau_mot_v[:N_volta]>0), color='red', alpha=0.15)
ax.fill_between(ang, tau_mot_v[:N_volta], 0,
                where=(tau_mot_v[:N_volta]<0), color='blue', alpha=0.15,
                label='Motor recupera (gen sec.)')
ax.set_xlabel('Angle [°]');  ax.set_ylabel('Parell [N·m]')
ax.set_title('Parells vs angle — 1 volta')
ax.set_xlim(0, 360);  ax.legend(fontsize=9);  ax.grid(True, alpha=0.3)

# 2. Potència instantània vs angle
ax = axes[1]
P_grav = tau_g_v[:N_volta]   * om
P_gen  = tau_gen_v[:N_volta] * om
P_mot  = tau_mot_v[:N_volta] * om
ax.plot(ang, P_grav, 'firebrick',   lw=1.5, label='P_grav')
ax.plot(ang, P_gen,  'forestgreen', lw=1.5, label='P_gen (extreta)')
ax.plot(ang, P_mot,  'steelblue',   lw=1.5, label='P_motor')
ax.axhline(0, color='gray', lw=0.5)
ax.axvline(180, color='purple', lw=1, ls='--')
ax.fill_between(ang, P_gen, 0, where=(P_gen<0), color='green', alpha=0.2)
ax.set_xlabel('Angle [°]');  ax.set_ylabel('Potència [W]')
ax.set_title('Potència instantània vs angle')
ax.set_xlim(0, 360);  ax.legend(fontsize=9);  ax.grid(True, alpha=0.3)

# 3. Energies acumulades (tot el temps)
Wei_c = np.cumsum(dW_elec_in)
Weo_c = np.cumsum(dW_elec_out)
Wfr_c = np.cumsum(np.abs(dW_fric_v))
Wg_c  = np.cumsum(dW_grav_v)
ax = axes[2]
ax.plot(t, Wei_c, 'red',    lw=1.5, label='W_el_in  (motor consumit)')
ax.plot(t, Weo_c, 'green',  lw=1.5, label='W_el_out (gen produït)')
ax.plot(t, Wfr_c, 'orange', lw=1,   ls='--', label='W_fric acum.')
ax.plot(t, Wg_c,  'blue',   lw=1,   ls=':',  label='W_grav net acum.')
ax.set_ylabel('Energia [J]');  ax.set_xlabel('Temps [s]')
ax.set_title(f'Energies acumulades ({N_cicles} voltes)')
ax.legend(fontsize=9);  ax.grid(True, alpha=0.3)

# 4. Balanç net per volta
nets = []
for c in range(N_cicles):
    s  = slice(c*N_volta, (c+1)*N_volta)
    ei = np.sum(dW_elec_in[s])
    eo = np.sum(dW_elec_out[s])
    nets.append(eo - ei)
ax = axes[3]
ax.bar(range(1, N_cicles+1), nets,
       color=['green' if n>0 else 'red' for n in nets], alpha=0.7)
ax.axhline(0, color='black', lw=1)
ax.axhline(np.mean(nets[5:]), color='purple', lw=2, ls='--',
           label=f'Promig règim estac.: {np.mean(nets[5:]):.5f} J/volta')
ax.set_xlabel('Volta');  ax.set_ylabel('Net elèctric [J]')
ax.set_title('Balanç net elèctric per volta')
ax.legend(fontsize=9);  ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/kilometre_v6.png', dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v6.py  |  kilometre_v6.png")
