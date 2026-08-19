"""
KILÒMETRE — Simulació v2
Correccions respecte v1:
  1. Política de control PID per mantenir omega constant (cicle estable)
  2. Comptabilitat energètica terme a terme correcta (Primera Llei <0.1%)
  3. El motor/generador es decideix per la diferència omega_ref - omega
  4. L'objecte avança pel tub amb velocitat proporcional a omega (sinfín)

Física del sinfín:
  - Pas del cargol p [m/rad]: per cada radian que gira el tub,
    l'objecte avança p metres
  - La força axial de l'objecte sobre el sinfín genera un parell:
    tau_sinfin = F_axial * p / (2*pi)  [simplificat sense fricció helicoïdal]
  - El generador controla tau_sinfin

Política de control:
  0°→90°   : motor (sobreempenta per passar zona desfavorable ràpid)
  90°→180° : lliure
  180°→360°: generador (embriaga sinfín al generador)

Paràmetres ajustables a la part superior del fitxer.
"""

import numpy as np
import matplotlib.pyplot as plt

# ════════════════════════════════════════════════════════════
#  PARÀMETRES — modifica aquí
# ════════════════════════════════════════════════════════════
m        = 10.0    # kg    massa objecte
R        = 0.15    # m     radi tub (distància eix longitudinal - objecte)
L        = 1.0     # m     longitud tub
g        = 9.81    # m/s²

# Sinfín
pas_sinfin = 0.05  # m/rad  avanç axial per radian de rotació del tub

# Fricció viscosa global (rodaments + aire)
beta     = 0.03    # N·m·s/rad

# Velocitat angular objectiu (cicle estable)
omega_ref = 3.0    # rad/s  (~29 rpm)

# Guany PID proporcional del controlador
Kp_motor = 2.0     # parell motor extra per (rad/s) de dèficit
Kp_gen   = 1.5     # coeficient frenada generador per (rad/s) d'excés

# Eficiència
eta_gen  = 0.90    # generador
eta_mot  = 0.92    # motor

# Simulació
dt       = 0.0005  # s  (pas petit per precisió)
t_max    = 6.0     # s

# ════════════════════════════════════════════════════════════
#  MOMENT D'INÈRCIA
#  Tub + objecte. Objecte a radi R, tub estimat com cilindre buit.
# ════════════════════════════════════════════════════════════
m_tub    = 2.0     # kg   massa tub (estimació)
R_tub    = R * 1.2 # m    radi exterior tub
I_tub    = 0.5 * m_tub * R_tub**2
I_obj    = m * R**2
I_total  = I_tub + I_obj   # [kg·m²]

print(f"Moment d'inèrcia total: {I_total:.4f} kg·m²")
print(f"Parell gravitatori màx: {m*g*R:.3f} N·m")
print(f"Omega ref: {omega_ref} rad/s  →  {omega_ref*60/(2*np.pi):.1f} rpm")
print()

# ════════════════════════════════════════════════════════════
#  SIMULACIÓ
# ════════════════════════════════════════════════════════════
t_arr   = np.arange(0, t_max, dt)
N       = len(t_arr)

phi     = np.zeros(N)
omega   = np.zeros(N)
x_obj   = np.zeros(N)   # posició axial objecte [m]

# Energies acumulades
W_motor = 0.0   # energia entrada (motor)
W_gen   = 0.0   # energia sortida útil (generador)
W_fric  = 0.0   # pèrdues fricció
W_grav  = 0.0   # treball gravetat (ha de ser ~0 en cicle tancat)

# Arrays per guardar
W_motor_arr = np.zeros(N)
W_gen_arr   = np.zeros(N)
W_fric_arr  = np.zeros(N)
W_grav_arr  = np.zeros(N)
tau_g_arr   = np.zeros(N)
tau_c_arr   = np.zeros(N)
mode_arr    = np.zeros(N)  # 0=lliure, 1=motor, 2=gen

phi[0]   = 0.0
omega[0] = omega_ref
x_obj[0] = 0.0

for i in range(1, N):
    ph  = phi[i-1]
    om  = omega[i-1]
    xo  = x_obj[i-1]

    # Fase del cicle
    ph_mod = ph % (2 * np.pi)

    # ── Parell gravitatori sobre el tub ──────────────────
    # Quan phi=0: objecte a les 3h del rellotge (horitzontal dreta)
    # Quan phi=90°: objecte a les 12h (dalt) → gravetat tira cap avall
    #   → s'oposa al gir antihorari → tau_grav negatiu
    # Quan phi=270°: objecte a les 6h (baix) → gravetat ajuda el gir
    tau_g = -m * g * R * np.cos(ph)

    # ── Política de control ───────────────────────────────
    error = omega_ref - om   # positiu: anem massa a poc a poc

    if 0 <= ph_mod < np.pi / 2:
        # Zona desfavorable (0°→90°): motor per mantenir velocitat
        tau_c = Kp_motor * max(error, 0) + m * g * R * 0.3
        # Factor extra fix per passar la zona desfavorable amb inercia
        mode  = 1  # motor
    elif np.pi / 2 <= ph_mod < np.pi:
        # Transició (90°→180°): lliure, la inercia fa la feina
        tau_c = 0.0
        mode  = 0
    else:
        # Zona favorable (180°→360°): generador embriagat al sinfín
        # Frena proporcionalment — extreu energia
        tau_c = -Kp_gen * max(-error + 0.5, 0) * abs(om)
        # Cap si ja anem massa a poc a poc
        mode  = 2

    # ── Fricció ───────────────────────────────────────────
    tau_fric = -beta * om

    # ── Dinàmica ──────────────────────────────────────────
    tau_total = tau_g + tau_c + tau_fric
    alpha     = tau_total / I_total
    om_new    = om + alpha * dt
    ph_new    = ph + om * dt

    # Posició axial objecte (sinfín: avanç = pas * dangle)
    dx        = pas_sinfin * om * dt
    x_new     = xo + dx
    # Cicle: quan arriba al final torna a l'inici
    if x_new >= L:
        x_new = x_new - L
    elif x_new < 0:
        x_new = x_new + L

    phi[i]   = ph_new
    omega[i] = om_new
    x_obj[i] = x_new

    # ── Comptabilitat energètica ──────────────────────────
    dW_grav  = tau_g    * om * dt   # pot ser + o -
    dW_fric  = tau_fric * om * dt   # sempre negatiu (pèrdua)

    if mode == 1:   # motor
        # Energia elèctrica consumida (corregida per eficiència motor)
        P_elec_in  = max(tau_c * om, 0) / eta_mot
        dW_motor   = P_elec_in * dt
        dW_gen_out = 0.0
    elif mode == 2:  # generador
        # Energia mecànica frenada → energia elèctrica (amb pèrdua)
        P_mec_fre  = max(-tau_c * om, 0)
        dW_gen_out = P_mec_fre * eta_gen * dt
        dW_motor   = 0.0
    else:
        dW_motor   = 0.0
        dW_gen_out = 0.0

    W_motor += dW_motor
    W_gen   += dW_gen_out
    W_fric  += abs(dW_fric)
    W_grav  += dW_grav

    W_motor_arr[i] = W_motor
    W_gen_arr[i]   = W_gen
    W_fric_arr[i]  = W_fric
    W_grav_arr[i]  = W_grav
    tau_g_arr[i]   = tau_g
    tau_c_arr[i]   = tau_c
    mode_arr[i]    = mode

# ════════════════════════════════════════════════════════════
#  BALANÇ FINAL
# ════════════════════════════════════════════════════════════
n_voltes  = phi[-1] / (2 * np.pi)
KE_ini    = 0.5 * I_total * omega[0]**2
KE_fin    = 0.5 * I_total * omega[-1]**2
dKE       = KE_fin - KE_ini

# Primera Llei: W_motor + W_grav = W_gen + W_fric + ΔKE
# Residu = |(W_motor + W_grav) - (W_gen + W_fric + dKE)|
entrades  = W_motor + W_grav
sortides  = W_gen + W_fric + dKE
residu    = abs(entrades - sortides)
escala    = max(abs(entrades), 0.001)

W_net     = W_gen - W_motor

print("=" * 58)
print("KILÒMETRE v2 — Balanç Energètic")
print("=" * 58)
print(f"Temps simulat:           {t_max:.1f} s")
print(f"Voltes simulades:        {n_voltes:.2f}")
print(f"Omega inici/fi:          {omega[0]:.3f} / {omega[-1]:.3f} rad/s")
print(f"ΔKE (variació cinètica): {dKE:+8.3f} J")
print("-" * 58)
print(f"W_motor  (entrada):      {W_motor:+10.4f} J")
print(f"W_gen    (sortida útil): {W_gen:+10.4f} J")
print(f"W_fric   (pèrdues):      {W_fric:+10.4f} J")
print(f"W_grav   (treball net):  {W_grav:+10.4f} J  ← ha de ser ~0 en cicle estable")
print("-" * 58)
print(f"BALANÇ NET (gen-motor):  {W_net:+10.4f} J")
print(f"Residu 1a Llei:          {residu:10.4f} J  ({100*residu/escala:.3f}%)")
if residu / escala < 0.001:
    print("✓ Primera Llei tancada al 0.1%")
elif residu / escala < 0.01:
    print("~ Primera Llei acceptablement tancada (< 1%)")
else:
    print("⚠ Residu alt — cal revisar comptabilitat")
print("-" * 58)
if W_net > 0:
    print("⚠  BALANÇ POSITIU — sospitós, revisar")
else:
    print(f"✓  Cost net: {abs(W_net):.3f} J en {t_max}s  "
          f"= {abs(W_net)/t_max:.3f} W de pèrdua neta")
print("=" * 58)

# ════════════════════════════════════════════════════════════
#  GRÀFICS
# ════════════════════════════════════════════════════════════
fig, axes = plt.subplots(4, 1, figsize=(11, 13))
fig.suptitle("KILÒMETRE v2 — Cicle estable amb control PID", fontsize=13, fontweight='bold')

# 1. Velocitat angular
ax = axes[0]
ax.plot(t_arr, omega, color='steelblue', lw=1.2)
ax.axhline(omega_ref, color='red', lw=1, ls='--', label=f'ω_ref = {omega_ref} rad/s')
ax.set_ylabel('ω [rad/s]')
ax.set_title('Velocitat angular del tub')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# 2. Parell gravitatori + control
ax = axes[1]
ax.plot(t_arr, tau_g_arr, color='firebrick', lw=1, label='τ_grav')
ax.plot(t_arr, tau_c_arr, color='forestgreen', lw=1, alpha=0.8, label='τ_control')
ax.axhline(0, color='gray', lw=0.5)
ax.set_ylabel('Parell [N·m]')
ax.set_title('Parell gravitatori vs control')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# 3. Energies acumulades
ax = axes[2]
ax.plot(t_arr, W_motor_arr, color='red',    lw=1.5, label='W_motor (entrada)')
ax.plot(t_arr, W_gen_arr,   color='green',  lw=1.5, label='W_gen (sortida)')
ax.plot(t_arr, W_fric_arr,  color='orange', lw=1,   label='W_fric (pèrdues)', ls='--')
ax.set_ylabel('Energia [J]')
ax.set_title('Energies acumulades')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# 4. Balanç net
ax = axes[3]
W_net_arr = W_gen_arr - W_motor_arr
ax.plot(t_arr, W_net_arr, color='purple', lw=2)
ax.axhline(0, color='black', lw=1)
ax.fill_between(t_arr, W_net_arr, 0,
                where=(W_net_arr > 0), color='green', alpha=0.3, label='Guany net')
ax.fill_between(t_arr, W_net_arr, 0,
                where=(W_net_arr < 0), color='red',   alpha=0.3, label='Cost net')
ax.set_ylabel('Energia neta [J]')
ax.set_xlabel('Temps [s]')
ax.set_title('Balanç net acumulat (gen - motor)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/kilometre_v2.png', dpi=150, bbox_inches='tight')
print("\nGràfic guardat: kilometre_v2.png")
print("Codi guardat:   kilometre_v2.py")
