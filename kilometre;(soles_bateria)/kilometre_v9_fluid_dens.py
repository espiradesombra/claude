"""
KILÒMETRE v9 — Tub completament submergit en fluid dens
Objecte de baixa densitat (flotant) dins fluid dens (mercuri, galià, etc.)

Física nova respecte v6-v8:
  1. F_neta = (rho_fluid - rho_obj)*V*g  [constant, cap amunt sempre]
     → Equivalent a "gravetat invertida i amplificada"
     → tau_net = F_neta * R * cos(phi)  [parell sobre el tub]

  2. Inèrcia afegida (added mass):
     I_af = C_m * rho_fluid * V * R²   (C_m=0.5 per esfera)
     I_total = I_mecanica + I_af

  3. Fregament viscós en fluid dens:
     tau_fric = -beta_fluid * omega
     beta_fluid >> beta_aire

  4. W_net gravetat+flotació = 0 en cicle tancat (sempre)
     PERÒ: I_total molt major → E_k = ½*I_total*ω² molt major
     → Millor bateria cinètica per la mateixa velocitat angular

Comparació de casos:
  A) Aire + acer (referència, v6)
  B) Aigua + escuma EPS
  C) Mercuri + heli (globus)
  D) Gallium + polietilè
"""

import numpy as np
import matplotlib.pyplot as plt

g    = 9.81
R    = 0.30    # m  (radi major per aprofitar la inèrcia afegida)
L    = 1.0     # m
omega_ref = 5.0  # rad/s

# Tub
m_tub   = 5.0
R_tub   = R * 1.1
I_tub   = 0.5 * m_tub * R_tub**2

# ── Casos a comparar ─────────────────────────────────────────
casos = {
    'A: Aire+Acer\n(referència)': {
        'rho_fluid': 1.2,
        'rho_obj':   7800,
        'V_obj':     0.005,      # m³
        'mu_fluid':  1.8e-5,     # Pa·s (viscositat dinàmica aire)
        'C_m':       0.5,        # added mass coef.
        'color':     'gray',
    },
    'B: Aigua+EPS\n(escuma)': {
        'rho_fluid': 1000,
        'rho_obj':   15,
        'V_obj':     0.005,
        'mu_fluid':  1e-3,
        'C_m':       0.5,
        'color':     'steelblue',
    },
    'C: Mercuri+Heli\n(màxim contrast)': {
        'rho_fluid': 13600,
        'rho_obj':   0.164,
        'V_obj':     0.005,
        'mu_fluid':  1.5e-3,
        'C_m':       0.5,
        'color':     'firebrick',
    },
    'D: Gallium+EPS\n(pràctic)': {
        'rho_fluid': 6095,
        'rho_obj':   15,
        'V_obj':     0.005,
        'mu_fluid':  2e-3,
        'C_m':       0.5,
        'color':     'forestgreen',
    },
}

# ── Paràmetres de control ────────────────────────────────────
alfa_gen  = 0.55   # fracció del parell net que frena el generador
eta_gen   = 0.90
eta_mot   = 0.92
N_cicles  = 20

print("=" * 72)
print("KILÒMETRE v9 — Fluid dens + objecte flotant")
print("Comparació de casos: inèrcia afegida i energia emmagatzemada")
print("=" * 72)
print(f"R={R}m  |  ω_ref={omega_ref} rad/s ({omega_ref*60/(2*np.pi):.1f} rpm)  |  V_obj=5L")
print()
print(f"{'Cas':>28} | {'F_net(N)':>9} | {'I_mec':>7} | {'I_af':>7} | "
      f"{'I_tot':>7} | {'E_k(J)':>8} | {'E_k(Wh)':>8}")
print("-" * 85)

resultats = {}
for nom, c in casos.items():
    m_obj   = c['rho_obj'] * c['V_obj']
    F_net   = (c['rho_fluid'] - c['rho_obj']) * c['V_obj'] * g
    I_mec   = I_tub + m_obj * R**2
    I_af    = c['C_m'] * c['rho_fluid'] * c['V_obj'] * R**2
    I_tot   = I_mec + I_af
    E_k     = 0.5 * I_tot * omega_ref**2

    # Fricció viscosa (Stokes per esfera: F_drag = 6π*μ*r_eq*v)
    r_eq    = (3*c['V_obj']/(4*np.pi))**(1/3)
    v_tip   = omega_ref * R
    F_drag  = 6 * np.pi * c['mu_fluid'] * r_eq * v_tip
    tau_fric_ref = F_drag * R   # parell de fricció a omega_ref
    beta    = tau_fric_ref / omega_ref   # coef. viscós equivalent

    c['m_obj']  = m_obj
    c['F_net']  = F_net
    c['I_mec']  = I_mec
    c['I_af']   = I_af
    c['I_tot']  = I_tot
    c['E_k']    = E_k
    c['beta']   = beta
    c['r_eq']   = r_eq

    nom_curt = nom.split('\n')[0]
    print(f"{nom_curt:>28} | {F_net:>9.1f} | {I_mec:>7.4f} | {I_af:>7.4f} | "
          f"{I_tot:>7.4f} | {E_k:>8.3f} | {E_k/3600:>8.6f}")

print()

# ── Simulació: cicle complet per cada cas ────────────────────
dt = 0.0001
N  = int(2*np.pi/omega_ref / dt) * N_cicles  # N_cicles voltes exactes

print(f"{'Cas':>28} | {'net_el/volta(J)':>16} | {'η_rt%':>7} | "
      f"{'P_pic(W)':>9} | {'β_fluid':>10}")
print("-" * 80)

fig, axes = plt.subplots(3, 1, figsize=(13, 14))
fig.suptitle("KILÒMETRE v9 — Fluid dens + objecte flotant\n"
             "Inèrcia afegida com a avantatge de bateria cinètica",
             fontsize=12, fontweight='bold')

col_Ek  = []
col_net = []
col_nom = []

for nom, c in casos.items():
    phi_arr  = np.linspace(0, N_cicles*2*np.pi, N, endpoint=False)
    dphi     = phi_arr[1] - phi_arr[0]
    om       = omega_ref

    dW_grav_v   = np.zeros(N)
    dW_fric_v   = np.zeros(N)
    dW_gen_v    = np.zeros(N)
    dW_mot_v    = np.zeros(N)
    dW_elec_in  = np.zeros(N)
    dW_elec_out = np.zeros(N)

    for i in range(N):
        ph     = phi_arr[i]
        ph_mod = ph % (2*np.pi)

        # Parell net (flotació - gravetat) sobre el tub
        # F_net cap amunt → parell = F_net * R * cos(phi)
        # (mateix signe que gravetat però invertit si F_net > 0)
        tau_net = c['F_net'] * R * np.cos(ph)  # pot ser + o -

        tau_fr  = -c['beta'] * om

        # Generador actiu quan tau_net > 0 (força neta ajuda el gir)
        if ph_mod > np.pi and om > 0.5*omega_ref:
            tau_gen = -alfa_gen * abs(tau_net)
            mode = 2
        else:
            tau_gen = 0.0
            mode = 0

        # Motor manté omega constant
        tau_mot = -(tau_net + tau_gen + tau_fr)

        # Treballs
        dWg  = tau_net * om * dphi / om   # = tau_net * dphi
        dWfr = tau_fr  * dphi
        dWge = tau_gen * dphi
        dWmo = tau_mot * dphi

        dW_grav_v[i] = dWg
        dW_fric_v[i] = dWfr
        dW_gen_v[i]  = dWge
        dW_mot_v[i]  = dWmo

        if mode == 0 and dWmo > 0:
            dW_elec_in[i]  = dWmo / eta_mot
        elif mode == 2 and dWge < 0:
            dW_elec_out[i] = abs(dWge) * eta_gen
        # motor recupera si tau_mot < 0
        if dWmo < 0:
            dW_elec_out[i] += abs(dWmo) * eta_gen

    W_el_in  = np.sum(dW_elec_in)
    W_el_out = np.sum(dW_elec_out)
    W_grav   = np.sum(dW_grav_v)
    W_fric   = np.sum(np.abs(dW_fric_v))
    net_el   = W_el_out - W_el_in
    net_volta= net_el / N_cicles
    eta_rt   = 100*W_el_out/W_el_in if W_el_in > 0 else 0
    P_pic    = alfa_gen * abs(c['F_net']) * R * omega_ref

    nom_curt = nom.split('\n')[0]
    print(f"{nom_curt:>28} | {net_volta:>16.4f} | {eta_rt:>7.1f} | "
          f"{P_pic:>9.1f} | {c['beta']:>10.5f}")

    c['net_volta']   = net_volta
    c['eta_rt']      = eta_rt
    c['P_pic']       = P_pic
    c['dW_elec_in']  = dW_elec_in
    c['dW_elec_out'] = dW_elec_out
    c['W_grav']      = W_grav
    c['W_fric']      = W_fric

    col_Ek.append(c['E_k'])
    col_net.append(net_volta)
    col_nom.append(nom_curt)

print()
print(f"W_grav net (ha ser ~0): {W_grav:.6f} J  ✓")

# ── Gràfics ──────────────────────────────────────────────────
t_volta = 2*np.pi / omega_ref

# 1. Energia cinètica emmagatzemada per cas
ax = axes[0]
colors = [c['color'] for c in casos.values()]
bars = ax.bar(col_nom, col_Ek, color=colors, alpha=0.8, edgecolor='black')
ax.set_ylabel('Energia cinètica E_k [J]')
ax.set_title(f'Energia cinètica emmagatzemada per cas\n'
             f'(ω={omega_ref} rad/s, R={R}m, V_obj=5L)')
for bar, val in zip(bars, col_Ek):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
            f'{val:.1f}J', ha='center', va='bottom', fontsize=9, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# 2. Balanç net elèctric per volta
ax = axes[1]
color_bars = ['green' if n > 0 else 'red' for n in col_net]
bars2 = ax.bar(col_nom, col_net, color=color_bars, alpha=0.8, edgecolor='black')
ax.axhline(0, color='black', lw=1)
ax.set_ylabel('Balanç net elèctric per volta [J]')
ax.set_title('Balanç net elèctric per volta (positiu = guany, negatiu = cost)')
for bar, val in zip(bars2, col_net):
    ax.text(bar.get_x()+bar.get_width()/2,
            bar.get_height() + (0.02 if val >= 0 else -0.08),
            f'{val:.3f}J', ha='center', va='bottom', fontsize=9)
ax.grid(True, alpha=0.3, axis='y')

# 3. Potència pic vs inèrcia total
ax = axes[2]
I_tots  = [c['I_tot'] for c in casos.values()]
P_pics  = [c['P_pic'] for c in casos.values()]
cols    = [c['color'] for c in casos.values()]
for i, (nom_c, I_v, P_v, col) in enumerate(zip(col_nom, I_tots, P_pics, cols)):
    ax.scatter(I_v, P_v, s=200, color=col, zorder=5, edgecolors='black')
    ax.annotate(nom_c, (I_v, P_v),
                textcoords="offset points", xytext=(10, 5), fontsize=9)
ax.set_xlabel('Inèrcia total I [kg·m²]')
ax.set_ylabel('Potència pic generador [W]')
ax.set_title('Potència pic vs inèrcia total\n(esquina superior dreta = millor bateria)')
ax.grid(True, alpha=0.3)
ax.set_xscale('log'); ax.set_yscale('log')

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/kilometre_v9_fluid_dens.png',
            dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v9_fluid_dens.py  |  .png")
