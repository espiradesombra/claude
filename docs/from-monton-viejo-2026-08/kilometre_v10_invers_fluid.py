"""
KILÒMETRE v10 — Estudi invers: fluid dens + objecte flotant
Comptabilitat corregida + estudi invers de dimensions

Física neta (omega constant, motor ideal):
  tau_net(phi) = (rho_f - rho_o)*V*g * R * cos(phi)   [parell flotació-gravetat]
  tau_fric     = -beta * omega
  tau_mot      = -(tau_net + tau_gen + tau_fric)        [motor compensa tot]

Energia per volta (analítica, sense simulació):
  W_net_mec = ∮ tau_net dphi = 0   (sempre, conservatiu)
  W_fric_volta = beta * omega * 2*pi
  W_gen_mec_volta = alfa_gen * |F_net| * R * ∮_{180→360} |cos(phi)| dphi
                  = alfa_gen * |F_net| * R * 2
  W_mot_volta = W_fric + W_gen_mec  (motor paga fricció + el que extreu el gen)

  W_el_out = W_gen_mec * eta_gen  +  W_mot_regen * eta_gen
  W_el_in  = W_mot_actiu / eta_mot

  NET = W_el_out - W_el_in  → ha de ser NEGATIU (motor paga les pèrdues)

Inèrcia efectiva (bateria):
  I_ef = I_mec + I_afegida = m_obj*R² + C_m*rho_f*V*R²
       = R² * (m_obj + C_m*rho_f*V)

Estudi invers:
  Donats: E_k objectiu, fluid, objecte
  → Quins R i omega cal?
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

g      = 9.81
C_m    = 0.5     # added mass coef. esfera
eta_gen = 0.90
eta_mot = 0.92
alfa_gen = 0.55  # fracció parell extret pel gen

# ── Fluids i objectes ────────────────────────────────────────
fluids = {
    'Agua':        {'rho': 1000,  'mu': 1.0e-3,  'color': 'steelblue',  'cost': 'molt baix'},
    'Salmuera 25%':{'rho': 1200,  'mu': 1.5e-3,  'color': 'cyan',       'cost': 'molt baix'},
    'Glicerol':    {'rho': 1260,  'mu': 1.4,     'color': 'orange',     'cost': 'baix'},
    'Gallium':     {'rho': 6095,  'mu': 2.0e-3,  'color': 'purple',     'cost': 'molt alt'},
    'Mercuri':     {'rho': 13600, 'mu': 1.5e-3,  'color': 'red',        'cost': 'alt+tòxic'},
}

objectes = {
    'Escuma EPS':  {'rho': 15,    'color': 'lightblue'},
    'Fusta balsa': {'rho': 120,   'color': 'peru'},
    'Heli (globus)':{'rho': 0.164,'color': 'yellow'},
    'Polietilè':   {'rho': 950,   'color': 'gray'},
}

# ── Funcions analítiques ─────────────────────────────────────
def I_efectiva(rho_f, rho_o, V, R):
    m_obj = rho_o * V
    I_mec = m_obj * R**2
    I_af  = C_m * rho_f * V * R**2
    return I_mec + I_af, I_mec, I_af

def F_neta(rho_f, rho_o, V):
    return (rho_f - rho_o) * V * g

def beta_viscosa(mu_f, V, omega, R):
    """Fricció viscosa via Stokes (esfera)"""
    r_eq  = (3*V/(4*np.pi))**(1/3)
    v_tip = omega * R
    F_d   = 6 * np.pi * mu_f * r_eq * v_tip
    return F_d * R / omega   # beta = tau_fric / omega

def balanc_volta(rho_f, rho_o, V, R, omega, mu_f):
    """
    Balanç energètic per volta — analític net.
    Retorna: W_el_in, W_el_out, net_el [J/volta]
    """
    Fn   = F_neta(rho_f, rho_o, V)
    beta = beta_viscosa(mu_f, V, omega, R)

    # Treball gen (zona 180-360, ∫|cos|dphi = 2)
    W_gen_mec  = alfa_gen * abs(Fn) * R * 2.0    # J/volta mecànic
    W_fric     = beta * omega * 2*np.pi            # J/volta fricció

    # Motor: paga fricció + gen - el que recupera quan tau_net ajuda
    # En zona 0-180: tau_net pot ser positiu (ajuda) → motor recupera
    # ∫_0^pi cos(phi)dphi = 0, però el motor no actua uniformement.
    # Simplificació: motor paga W_fric + W_gen_mec net
    # (comptabilitat conservadora, no compta regeneració motor)
    W_mot_mec  = W_fric + W_gen_mec

    W_el_in    = W_mot_mec / eta_mot
    W_el_out   = W_gen_mec * eta_gen
    net_el     = W_el_out - W_el_in
    return W_el_in, W_el_out, net_el, W_fric, W_gen_mec

# ── Verificació: el net ha de ser negatiu ────────────────────
print("=" * 70)
print("VERIFICACIÓ COMPTABILITAT (analítica)")
print("Net elèctric per volta — ha de ser NEGATIU (motor paga pèrdues)")
print("=" * 70)
R_ref = 0.30; omega_ref = 5.0; V_ref = 0.005

for f_nom, fl in fluids.items():
    for o_nom, ob in objectes.items():
        rho_f = fl['rho']; rho_o = ob['rho']
        if rho_f <= rho_o:
            continue   # objecte no flota
        Wi, Wo, net, Wf, Wg = balanc_volta(rho_f, rho_o, V_ref,
                                            R_ref, omega_ref, fl['mu'])
        if abs(net) > 0.1:
            print(f"  {f_nom:>16}+{o_nom:<14}: net={net:+8.3f} J/v  "
                  f"(gen={Wg:.2f} fric={Wf:.4f})")

# ── Estudi invers ────────────────────────────────────────────
print()
print("=" * 75)
print("ESTUDI INVERS: de l'objectiu E_k als paràmetres R i ω")
print("E_k = ½ * I_ef * ω²  |  I_ef = R²*(m_obj + C_m*rho_f*V)")
print("=" * 75)

objectius_Ek = {
    '100 J  (~0.03 Wh)': 100,
    '1 kJ   (~0.28 Wh)': 1e3,
    '10 kJ  (~2.8 Wh)':  10e3,
    '100 kJ (~28 Wh)':   100e3,
    '1 MJ   (~278 Wh)':  1e6,
}

V_obj = 0.010   # 10 litres (objecte)
Rs    = np.array([0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0])
omeges_max = {   # límit pràctic per fluid (cavitació, contenció)
    'Agua':        300,   # rad/s
    'Salmuera 25%':300,
    'Glicerol':    100,   # viscós, lent
    'Gallium':     200,
    'Mercuri':     150,   # contenció + seguretat
}

# Cas principal: Mercuri + Heli i Gallium + EPS
casos_inv = [
    ('Mercuri', 'Heli (globus)', fluids['Mercuri'], objectes['Heli (globus)']),
    ('Gallium', 'Escuma EPS',   fluids['Gallium'],  objectes['Escuma EPS']),
    ('Salmuera 25%', 'Escuma EPS', fluids['Salmuera 25%'], objectes['Escuma EPS']),
    ('Agua', 'Escuma EPS',     fluids['Agua'],      objectes['Escuma EPS']),
]

for E_nom, E_obj in objectius_Ek.items():
    print(f"\n{'─'*75}")
    print(f"OBJECTIU: {E_nom}  ({E_obj:.2e} J)")
    print(f"{'Fluid+Obj':>22} | {'R(m)':>5} | {'ω(r/s)':>7} | {'rpm':>6} | "
          f"{'ω_max?':>7} | {'P_pic(W)':>9} | {'net/v(J)':>9}")
    print("-" * 75)

    for f_nom, o_nom, fl, ob in casos_inv:
        rho_f = fl['rho']; rho_o = ob['rho']
        if rho_f <= rho_o:
            continue
        om_max_fl = omeges_max.get(f_nom, 200)

        for R in [0.3, 0.5, 1.0, 1.5]:
            I_ef, I_m, I_a = I_efectiva(rho_f, rho_o, V_obj, R)
            # ω necessària per assolir E_obj
            om_need = np.sqrt(2 * E_obj / I_ef)
            rpm     = om_need * 60 / (2*np.pi)
            viable  = om_need <= om_max_fl
            P_pic   = alfa_gen * F_neta(rho_f,rho_o,V_obj) * R * om_need
            Wi, Wo, net, Wf, Wg = balanc_volta(rho_f, rho_o, V_obj,
                                                R, om_need, fl['mu'])
            mark = "✓" if viable else f"✗(max{om_max_fl:.0f})"
            nom_cas = f"{f_nom}+{o_nom[:3]}"
            print(f"{nom_cas:>22} | {R:>5.1f} | {om_need:>7.1f} | {rpm:>6.0f} | "
                  f"{mark:>7} | {P_pic:>9.1f} | {net:>9.3f}")

# ── Taula resum: densitat energètica ─────────────────────────
print()
print("=" * 70)
print("DENSITAT ENERGÈTICA MÀXIMA (Wh/kg fluid+objecte)")
print(f"V_obj={V_obj*1000:.0f}L  |  R=0.5m")
print("=" * 70)
R_d = 0.5
print(f"{'Fluid':>16} | {'Objecte':>14} | {'I_ef(kg·m²)':>12} | "
      f"{'ω_max':>7} | {'E_max(Wh)':>10} | {'m_fluid(kg)':>12} | {'Wh/kg':>8}")
print("-" * 90)

for f_nom, fl in fluids.items():
    for o_nom, ob in objectes.items():
        rho_f = fl['rho']; rho_o = ob['rho']
        if rho_f <= rho_o: continue
        I_ef, _, _ = I_efectiva(rho_f, rho_o, V_obj, R_d)
        om_max = omeges_max.get(f_nom, 200)
        E_max  = 0.5 * I_ef * om_max**2 / 3600   # Wh
        # Massa total sistema (fluid en tub cilíndric R=0.5m, L=1m)
        V_tub  = np.pi * R_d**2 * 1.0
        m_tot  = fl['rho'] * V_tub + ob['rho'] * V_obj
        Wh_kg  = E_max / m_tot if m_tot > 0 else 0
        if Wh_kg > 0.01:
            print(f"{f_nom:>16} | {o_nom:>14} | {I_ef:>12.4f} | "
                  f"{om_max:>7.0f} | {E_max:>10.4f} | {m_tot:>12.1f} | {Wh_kg:>8.5f}")

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 14))
fig.suptitle("KILÒMETRE v10 — Estudi invers: fluid dens + objecte flotant\n"
             "De l'objectiu energètic als paràmetres de disseny",
             fontsize=12, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

# 1. I_ef vs R per cada combinació
ax1 = fig.add_subplot(gs[0, 0])
Rs_plot = np.linspace(0.1, 2.0, 100)
for f_nom, o_nom, fl, ob in casos_inv:
    rho_f = fl['rho']; rho_o = ob['rho']
    if rho_f <= rho_o: continue
    I_ef_v = Rs_plot**2 * (rho_o*V_obj + C_m*rho_f*V_obj)
    ax1.plot(Rs_plot, I_ef_v, lw=2, color=fl['color'],
             label=f'{f_nom[:4]}+{o_nom[:3]}')
ax1.set_xlabel('R [m]'); ax1.set_ylabel('I_efectiva [kg·m²]')
ax1.set_title('Inèrcia efectiva vs radi\n(inclou inèrcia afegida del fluid)')
ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)

# 2. ω necessària per E_k objectiu vs R
ax2 = fig.add_subplot(gs[0, 1])
for E_nom_p, E_p in [('100J', 100), ('1kJ', 1e3), ('10kJ', 10e3)]:
    for f_nom, o_nom, fl, ob in casos_inv[:2]:  # Hg i Ga
        rho_f = fl['rho']; rho_o = ob['rho']
        if rho_f <= rho_o: continue
        I_ef_v = Rs_plot**2 * (rho_o*V_obj + C_m*rho_f*V_obj)
        om_v   = np.sqrt(2*E_p / I_ef_v)
        om_v[om_v > omeges_max.get(f_nom,200)] = np.nan
        ax2.plot(Rs_plot, om_v, lw=1.5, color=fl['color'],
                 ls='-' if 'Hg' in f_nom or 'Merc' in f_nom else '--',
                 label=f'{f_nom[:4]} {E_nom_p}')
ax2.axhline(150, color='red', lw=1, ls=':', label='ω_max Hg')
ax2.axhline(200, color='purple', lw=1, ls=':', label='ω_max Ga')
ax2.set_xlabel('R [m]'); ax2.set_ylabel('ω necessària [rad/s]')
ax2.set_title('ω per assolir E_k objectiu\n(NaN = supera límit material/fluid)')
ax2.set_ylim(0, 400); ax2.legend(fontsize=7); ax2.grid(True, alpha=0.3)

# 3. Mapa E_k vs R i ω (Mercuri+Heli)
ax3 = fig.add_subplot(gs[0, 2])
rho_f_hg = fluids['Mercuri']['rho']
rho_o_he = objectes['Heli (globus)']['rho']
Rs_m   = np.linspace(0.1, 2.0, 60)
oms_m  = np.linspace(10, 150, 60)
RR, OO = np.meshgrid(Rs_m, oms_m)
I_ef_m = RR**2 * (rho_o_he*V_obj + C_m*rho_f_hg*V_obj)
E_k_m  = 0.5 * I_ef_m * OO**2 / 3600   # Wh
cf = ax3.contourf(RR, OO, E_k_m, levels=20, cmap='viridis')
plt.colorbar(cf, ax=ax3, label='E_k [Wh]')
# Contorns objectiu
cs = ax3.contour(RR, OO, E_k_m,
                 levels=[0.028, 0.28, 2.8, 28],
                 colors='white', linewidths=1.5)
ax3.clabel(cs, fmt='%.2g Wh', fontsize=8)
ax3.set_xlabel('R [m]'); ax3.set_ylabel('ω [rad/s]')
ax3.set_title('Mercuri+Heli: mapa E_k(R,ω)\n[Wh]')

# 4. Potència pic vs E_k (diagrama Ragone simplificat)
ax4 = fig.add_subplot(gs[1, 0])
oms_r = np.linspace(10, 150, 50)
for f_nom, o_nom, fl, ob in casos_inv:
    rho_f = fl['rho']; rho_o = ob['rho']
    if rho_f <= rho_o: continue
    Fn    = F_neta(rho_f, rho_o, V_obj)
    for R_r in [0.3, 0.5, 1.0]:
        I_ef_r = R_r**2 * (rho_o*V_obj + C_m*rho_f*V_obj)
        E_k_r  = 0.5 * I_ef_r * oms_r**2 / 3600
        P_r    = alfa_gen * abs(Fn) * R_r * oms_r
        ax4.loglog(E_k_r, P_r, lw=1.5, color=fl['color'],
                   label=f'{f_nom[:4]} R={R_r}m')
ax4.set_xlabel('E_k [Wh]'); ax4.set_ylabel('P_pic [W]')
ax4.set_title('Diagrama Ragone simplificat\n(cada línia = fluid+R fix, ω variable)')
ax4.legend(fontsize=7); ax4.grid(True, alpha=0.3, which='both')

# 5. Net elèctric per volta vs omega (Mercuri vs Gallium vs Salmuera)
ax5 = fig.add_subplot(gs[1, 1])
oms_net = np.linspace(5, 150, 100)
for f_nom, o_nom, fl, ob in casos_inv:
    rho_f = fl['rho']; rho_o = ob['rho']
    if rho_f <= rho_o: continue
    nets = []
    for om_n in oms_net:
        _, _, net, _, _ = balanc_volta(rho_f, rho_o, V_obj,
                                       0.5, om_n, fl['mu'])
        nets.append(net)
    ax5.plot(oms_net, nets, lw=2, color=fl['color'],
             label=f'{f_nom[:8]}+{o_nom[:3]}')
ax5.axhline(0, color='black', lw=1)
ax5.set_xlabel('ω [rad/s]'); ax5.set_ylabel('Net elèctric [J/volta]')
ax5.set_title('Net elèctric per volta vs ω\n(R=0.5m, V=10L)')
ax5.legend(fontsize=8); ax5.grid(True, alpha=0.3)

# 6. Comparació densitat energètica: km v10 vs tecnologies
ax6 = fig.add_subplot(gs[1, 2])
# Regions tecnologies (Wh/kg vs W/kg)
tecs = [
    ('Li-ion',      (50, 300),   (100, 500),   'blue',   0.25),
    ('Supercap',    (1, 20),     (1000,10000), 'orange', 0.25),
    ('Volant CFRP', (10, 100),   (500, 5000),  'green',  0.25),
    ('Pb-àcid',     (20, 50),    (50, 200),    'brown',  0.25),
]
for nom_t, E_r, P_r, col, al in tecs:
    ax6.fill_between(E_r, P_r[0], P_r[1], alpha=al, color=col, label=nom_t)

# Kilòmetre v10 (Hg+He, R=0.5m, omega variable)
rho_f_h = fluids['Mercuri']['rho']
rho_o_h = objectes['Heli (globus)']['rho']
V_tub6  = np.pi * 0.5**2 * 1.0
m_tot6  = rho_f_h * V_tub6 + rho_o_h * V_obj
oms6    = np.linspace(10, 150, 50)
I_ef6   = 0.5**2 * (rho_o_h*V_obj + C_m*rho_f_h*V_obj)
Ek6     = 0.5 * I_ef6 * oms6**2 / 3600 / m_tot6   # Wh/kg
Fn6     = F_neta(rho_f_h, rho_o_h, V_obj)
Pp6     = alfa_gen * abs(Fn6) * 0.5 * oms6 / m_tot6  # W/kg
ax6.plot(Ek6, Pp6, 'r-', lw=3, label='Km v10 Hg+He R=0.5m')

# Gallium+EPS
rho_f_g = fluids['Gallium']['rho']
rho_o_g = objectes['Escuma EPS']['rho']
m_tot_g = rho_f_g*V_tub6 + rho_o_g*V_obj
I_ef_g  = 0.5**2*(rho_o_g*V_obj + C_m*rho_f_g*V_obj)
oms_g   = np.linspace(10, 200, 50)
Ek_g    = 0.5*I_ef_g*oms_g**2/3600/m_tot_g
Fn_g    = F_neta(rho_f_g, rho_o_g, V_obj)
Pp_g    = alfa_gen*abs(Fn_g)*0.5*oms_g/m_tot_g
ax6.plot(Ek_g, Pp_g, color='purple', lw=3, label='Km v10 Ga+EPS R=0.5m')

ax6.set_xscale('log'); ax6.set_yscale('log')
ax6.set_xlabel('Densitat energètica [Wh/kg]')
ax6.set_ylabel('Densitat de potència [W/kg]')
ax6.set_title('Diagrama de Ragone\nposicionament tecnològic')
ax6.legend(fontsize=7); ax6.grid(True, alpha=0.3, which='both')

plt.savefig('/mnt/user-data/outputs/kilometre_v10_invers_fluid.png',
            dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v10_invers_fluid.py  |  .png")
