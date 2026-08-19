"""
KILÒMETRE v15 — Fluids densos: Mercuri+Heli vs Fe+Oli+Heli
Comparació completa amb salmuera de referència

Fe+Oli: ferrofluid o simplement partícules de ferro en oli
  - ρ_mescla depèn de la fracció volumètrica de Fe
  - Fe pur: 7874 kg/m³
  - Oli mineral: 870 kg/m³
  - Mescla 30% Fe + 70% oli: ρ ≈ 0.3×7874 + 0.7×870 = 2971 kg/m³
  - Mescla 60% Fe + 40% oli: ρ ≈ 0.6×7874 + 0.4×870 = 5072 kg/m³
  - Cost: molt baix (Fe ~0.5$/kg, oli ~1$/kg)
  - No tòxic, no volàtil, recuperable
  - Desavantatge: sedimentació (necessita agitació o ferrofluid)

Mercuri:
  - ρ = 13.600 kg/m³ → Δρ amb heli = 13.600 kg/m³ (màxim!)
  - Cost: ~$100/kg → $1.360.000/m³ (molt car)
  - Molt tòxic (mercuri elemental)
  - Avantatge: liquid a temperatura ambient, no solidifica
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

g = 9.81; C_m = 0.5; alfa = 0.55; eta = 0.90
rho_acer = 7800; sigma_acer = 250e6

# ── Definició de fluids ──────────────────────────────────────
fluids = {
    'Salmuera 25%\n(referència)': {
        'rho':      1200,
        'mu':       1.5e-3,
        'cost_m3':  0.6,       # $/m³ (sal + aigua)
        'toxic':    False,
        'color':    'steelblue',
        'nom':      'Sal',
    },
    'Fe+Oli 30%Fe\n(pràctic)': {
        'rho':      int(0.30*7874 + 0.70*870),  # ~2671 kg/m³
        'mu':       50e-3,     # Pa·s (oli viscós amb Fe)
        'cost_m3':  150,       # $/m³ (Fe~0.5$/kg × 2000kg + oli)
        'toxic':    False,
        'color':    'peru',
        'nom':      'Fe30',
    },
    'Fe+Oli 60%Fe\n(dens)': {
        'rho':      int(0.60*7874 + 0.40*870),  # ~5072 kg/m³
        'mu':       200e-3,    # Pa·s (molt viscós)
        'cost_m3':  300,       # $/m³
        'toxic':    False,
        'color':    'saddlebrown',
        'nom':      'Fe60',
    },
    'Gallium líquid\n(industrial)': {
        'rho':      6095,
        'mu':       2e-3,
        'cost_m3':  6095*220,  # $220/kg → $1.34M/m³ !!!
        'toxic':    False,
        'color':    'purple',
        'nom':      'Ga',
    },
    'Mercuri\n(tòxic, ref)': {
        'rho':      13600,
        'mu':       1.5e-3,
        'cost_m3':  13600*100, # $100/kg → $1.36M/m³
        'toxic':    True,
        'color':    'red',
        'nom':      'Hg',
    },
}

objectes = {
    'Heli (ρ=0.164)': {'rho': 0.164, 'cost_m3': 50,  'nom': 'He'},
    'Buit (ρ≈0)':     {'rho': 0.001, 'cost_m3': 200, 'nom': 'Buit'},
}

# ── Funció d'anàlisi ─────────────────────────────────────────
def analitza_complet(R, L_mult, rho_f, rho_o, mu_f,
                     cost_m3_f, cost_m3_o, frac_vol=0.05):
    L      = R * L_mult
    V_tub  = np.pi * R**2 * L
    V_obj  = frac_vol * V_tub
    m_obj  = rho_o * V_obj
    m_fl   = rho_f * V_tub

    t_p    = max(R/50, 0.005)
    m_cont = rho_acer * 2*np.pi*R * L * t_p
    I_mec  = m_obj*R**2 + 0.5*m_cont*R**2
    I_af   = C_m * rho_f * V_obj * R**2
    I_tot  = I_mec + I_af

    om_max = np.sqrt(sigma_acer / (rho_acer * R**2))
    Fn     = (rho_f - rho_o) * V_obj * g
    E_k    = 0.5 * I_tot * om_max**2
    P_nom  = alfa * Fn * R * om_max * 0.5
    t_sol  = E_k * 0.8 / P_nom if P_nom > 0 else 0

    # Fricció viscosa (Stokes)
    r_eq   = (3*V_obj/(4*np.pi))**(1/3)
    beta   = 6*np.pi*mu_f*r_eq*R / 1.0   # tau_fric/omega
    W_fric_volta = beta * om_max * 2*np.pi

    # Balanç net per volta
    W_gen  = alfa * abs(Fn) * R * 2.0
    W_mot  = W_fric_volta + W_gen
    net_v  = W_gen*eta - W_mot/eta

    # Costs
    cost_f = cost_m3_f * V_tub
    cost_o = cost_m3_o * V_obj
    cost_c = m_cont * 3.0
    cost_g = E_k/3.6e6 * 200
    cost_t = cost_f + cost_o + cost_c + cost_g + 50000
    cost_kWh = cost_t / (E_k/3.6e6) if E_k > 0 else 0

    m_tot  = m_fl + m_obj + m_cont
    Wh_kg  = E_k/3.6e3 / m_tot if m_tot > 0 else 0

    return {
        'R': R, 'L': L, 'V_tub': V_tub,
        'Drho': rho_f - rho_o,
        'Fn_MN': Fn/1e6,
        'I_tot': I_tot, 'I_af': I_af,
        'om_max': om_max, 'rpm': om_max*60/(2*np.pi),
        'E_MWh': E_k/3.6e9,
        'P_MW': P_nom/1e6,
        't_h': t_sol/3600,
        'Wh_kg': Wh_kg,
        'cost_kWh': cost_kWh,
        'cost_tot_M': cost_t/1e6,
        'cost_fluid_M': cost_f/1e6,
        'W_fric': W_fric_volta,
        'net_volta': net_v,
        'm_tot_t': m_tot/1000,
        'mu': mu_f,
    }

# ── Taula comparativa R=10m, L=40m ───────────────────────────
print("=" * 90)
print("KILÒMETRE v15 — Fluids densos: Mercuri+Heli vs Fe+Oli+Heli")
print("Cas base: R=10m, L=40m (L=4R)")
print("=" * 90)

print(f"\n{'Fluid':>22} {'Obj':>5} {'Δρ':>7} {'Fn(MN)':>8} "
      f"{'E(MWh)':>8} {'P(MW)':>8} {'t(h)':>6} "
      f"{'Wh/kg':>7} {'$/kWh':>8} {'$fluid(M)':>10}")
print("-" * 95)

R_ref = 10.0; L_mult_ref = 4
resultats_ref = {}

for f_nom, fl in fluids.items():
    for o_nom, ob in objectes.items():
        if fl['rho'] <= ob['rho']:
            continue
        r = analitza_complet(R_ref, L_mult_ref,
                             fl['rho'], ob['rho'],
                             fl['mu'], fl['cost_m3'], ob['cost_m3'])
        key = f"{fl['nom']}+{ob['nom']}"
        resultats_ref[key] = {**r, 'color': fl['color'],
                               'f_nom': f_nom, 'o_nom': o_nom}
        tox = " ⚠TÒXIC" if fl['toxic'] else ""
        print(f"{f_nom.replace(chr(10),' '):>22} {ob['nom']:>5} "
              f"{r['Drho']:>7.0f} {r['Fn_MN']:>8.3f} "
              f"{r['E_MWh']:>8.2f} {r['P_MW']:>8.2f} {r['t_h']:>6.2f} "
              f"{r['Wh_kg']:>7.4f} {r['cost_kWh']:>8.0f} "
              f"{r['cost_fluid_M']:>10.2f}{tox}")

# ── Factor d'escala vs salmuera ───────────────────────────────
print(f"\n{'─'*90}")
print("FACTORS D'ESCALA vs Salmuera+Buit (referència)")
r_ref_sal = resultats_ref['Sal+Buit']
print(f"{'Cas':>20} | {'×E_k':>7} | {'×P_nom':>7} | {'×Fn':>7} | "
      f"{'×cost_fluid':>12} | {'×cost/kWh':>10}")
print("-" * 70)
for key, r in resultats_ref.items():
    if 'Buit' not in key:
        continue
    fE  = r['E_MWh']       / r_ref_sal['E_MWh']
    fP  = r['P_MW']        / r_ref_sal['P_MW']
    fFn = r['Fn_MN']       / r_ref_sal['Fn_MN']
    fCf = r['cost_fluid_M']/ r_ref_sal['cost_fluid_M']
    fCk = r['cost_kWh']    / r_ref_sal['cost_kWh']
    print(f"{key:>20} | {fE:>7.2f}× | {fP:>7.2f}× | {fFn:>7.2f}× | "
          f"{fCf:>12.1f}× | {fCk:>10.2f}×")

# ── Escombrat R per cada fluid (objecte=Buit) ────────────────
print(f"\n{'─'*90}")
print("ESCOMBRAT R: 1m → 50m (objecte=Buit, L=4R)")
Rs = np.logspace(0, np.log10(50), 80)

resultats_R = {}
for f_nom, fl in fluids.items():
    key = fl['nom']
    Es, Ps, ts, costs = [], [], [], []
    for R in Rs:
        r = analitza_complet(R, 4, fl['rho'], 0.001,
                             fl['mu'], fl['cost_m3'], 200)
        Es.append(r['E_MWh'])
        Ps.append(r['P_MW'])
        ts.append(r['t_h'])
        costs.append(r['cost_kWh'])
    resultats_R[key] = {
        'Rs': Rs, 'Es': np.array(Es), 'Ps': np.array(Ps),
        'ts': np.array(ts), 'costs': np.array(costs),
        'color': fl['color'], 'nom': fl['nom'],
        'toxic': fl['toxic'],
    }

# ── Flota: quants per 220MW / 1 dia ──────────────────────────
P_parc = 220.0   # MW
t_obj  = 24.0    # hores

print(f"\nFLOTA per {P_parc}MW / {t_obj}h = {P_parc*t_obj:.0f}MWh")
print(f"{'Fluid+Obj':>20} | {'R(m)':>5} | {'N_units':>8} | "
      f"{'E_tot(MWh)':>11} | {'cost_tot($M)':>13} | {'$/kWh':>8}")
print("-" * 75)

for R_fl in [5.0, 10.0]:
    for f_nom, fl in fluids.items():
        r1 = analitza_complet(R_fl, 4, fl['rho'], 0.001,
                              fl['mu'], fl['cost_m3'], 200)
        if r1['P_MW'] <= 0:
            continue
        N  = int(np.ceil(P_parc * t_obj / r1['E_MWh']))
        Et = N * r1['E_MWh']
        Ct = N * r1['cost_tot_M']
        Ck = Ct*1e6 / (Et*1e3)
        tox = " ⚠" if fl['toxic'] else ""
        print(f"{fl['nom']+'+Buit':>20} | {R_fl:>5.1f} | {N:>8} | "
              f"{Et:>11.1f} | {Ct:>13.1f} | {Ck:>8.0f}{tox}")
    print()

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
fig.suptitle("KILÒMETRE v15 — Comparació fluids densos\n"
             "Mercuri+Heli  vs  Fe+Oli+Heli  vs  Salmuera+Buit",
             fontsize=13, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.35)

# 1. Energia vs R per fluid
ax = fig.add_subplot(gs[0, 0])
for key, r in resultats_R.items():
    ls = '--' if fluids[[f for f in fluids if fluids[f]['nom']==key][0]]['toxic'] else '-'
    ax.loglog(r['Rs'], r['Es'], lw=2.5, color=r['color'],
              ls=ls, label=key+(' ⚠' if fluids[[f for f in fluids
              if fluids[f]['nom']==key][0]]['toxic'] else ''))
ax.axhline(1000, color='gray', lw=1, ls=':', label='1 GWh')
ax.axhline(1,    color='gray', lw=1, ls='--', label='1 MWh')
ax.set_xlabel('R [m]'); ax.set_ylabel('E_k [MWh]')
ax.set_title('Energia vs R (L=4R, obj=Buit)\n(-- = tòxic)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

# 2. Força neta Fn vs R
ax = fig.add_subplot(gs[0, 1])
for f_nom, fl in fluids.items():
    Fns = [(fl['rho']-0.001)*0.05*np.pi*R**2*4*R*g/1e6 for R in Rs]
    ls  = '--' if fl['toxic'] else '-'
    ax.loglog(Rs, Fns, lw=2.5, color=fl['color'], ls=ls, label=fl['nom'])
ax.set_xlabel('R [m]'); ax.set_ylabel('F_neta [MN]')
ax.set_title('Força neta flotació vs R\n(Fn ∝ Δρ·R³)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 3. Cost fluid vs energia (el dilema)
ax = fig.add_subplot(gs[0, 2])
for key, r in resultats_R.items():
    fl_data = fluids[[f for f in fluids if fluids[f]['nom']==key][0]]
    cost_fl = [fl_data['cost_m3']*np.pi*R**2*4*R/1e6 for R in r['Rs']]
    ls = '--' if fl_data['toxic'] else '-'
    ax.loglog(r['Es'], cost_fl, lw=2.5, color=r['color'],
              ls=ls, label=key)
ax.set_xlabel('E_k [MWh]')
ax.set_ylabel('Cost fluid [M$]')
ax.set_title('Cost del fluid vs energia\n(el dilema: Hg i Ga car!)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

# 4. $/kWh vs R
ax = fig.add_subplot(gs[1, 0])
for key, r in resultats_R.items():
    fl_data = fluids[[f for f in fluids if fluids[f]['nom']==key][0]]
    ls = '--' if fl_data['toxic'] else '-'
    ax.loglog(r['Rs'], r['costs'], lw=2.5, color=r['color'],
              ls=ls, label=key)
ax.axhline(120,  color='blue',  lw=2, ls=':', label='Li-ion $120/kWh')
ax.axhline(80,   color='navy',  lw=1, ls=':', label='Bombeo $80/kWh')
ax.axhline(2000, color='green', lw=1, ls=':', label='Volant CFRP $2000/kWh')
ax.set_xlabel('R [m]'); ax.set_ylabel('$/kWh')
ax.set_title('Cost per kWh vs R\n(horitzontal = tecnologies ref.)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

# 5. Diagrama E vs cost (Ragone econòmic)
ax = fig.add_subplot(gs[1, 1])
for key, r in resultats_R.items():
    fl_data = fluids[[f for f in fluids if fluids[f]['nom']==key][0]]
    ls = '--' if fl_data['toxic'] else '-'
    ax.loglog(r['Es'], r['costs'], lw=2.5, color=r['color'],
              ls=ls, label=key)
    # Marca R=5 i R=10
    for R_m in [5, 10]:
        idx = np.argmin(np.abs(r['Rs']-R_m))
        ax.plot(r['Es'][idx], r['costs'][idx], 'o',
                color=r['color'], ms=7)
        if key == 'Fe30':
            ax.annotate(f'R={R_m}m',
                        (r['Es'][idx], r['costs'][idx]),
                        xytext=(5,5), textcoords='offset points', fontsize=7)
ax.axhline(120, color='blue', lw=2, ls='--', label='Li-ion $120/kWh')
ax.set_xlabel('E_k [MWh]'); ax.set_ylabel('$/kWh')
ax.set_title('E vs $/kWh — tots els fluids\n(baix-dreta = millor)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

# 6. Taula resum R=10m
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
headers = ['Fluid', 'Δρ', 'E(MWh)', 'P(MW)', '$/kWh', 'Tòxic']
rows = []
for key, r in resultats_ref.items():
    if 'Buit' not in key:
        continue
    fl_data = fluids[[f for f in fluids
                      if fluids[f]['nom']==key.split('+')[0]][0]]
    rows.append([
        fl_data['nom'],
        f"{r['Drho']:.0f}",
        f"{r['E_MWh']:.1f}",
        f"{r['P_MW']:.0f}",
        f"${r['cost_kWh']:.0f}",
        '⚠ Sí' if fl_data['toxic'] else '✓ No',
    ])
t = ax.table(cellText=rows, colLabels=headers,
             cellLoc='center', loc='center')
t.auto_set_font_size(False)
t.set_fontsize(9)
t.scale(1.15, 2.0)
ax.set_title('Resum R=10m, L=40m, obj=Buit',
             fontsize=10, fontweight='bold', pad=20)

plt.savefig('/mnt/user-data/outputs/kilometre_v15_fluids_densos.png',
            dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v15_fluids_densos.py  |  .png")
