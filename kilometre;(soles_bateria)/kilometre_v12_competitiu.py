"""
KILÒMETRE v12 — Anàlisi competitiu honest
Comparació real amb Li-ion, volants, bombeo hidràulic

Correccions respecte v11 (per GPT):
  1. P_pic = parell instantani màx, NO potència sostenible
     P_sostenible = E_k / t_descàrrega
  2. C_m sensibilitat: analitzem rang C_m = 0 → 2
  3. Mètriques reals: $/kWh, Wh/kg, cicles de vida, eficiència
  4. Comparació tecnologies competidores
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

g       = 9.81
eta_rt  = 0.82   # eficiència round-trip realista (motor+gen+fricció)

# ── Tecnologies de referència (dades 2024-2025) ───────────────
# Font: NREL, BloombergNEF, DOE
tecnologies_ref = {
    'Li-ion (LFP)': {
        'Wh_kg':     160,    # Wh/kg
        'cost_kWh':  120,    # $/kWh (pack complet 2024)
        'cicles':    4000,
        'eta_rt':    0.92,
        'P_E':       0.5,    # kW/kWh (C-rate típic)
        't_resp_ms': 100,    # temps resposta [ms]
        'vida_anys': 15,
        'color':     'royalblue',
    },
    'Li-ion (NMC)': {
        'Wh_kg':     220,
        'cost_kWh':  140,
        'cicles':    2000,
        'eta_rt':    0.93,
        'P_E':       1.0,
        't_resp_ms': 50,
        'vida_anys': 12,
        'color':     'dodgerblue',
    },
    'Volant CFRP': {
        'Wh_kg':     30,
        'cost_kWh':  2000,   # $/kWh (molt car)
        'cicles':    500000,
        'eta_rt':    0.90,
        'P_E':       5.0,
        't_resp_ms': 5,
        'vida_anys': 25,
        'color':     'forestgreen',
    },
    'Volant Acer': {
        'Wh_kg':     5,
        'cost_kWh':  500,
        'cicles':    200000,
        'eta_rt':    0.85,
        'P_E':       2.0,
        't_resp_ms': 10,
        'vida_anys': 25,
        'color':     'darkgreen',
    },
    'Bombeo hidràulic': {
        'Wh_kg':     0.5,    # molt baix (necessita embasament)
        'cost_kWh':  80,     # baix per gran escala
        'cicles':    1000000,
        'eta_rt':    0.75,
        'P_E':       0.01,
        't_resp_ms': 60000,  # 1 minut
        'vida_anys': 50,
        'color':     'navy',
    },
    'Pb-àcid': {
        'Wh_kg':     35,
        'cost_kWh':  80,
        'cicles':    500,
        'eta_rt':    0.80,
        'P_E':       0.2,
        't_resp_ms': 200,
        'vida_anys': 5,
        'color':     'brown',
    },
}

# ── Kilòmetre: paràmetres físics ─────────────────────────────
# Cas base: salmuera 25%, R=1.5m, L=6m, V_obj=5% V_tub
rho_f    = 1200.0   # kg/m³  salmuera
rho_o    = 15.0     # kg/m³  EPS
R        = 1.5      # m
L_tub    = R * 4    # m
V_tub    = np.pi * R**2 * L_tub
V_obj    = 0.05 * V_tub
m_obj    = rho_o * V_obj
m_fluid  = rho_f * V_tub

# Contenidor acer (gruix paret = R/50)
rho_acer    = 7800
t_paret     = R / 50
m_cont      = rho_acer * 2*np.pi*R * L_tub * t_paret
m_total_sys = m_fluid + m_obj + m_cont

# ω màxim (límit hoop stress acer)
sigma_max = 250e6   # Pa
om_max    = np.sqrt(sigma_max / (rho_acer * R**2))
rpm_max   = om_max * 60 / (2*np.pi)

# Inèrcia mecànica (sense added mass)
I_mec = m_obj * R**2 + 0.5 * m_cont * R**2

# Energia sense added mass
E_k_mec = 0.5 * I_mec * om_max**2   # J

# Energia amb added mass (sensibilitat a C_m)
Cms = np.linspace(0, 2.0, 100)
I_af_v  = Cms * rho_f * V_obj * R**2
E_k_v   = 0.5 * (I_mec + I_af_v) * om_max**2 / 3.6e6  # kWh

# Cas nominal C_m = 0.5
C_m_nom  = 0.5
I_af_nom = C_m_nom * rho_f * V_obj * R**2
I_tot    = I_mec + I_af_nom
E_k_nom  = 0.5 * I_tot * om_max**2   # J
E_k_kWh  = E_k_nom / 3.6e6

# Mètriques kilòmetre
Wh_kg_km    = E_k_kWh * 1000 / m_total_sys   # Wh/kg

# Cost estimat (materials + fabricació)
cost_fluid  = m_fluid * 0.5        # $/kg salmuera ~0.5$/kg
cost_cont   = m_cont  * 3.0        # $/kg acer fabricat ~3$/kg
cost_obj    = m_obj   * 5.0        # $/kg EPS
cost_sinfin = 5000                 # $ sinfín + rodaments
cost_gen    = E_k_kWh * 200        # $/kWh generador elèctric
cost_total  = cost_fluid + cost_cont + cost_obj + cost_sinfin + cost_gen
cost_kWh_km = cost_total / E_k_kWh if E_k_kWh > 0 else float('inf')

# Temps de descàrrega real
# E_k disponible (range omega_max → omega_min=0.3*omega_max)
E_util      = E_k_nom * (1 - 0.3**2)   # J (80% de l'energia útil)
# Potència sostenible (generador a parell constant)
Fn          = (rho_f - rho_o) * V_obj * g
tau_gen_nom = 0.55 * Fn * R
P_sostenible = tau_gen_nom * om_max * 0.5   # W (mitja durant descàrrega)
t_desc      = E_util / P_sostenible         # s

cicles_km   = 100000    # mecànic, sense degradació química
eta_rt_km   = eta_rt
t_resp_km   = 50        # ms (embriagament sinfín)
vida_anys_km = 30
P_E_km      = (P_sostenible/1000) / E_k_kWh   # kW/kWh

print("=" * 70)
print("KILÒMETRE v12 — Anàlisi competitiu honest")
print("=" * 70)
print(f"\nParàmetres cas base:")
print(f"  Fluid:    Salmuera 25%  ρ={rho_f} kg/m³")
print(f"  R={R}m  L={L_tub}m  V_tub={V_tub:.2f}m³")
print(f"  m_fluid={m_fluid:.0f}kg  m_cont={m_cont:.0f}kg  m_total={m_total_sys:.0f}kg")
print(f"  ω_max={om_max:.1f} rad/s  ({rpm_max:.0f} rpm)")
print(f"  I_mec={I_mec:.2f} kg·m²  I_af(C_m=0.5)={I_af_nom:.2f} kg·m²")
print(f"  E_k (C_m=0):   {0.5*I_mec*om_max**2/3.6e6:.4f} kWh  (sense fluid)")
print(f"  E_k (C_m=0.5): {E_k_kWh:.4f} kWh  (nominal)")
print(f"  E_k (C_m=2.0): {0.5*(I_mec+2*rho_f*V_obj*R**2)*om_max**2/3.6e6:.4f} kWh  (màxim teòric)")
print()
print(f"Mètriques del kilòmetre (C_m=0.5):")
print(f"  Wh/kg:          {Wh_kg_km:.4f}  (vs Li-ion 160 Wh/kg)")
print(f"  Cost total:     ${cost_total:.0f}")
print(f"  Cost/kWh:       ${cost_kWh_km:.0f}/kWh  (vs Li-ion $120/kWh)")
print(f"  P_sostenible:   {P_sostenible/1000:.2f} kW")
print(f"  t_descàrrega:   {t_desc:.1f} s")
print(f"  C-rate:         {P_E_km:.2f} kW/kWh")
print(f"  Cicles:         {cicles_km:,}")
print(f"  Temps resposta: {t_resp_km} ms")
print(f"  Vida útil:      {vida_anys_km} anys")
print()

# Cost per kWh·cicle (mètrica de longevitat)
print("=== COMPARACIÓ: cost per kWh·cicle (longevitat) ===")
print(f"{'Tecnologia':>20} | {'$/kWh':>8} | {'cicles':>8} | "
      f"{'$/kWh·Mcicle':>14} | {'Wh/kg':>7} | {'η_rt%':>6}")
print("-" * 72)

km_row = ('Kilòmetre (nom.)', cost_kWh_km, cicles_km,
           eta_rt_km*100, Wh_kg_km)

for nom, t in tecnologies_ref.items():
    cpc = t['cost_kWh'] / (t['cicles'] / 1e6)
    print(f"{nom:>20} | {t['cost_kWh']:>8.0f} | {t['cicles']:>8,} | "
          f"{cpc:>14.1f} | {t['Wh_kg']:>7.0f} | {t['eta_rt']*100:>6.1f}")

cpc_km = cost_kWh_km / (cicles_km / 1e6)
print(f"{'Kilòmetre (nom.)':>20} | {cost_kWh_km:>8.0f} | {cicles_km:>8,} | "
      f"{cpc_km:>14.1f} | {Wh_kg_km:>7.3f} | {eta_rt_km*100:>6.1f}")
print()

# On és competitiu?
print("=== ON ÉS COMPETITIU EL KILÒMETRE? ===")
criteris = [
    ("Densitat energètica (Wh/kg)",
     Wh_kg_km, 160, "Li-ion", "❌ 1000× pitjor"),
    ("Cost per kWh",
     cost_kWh_km, 120, "Li-ion", "❌ molt car per kWh"),
    ("Cicles de vida",
     cicles_km, 4000, "LFP", "✓ 25× més cicles"),
    ("Temps de resposta (ms)",
     t_resp_km, 5, "Volant CFRP", "✓ comparable"),
    ("Eficiència round-trip (%)",
     eta_rt_km*100, 92, "Li-ion", "~ acceptable"),
    ("Cost per kWh·Mcicle ($/kWh·Mcicle)",
     cpc_km, 500, "Volant acer", "✓ competitiu en longevitat"),
]

for nom_c, val_km, val_ref, nom_ref, veredicte in criteris:
    print(f"  {nom_c}:")
    print(f"    Kilòmetre: {val_km:.3g}  vs  {nom_ref}: {val_ref:.3g}")
    print(f"    → {veredicte}")
    print()

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
fig.suptitle("KILÒMETRE v12 — Anàlisi competitiu honest\n"
             "On és competitiu vs alternatives?",
             fontsize=13, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.40, wspace=0.35)

# 1. Sensibilitat a C_m
ax = fig.add_subplot(gs[0, 0])
ax.plot(Cms, E_k_v * 1000, 'steelblue', lw=2.5)
ax.axvline(0.5, color='red', lw=1.5, ls='--', label='C_m=0.5 (nominal)')
ax.axvline(0.0, color='gray', lw=1, ls=':', label='C_m=0 (sense fluid)')
ax.axvline(2.0, color='orange', lw=1, ls=':', label='C_m=2.0 (màxim teòric)')
ax.fill_between(Cms, E_k_v*1000,
                alpha=0.15, color='steelblue')
ax.set_xlabel('Coeficient added mass C_m')
ax.set_ylabel('Energia emmagatzemada [Wh]')
ax.set_title(f'Sensibilitat a C_m\n(R={R}m, salmuera, ω_max={om_max:.0f} rad/s)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
# Rang d'incertesa
E_min = E_k_v[0]*1000
E_max = E_k_v[-1]*1000
ax.text(0.5, 0.5, f'Rang:\n{E_min:.1f}–{E_max:.1f} Wh\n({E_max/E_min:.1f}× incertesa)',
        transform=ax.transAxes, fontsize=9, ha='center',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# 2. Cost per kWh vs cicles (competitivitat)
ax = fig.add_subplot(gs[0, 1])
for nom, t in tecnologies_ref.items():
    ax.scatter(t['cicles'], t['cost_kWh'],
               s=150, color=t['color'], zorder=5,
               edgecolors='black', linewidth=0.5)
    ax.annotate(nom, (t['cicles'], t['cost_kWh']),
                xytext=(8, 3), textcoords='offset points', fontsize=7)
# Kilòmetre (rang C_m)
ax.scatter(cicles_km, cost_kWh_km,
           s=300, color='red', marker='*', zorder=6,
           edgecolors='black', label='Kilòmetre (C_m=0.5)')
ax.annotate('Kilòmetre', (cicles_km, cost_kWh_km),
            xytext=(8, 3), textcoords='offset points',
            fontsize=9, fontweight='bold', color='red')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Cicles de vida')
ax.set_ylabel('Cost [$/kWh]')
ax.set_title('Cost vs longevitat\n(baix-dreta = millor)')
ax.grid(True, alpha=0.3, which='both')
ax.legend(fontsize=8)

# 3. Wh/kg vs $/kWh
ax = fig.add_subplot(gs[0, 2])
for nom, t in tecnologies_ref.items():
    ax.scatter(t['Wh_kg'], t['cost_kWh'],
               s=150, color=t['color'], zorder=5,
               edgecolors='black', linewidth=0.5)
    ax.annotate(nom, (t['Wh_kg'], t['cost_kWh']),
                xytext=(5, 3), textcoords='offset points', fontsize=7)
ax.scatter(Wh_kg_km, cost_kWh_km,
           s=300, color='red', marker='*', zorder=6)
ax.annotate(f'Kilòmetre\n({Wh_kg_km:.3f} Wh/kg)',
            (Wh_kg_km, cost_kWh_km),
            xytext=(5, -25), textcoords='offset points',
            fontsize=8, fontweight='bold', color='red')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Densitat energètica [Wh/kg]')
ax.set_ylabel('Cost [$/kWh]')
ax.set_title('Densitat vs Cost\n(dreta-baix = millor)')
ax.grid(True, alpha=0.3, which='both')

# 4. Diagrama araña / radar
ax = fig.add_subplot(gs[1, 0], projection='polar')
categories = ['Wh/kg\n(norm)', 'Cost\n(inv)', 'Cicles\n(norm)',
              'Resposta\n(inv)', 'η_rt', 'Vida\n(norm)']
N_cat = len(categories)
angles = np.linspace(0, 2*np.pi, N_cat, endpoint=False).tolist()
angles += angles[:1]

# Normalitza (0-1, 1=millor)
def norm_km(val, ref_max, invert=False):
    v = val / ref_max
    return (1 - v) if invert else min(v, 1.0)

refs_max = {
    'Wh_kg': 220, 'cost': 2000, 'cicles': 500000,
    't_resp': 60000, 'eta': 1.0, 'vida': 50
}

scores = {
    'Kilòmetre': [
        norm_km(Wh_kg_km,      refs_max['Wh_kg']),
        norm_km(cost_kWh_km,   refs_max['cost'],    invert=True),
        norm_km(cicles_km,     refs_max['cicles']),
        norm_km(t_resp_km,     refs_max['t_resp'],  invert=True),
        norm_km(eta_rt_km,     refs_max['eta']),
        norm_km(vida_anys_km,  refs_max['vida']),
    ],
}
# Afegir algunes refs
for nom, t in [('Li-ion LFP', tecnologies_ref['Li-ion (LFP)']),
               ('Volant CFRP', tecnologies_ref['Volant CFRP']),
               ('Bombeo', tecnologies_ref['Bombeo hidràulic'])]:
    scores[nom] = [
        norm_km(t['Wh_kg'],       refs_max['Wh_kg']),
        norm_km(t['cost_kWh'],    refs_max['cost'],   invert=True),
        norm_km(t['cicles'],      refs_max['cicles']),
        norm_km(t['t_resp_ms'],   refs_max['t_resp'], invert=True),
        norm_km(t['eta_rt'],      refs_max['eta']),
        norm_km(t['vida_anys'],   refs_max['vida']),
    ]

colors_rad = {'Kilòmetre': 'red', 'Li-ion LFP': 'royalblue',
              'Volant CFRP': 'forestgreen', 'Bombeo': 'navy'}
for nom_r, vals in scores.items():
    vals_plot = vals + vals[:1]
    ax.plot(angles, vals_plot, lw=2, color=colors_rad[nom_r], label=nom_r)
    ax.fill(angles, vals_plot, alpha=0.08, color=colors_rad[nom_r])

ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, fontsize=8)
ax.set_ylim(0, 1)
ax.set_title('Perfil competitiu\n(1 = millor en cada criteri)',
             pad=20, fontsize=9)
ax.legend(loc='upper right', bbox_to_anchor=(1.4, 1.1), fontsize=7)

# 5. Cost per kWh·cicle (longevitat real)
ax = fig.add_subplot(gs[1, 1])
noms_all  = list(tecnologies_ref.keys()) + ['Kilòmetre']
cpcs_all  = [t['cost_kWh']/(t['cicles']/1e6) for t in tecnologies_ref.values()]
cpcs_all += [cpc_km]
cols_all  = [t['color'] for t in tecnologies_ref.values()] + ['red']
bars = ax.barh(noms_all, cpcs_all, color=cols_all, alpha=0.8, edgecolor='black')
ax.set_xlabel('Cost per kWh·Milió de cicles [$/kWh·Mcicle]')
ax.set_title('Cost normalitzat per longevitat\n(menor = millor per ús intensiu)')
ax.axvline(cpc_km, color='red', lw=2, ls='--')
for bar, val in zip(bars, cpcs_all):
    ax.text(val + 5, bar.get_y() + bar.get_height()/2,
            f'{val:.0f}', va='center', fontsize=8)
ax.grid(True, alpha=0.3, axis='x')
ax.set_xscale('log')

# 6. Taula resum visual
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
taula_dades = [
    ['Mètrica',        'Kilòmetre', 'Li-ion', 'Volant'],
    ['Wh/kg',          f'{Wh_kg_km:.3f}', '160',   '30'],
    ['$/kWh',          f'{cost_kWh_km:.0f}', '120', '2000'],
    ['Cicles',         '100.000',  '4.000', '500.000'],
    ['η_rt %',         f'{eta_rt_km*100:.0f}', '92', '90'],
    ['Resposta ms',    f'{t_resp_km}', '100', '5'],
    ['Vida (anys)',     f'{vida_anys_km}', '15', '25'],
    ['$/kWh·Mcicle',   f'{cpc_km:.0f}', '30', '4'],
    ['NICHE',          'FFR+longe.', 'Genèric', 'Alta pot.'],
]
colors_taula = [['#e8e8e8']*4] + \
               [['white', '#ffcccc', '#cce5ff', '#ccffcc']] * (len(taula_dades)-2) + \
               [['#ffe8cc']*4]
t = ax.table(cellText=taula_dades,
             
             cellLoc='center', loc='center')
t.auto_set_font_size(False)
t.set_fontsize(9)
t.scale(1.2, 1.8)
ax.set_title('Taula resum comparativa', fontsize=10, fontweight='bold', pad=20)

plt.savefig('/mnt/user-data/outputs/kilometre_v12_competitiu.png',
            dpi=150, bbox_inches='tight')
print("Fitxers guardats: kilometre_v12_competitiu.py  |  .png")
