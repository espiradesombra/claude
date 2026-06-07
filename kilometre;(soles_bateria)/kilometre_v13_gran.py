"""
KILÒMETRE v13 — Escala gran: D=100m, fluid dens, diferència de densitats màxima

Idea clau:
  - R gran → I_ef ∝ R⁵ → energia enorme a ω baixa
  - Δρ = ρ_fluid - ρ_obj gran → F_neta gran → parell gran → P_pic gran
  - L llarga → V_fluid gran → més inèrcia afegida → més temps de descàrrega
  - ω baixa (R gran) → seguretat mecànica, no cal CFRP

Paràmetres explorats:
  R: 10m, 25m, 50m  (D=20m, 50m, 100m)
  L: R*2, R*4, R*8  (proporcional o independent)
  Fluid: salmuera, gallium, mercuri
  Objecte: heli/aire atrapat, EPS, buit (ρ→0)

Mètriques:
  E_k [MWh]       → energia emmagatzemada
  t_desc [hores]  → temps de descàrrega a P_nom
  P_nom [MW]      → potència nominal sostenible
  $/kWh           → cost estimat
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

g       = 9.81
C_m     = 0.5
eta_rt  = 0.82
alfa_gen = 0.55
rho_acer = 7800
sigma_acer = 250e6  # Pa

# ── Fluids ────────────────────────────────────────────────────
fluids = {
    'Salmuera':  {'rho': 1200,  'mu': 1.5e-3, 'color': 'cyan',
                  'cost_m3': 0.5,   'nom': 'Sal'},
    'Gallium':   {'rho': 6095,  'mu': 2.0e-3, 'color': 'purple',
                  'cost_m3': 300000,'nom': 'Ga'},   # ~$220/kg × 6095 kg/m³
    'Mercuri':   {'rho': 13600, 'mu': 1.5e-3, 'color': 'red',
                  'cost_m3': 150000,'nom': 'Hg'},   # tòxic però il·lustratiu
}

# Objecte: buit ρ→0 (esfera buida pressuritzada) — màxim Δρ
rho_obj_casos = {
    'Buit (ρ≈0)':  0.001,   # kg/m³
    'Heli':        0.164,
    'EPS escuma':  15.0,
}

def analitza(R, L, rho_f, rho_o, mu_f, cost_m3_f,
             frac_vol=0.05, nom_fl='', nom_ob=''):
    """Calcula totes les mètriques per una configuració."""
    V_tub   = np.pi * R**2 * L
    V_obj   = frac_vol * V_tub
    m_obj   = rho_o * V_obj
    m_fluid = rho_f * V_tub

    # Contenidor: gruix paret hoop stress
    # σ = P_hidro + ρ_f·ω²·R²/2  ≈ ρ_f·ω²·R²  (dominant a alta ω)
    # Però a R gran, ω és baixa → podem usar acer estàndard
    # gruix mínim: t = ρ_f·ω²·R³ / σ_acer
    om_max  = np.sqrt(sigma_acer / (rho_acer * R**2))   # límit contenidor
    t_paret = max(rho_f * om_max**2 * R**2 * R / sigma_acer, 0.01)  # m, mínim 1cm
    m_cont  = rho_acer * 2*np.pi*R * L * t_paret
    m_total = m_fluid + m_obj + m_cont

    # Inèrcia efectiva
    I_mec   = m_obj * R**2 + 0.5 * m_cont * R**2
    I_af    = C_m * rho_f * V_obj * R**2
    I_tot   = I_mec + I_af

    # Energia
    E_k     = 0.5 * I_tot * om_max**2   # J
    E_kWh   = E_k / 3.6e6               # kWh
    E_MWh   = E_kWh / 1000              # MWh

    # Força neta i potència
    Fn      = (rho_f - rho_o) * V_obj * g
    tau_nom = alfa_gen * Fn * R
    P_nom   = tau_nom * om_max * 0.5    # W (mitja descàrrega)
    P_MW    = P_nom / 1e6               # MW

    # Temps descàrrega (80% E_k útil)
    E_util  = E_k * 0.80
    t_desc  = E_util / P_nom if P_nom > 0 else 0   # s
    t_hores = t_desc / 3600

    # Cost
    cost_fluid = m_fluid * cost_m3_f / (rho_f)  # $/kg → $/m³ ja aplicat
    cost_cont  = m_cont  * 3.0
    cost_obj   = m_obj   * 10.0
    cost_gen   = E_kWh   * 200   # $/kWh generador
    cost_tot   = cost_fluid*rho_f + cost_cont + cost_obj + cost_gen + 50000
    # Nota: cost_fluid en $/m³: cost_m3_f * V_tub
    cost_fluid_real = cost_m3_f * V_tub
    cost_tot   = cost_fluid_real + cost_cont + cost_obj + cost_gen + 50000
    cost_kWh   = cost_tot / E_kWh if E_kWh > 0 else float('inf')

    # Wh/kg
    Wh_kg   = E_kWh * 1000 / m_total if m_total > 0 else 0

    return {
        'R': R, 'L': L, 'V_tub': V_tub, 'V_obj': V_obj,
        'om_max': om_max, 'rpm': om_max*60/(2*np.pi),
        'm_total': m_total, 'm_fluid': m_fluid, 'm_cont': m_cont,
        'I_tot': I_tot, 'I_mec': I_mec, 'I_af': I_af,
        'E_MWh': E_MWh, 'E_kWh': E_kWh,
        'Fn': Fn, 'P_MW': P_MW,
        't_hores': t_hores,
        'Wh_kg': Wh_kg, 'cost_kWh': cost_kWh, 'cost_tot': cost_tot,
        'Delta_rho': rho_f - rho_o,
        'nom_fl': nom_fl, 'nom_ob': nom_ob,
        't_paret': t_paret,
    }

# ── Escombrat principal ──────────────────────────────────────
print("=" * 90)
print("KILÒMETRE v13 — Escala gran: R=10m → 50m, fluid dens, Δρ màxim")
print("=" * 90)

Rs_gran  = [5, 10, 25, 50]
Ls_mult  = [2, 4, 8]        # L = mult × R

# Taula principals
print(f"\n{'R(m)':>5} {'L(m)':>6} {'Fluid':>10} {'Objecte':>12} "
      f"{'E(MWh)':>8} {'P(MW)':>7} {'t(h)':>7} "
      f"{'rpm':>6} {'$/kWh':>8} {'t_paret(m)':>11}")
print("-" * 95)

resultats_gran = []
for R in Rs_gran:
    for L_mult in [4]:   # L = 4R per defecte
        L = R * L_mult
        for f_nom in ['Salmuera', 'Gallium']:
            fl = fluids[f_nom]
            for o_nom, rho_o in [('Buit', 0.001), ('EPS', 15.0)]:
                r = analitza(R, L, fl['rho'], rho_o,
                             fl['mu'], fl['cost_m3'],
                             nom_fl=f_nom, nom_ob=o_nom)
                resultats_gran.append(r)
                print(f"{R:>5.0f} {L:>6.0f} {f_nom:>10} {o_nom:>12} "
                      f"{r['E_MWh']:>8.3f} {r['P_MW']:>7.2f} {r['t_hores']:>7.2f} "
                      f"{r['rpm']:>6.1f} {r['cost_kWh']:>8.0f} {r['t_paret']:>11.3f}")

# ── Efecte longitud L (R=10m fix) ───────────────────────────
print(f"\n{'─'*90}")
print("Efecte de la LONGITUD L (R=10m, Salmuera, Buit)")
print(f"{'L(m)':>6} {'E(MWh)':>9} {'P(MW)':>8} {'t(h)':>8} "
      f"{'m_fluid(t)':>12} {'$/kWh':>8}")
print("-" * 60)
fl = fluids['Salmuera']
for L_m in [20, 40, 80, 160, 400, 1000]:
    r = analitza(10, L_m, fl['rho'], 0.001, fl['mu'], fl['cost_m3'])
    print(f"{L_m:>6.0f} {r['E_MWh']:>9.3f} {r['P_MW']:>8.2f} "
          f"{r['t_hores']:>8.2f} {r['m_fluid']/1000:>12.1f} {r['cost_kWh']:>8.0f}")

# ── Efecte Δρ (R=10m, L=40m) ────────────────────────────────
print(f"\n{'─'*90}")
print("Efecte de la DIFERÈNCIA DE DENSITATS Δρ (R=10m, L=40m)")
print(f"{'Fluid':>12} {'Objecte':>12} {'Δρ':>8} {'F_neta(MN)':>12} "
      f"{'E(MWh)':>8} {'P(MW)':>8} {'t(h)':>8}")
print("-" * 75)
for f_nom, fl in fluids.items():
    for o_nom, rho_o in [('Buit', 0.001), ('EPS', 15.0)]:
        r = analitza(10, 40, fl['rho'], rho_o, fl['mu'], fl['cost_m3'])
        print(f"{f_nom:>12} {o_nom:>12} {r['Delta_rho']:>8.0f} "
              f"{r['Fn']/1e6:>12.3f} {r['E_MWh']:>8.3f} "
              f"{r['P_MW']:>8.2f} {r['t_hores']:>8.2f}")

# ── Cas estrella: R=50m, L=200m, Salmuera, Buit ─────────────
print(f"\n{'='*90}")
print("CAS ESTRELLA: R=50m, L=200m, Salmuera 25%, objecte buit")
print(f"{'='*90}")
fl = fluids['Salmuera']
r_star = analitza(50, 200, fl['rho'], 0.001, fl['mu'], fl['cost_m3'])
print(f"""
  Geometria:
    Diàmetre: {2*r_star['R']:.0f}m  |  Longitud: {r_star['L']:.0f}m
    V_tub:    {r_star['V_tub']:.0f} m³  |  V_obj: {r_star['V_obj']:.0f} m³ ({r_star['V_obj']/r_star['V_tub']*100:.0f}% del tub)
    Gruix paret contenidor: {r_star['t_paret']:.3f}m

  Mecànica:
    ω_max: {r_star['om_max']:.3f} rad/s  ({r_star['rpm']:.2f} rpm)  ← molt lent!
    I_tot: {r_star['I_tot']:.2e} kg·m²
    I_mec: {r_star['I_mec']:.2e}  |  I_af: {r_star['I_af']:.2e}

  Masses:
    m_fluid:  {r_star['m_fluid']/1e6:.3f} Mt  ({r_star['m_fluid']:.2e} kg)
    m_cont:   {r_star['m_cont']/1000:.0f} t
    m_total:  {r_star['m_total']/1e6:.3f} Mt

  Energia i potència:
    E_k:      {r_star['E_MWh']:.1f} MWh  ({r_star['E_MWh']/1000:.4f} GWh)
    F_neta:   {r_star['Fn']/1e6:.2f} MN  ({r_star['Fn']/g/1000:.0f} tones equivalents)
    P_nom:    {r_star['P_MW']:.2f} MW
    t_desc:   {r_star['t_hores']:.2f} hores

  Economia:
    Cost total: ${r_star['cost_tot']/1e6:.1f}M
    $/kWh:      ${r_star['cost_kWh']:.0f}/kWh
    Wh/kg:      {r_star['Wh_kg']:.4f}
""")

# Comparació directa amb centrals de bombeo
print("  Comparació amb bombeo hidràulic equivalent:")
E_bombeo = r_star['E_MWh']
# Embasament: E = ρ·g·V·h, h=100m típic
V_emb = E_bombeo * 3.6e9 / (1000 * g * 100)
print(f"    Per {E_bombeo:.0f} MWh de bombeo (h=100m):")
print(f"    V_embasament = {V_emb:.0f} m³  ({V_emb/1e6:.3f} hm³)")
print(f"    Superfície (10m profunditat): {V_emb/10/1e4:.1f} hectàrees")

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
fig.suptitle("KILÒMETRE v13 — Escala gran (R fins 50m)\n"
             "Energia, potència i temps de descàrrega vs dimensions",
             fontsize=13, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

# Dades per gràfics: R variable, L=4R, Salmuera, Buit
Rs_plot = np.logspace(0, 2, 100)   # 1m → 100m
fl_s    = fluids['Salmuera']
E_MWh_v, P_MW_v, t_h_v, rpm_v, cst_v = [], [], [], [], []
for R_p in Rs_plot:
    r_p = analitza(R_p, R_p*4, fl_s['rho'], 0.001,
                   fl_s['mu'], fl_s['cost_m3'])
    E_MWh_v.append(r_p['E_MWh'])
    P_MW_v.append(r_p['P_MW'])
    t_h_v.append(r_p['t_hores'])
    rpm_v.append(r_p['rpm'])
    cst_v.append(r_p['cost_kWh'])

E_MWh_v = np.array(E_MWh_v)
P_MW_v  = np.array(P_MW_v)
t_h_v   = np.array(t_h_v)
rpm_v   = np.array(rpm_v)
cst_v   = np.array(cst_v)

# 1. Energia vs R
ax = fig.add_subplot(gs[0, 0])
ax.loglog(Rs_plot, E_MWh_v, 'steelblue', lw=2.5, label='Salmuera+Buit')
# Altres fluids
for f_nom, fl_p in [('Gallium', fluids['Gallium']),
                     ('Mercuri', fluids['Mercuri'])]:
    E_v2 = [analitza(R_p, R_p*4, fl_p['rho'], 0.001,
                      fl_p['mu'], fl_p['cost_m3'])['E_MWh']
             for R_p in Rs_plot]
    ax.loglog(Rs_plot, E_v2, lw=2, color=fl_p['color'],
              label=f'{f_nom}+Buit', alpha=0.8)
ax.axvline(10, color='gray', lw=1, ls='--', alpha=0.5)
ax.axvline(50, color='gray', lw=1, ls='--', alpha=0.5)
ax.axhline(1000, color='orange', lw=1, ls=':', label='1 GWh')
ax.axhline(1,    color='green',  lw=1, ls=':', label='1 MWh')
ax.set_xlabel('R [m]'); ax.set_ylabel('E_k [MWh]')
ax.set_title('Energia emmagatzemada vs R\n(L=4R, objecte buit)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')
ax.text(0.05, 0.85, 'E ∝ R³', transform=ax.transAxes,
        fontsize=12, fontweight='bold', color='steelblue')

# 2. Temps de descàrrega vs R
ax = fig.add_subplot(gs[0, 1])
ax.loglog(Rs_plot, t_h_v, 'steelblue', lw=2.5, label='Salmuera')
ax.axhline(1,   color='green',  lw=1, ls='--', label='1 hora')
ax.axhline(8,   color='orange', lw=1, ls='--', label='8 hores (jornada)')
ax.axhline(24,  color='red',    lw=1, ls='--', label='24 hores')
for R_m in [10, 50]:
    idx = np.argmin(np.abs(Rs_plot - R_m))
    ax.annotate(f'R={R_m}m\n{t_h_v[idx]:.1f}h',
                (Rs_plot[idx], t_h_v[idx]),
                xytext=(10, 10), textcoords='offset points', fontsize=8)
ax.set_xlabel('R [m]'); ax.set_ylabel('Temps descàrrega [h]')
ax.set_title('Temps de descàrrega vs R\n(a potència nominal)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 3. rpm vs R
ax = fig.add_subplot(gs[0, 2])
ax.loglog(Rs_plot, rpm_v, 'firebrick', lw=2.5)
ax.axhline(1000, color='green',  lw=1, ls='--', label='1000 rpm (industrial)')
ax.axhline(10,   color='orange', lw=1, ls='--', label='10 rpm (molt lent)')
ax.axhline(1,    color='red',    lw=1, ls='--', label='1 rpm')
for R_m in [10, 50]:
    idx = np.argmin(np.abs(Rs_plot - R_m))
    ax.annotate(f'R={R_m}m\n{rpm_v[idx]:.2f}rpm',
                (Rs_plot[idx], rpm_v[idx]),
                xytext=(10, 5), textcoords='offset points', fontsize=8)
ax.set_xlabel('R [m]'); ax.set_ylabel('rpm màximes')
ax.set_title('Velocitat de rotació vs R\n(límit hoop stress acer)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 4. Efecte Δρ fix R=10m
ax = fig.add_subplot(gs[1, 0])
Delta_rhos = np.linspace(10, 14000, 200)
E_delta, P_delta, t_delta = [], [], []
for Drho in Delta_rhos:
    rho_f_d = Drho + 0.001
    r_d = analitza(10, 40, rho_f_d, 0.001,
                   1.5e-3, 0.5*rho_f_d)
    E_delta.append(r_d['E_MWh'])
    P_delta.append(r_d['P_MW'])
    t_delta.append(r_d['t_hores'])

ax2 = ax.twinx()
ax.plot(Delta_rhos, E_delta, 'steelblue', lw=2, label='E_k (MWh)')
ax2.plot(Delta_rhos, t_delta, 'orange', lw=2, ls='--', label='t_desc (h)')
ax.axvline(200,   color='cyan',   lw=1, ls=':', label='Sal Δρ≈200')
ax.axvline(6080,  color='purple', lw=1, ls=':', label='Gallium Δρ≈6080')
ax.axvline(13600, color='red',    lw=1, ls=':', label='Hg Δρ≈13600')
ax.set_xlabel('Δρ = ρ_fluid - ρ_obj [kg/m³]')
ax.set_ylabel('E_k [MWh]', color='steelblue')
ax2.set_ylabel('t_desc [h]', color='orange')
ax.set_title('Efecte Δρ\n(R=10m, L=40m, V_obj=5%)')
lines1, labs1 = ax.get_legend_handles_labels()
lines2, labs2 = ax2.get_legend_handles_labels()
ax.legend(lines1+lines2, labs1+labs2, fontsize=7)
ax.grid(True, alpha=0.3)

# 5. E vs L (R=10m fix)
ax = fig.add_subplot(gs[1, 1])
Ls_plot = np.logspace(1, 4, 100)   # 10m → 10km
fl_s = fluids['Salmuera']
E_L, P_L, t_L = [], [], []
for L_p in Ls_plot:
    r_l = analitza(10, L_p, fl_s['rho'], 0.001,
                   fl_s['mu'], fl_s['cost_m3'])
    E_L.append(r_l['E_MWh'])
    P_L.append(r_l['P_MW'])
    t_L.append(r_l['t_hores'])

ax.loglog(Ls_plot, E_L, 'steelblue', lw=2, label='E_k (MWh)')
ax.loglog(Ls_plot, t_L, 'orange',    lw=2, ls='--', label='t_desc (h)')
ax.loglog(Ls_plot, P_L, 'green',     lw=2, ls=':',  label='P_nom (MW)')
ax.axhline(1000, color='gray', lw=1, ls='--', label='1 GWh / 1GW')
ax.set_xlabel('Longitud L [m]')
ax.set_ylabel('MWh / MW / hores')
ax.set_title('Efecte longitud L\n(R=10m, Salmuera, Buit)\nE∝L, P∝L, t∝cte!')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')
ax.text(0.05, 0.85, 'E ∝ L\nP ∝ L\nt ≈ cte!', transform=ax.transAxes,
        fontsize=11, fontweight='bold', color='darkred')

# 6. Diagrama final: E vs t_desc per totes combinacions
ax = fig.add_subplot(gs[1, 2])
# Zones d'aplicació
apps = [
    ('FFR / Regulació\nfreqüència',  (0.001,1),    (0.001,0.1), 'lightgreen'),
    ('Suport de xarxa\nintradiari',  (1,1000),     (0.1,8),     'lightyellow'),
    ('Emmagatzematge\nbulk diari',   (100,100000), (4,24),      'lightblue'),
    ('Estacional',                   (1000,1e6),   (100,8760),  'lightsalmon'),
]
for nom_a, E_r, t_r, col_a in apps:
    ax.fill_between(E_r, t_r[0], t_r[1], alpha=0.35, color=col_a)
    ax.text(np.sqrt(E_r[0]*E_r[1]), np.sqrt(t_r[0]*t_r[1]),
            nom_a, ha='center', va='center', fontsize=7,
            color='dimgray', fontweight='bold')

# Trajectòries kilòmetre
for f_nom, fl_p in fluids.items():
    E_tr, t_tr = [], []
    for R_p in np.logspace(0, 2, 50):
        r_tr = analitza(R_p, R_p*4, fl_p['rho'], 0.001,
                        fl_p['mu'], fl_p['cost_m3'])
        E_tr.append(r_tr['E_MWh'])
        t_tr.append(r_tr['t_hores'])
    ax.loglog(E_tr, t_tr, lw=2.5, color=fl_p['color'],
              label=f'Km {fl_p["nom"]}+Buit')
    # Marca R=10 i R=50
    for R_m in [1, 5, 10, 50]:
        idx = np.argmin(np.abs(np.logspace(0,2,50) - R_m))
        ax.plot(E_tr[idx], t_tr[idx], 'o',
                color=fl_p['color'], ms=6)
        if f_nom == 'Salmuera':
            ax.annotate(f'R={R_m}m', (E_tr[idx], t_tr[idx]),
                        xytext=(5,3), textcoords='offset points', fontsize=7)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('E_k [MWh]')
ax.set_ylabel('Temps descàrrega [hores]')
ax.set_title('On cau el kilòmetre gran?\nE vs t_desc per mida')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

plt.savefig('/mnt/user-data/outputs/kilometre_v13_gran.png',
            dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v13_gran.py  |  .png")
