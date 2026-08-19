"""
KILÒMETRE — Estudi invers: de l'objectiu d'energia als paràmetres físics
Per cada objectiu (1 kWh, 10 kWh, 100 kWh, 1 MWh):
  - Massa necessària per diversos R i ω
  - Tensió mecànica centrífuga (límit de materials)
  - Comparació amb volants industrials
  - Niche: on competeix el kilòmetre?
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

g = 9.81

# ── Tensió centrífuga (límit de materials) ───────────────────
# σ = ρ·ω²·R²  (tensió a la perifèria d'un rotor sòlid)
# Límits orientatius:
#   Acer estructural:  σ_max ~ 400 MPa,  ρ ~ 7800 kg/m³
#   Fibra de carboni:  σ_max ~ 1500 MPa, ρ ~ 1600 kg/m³
#   Alumini 7075:      σ_max ~ 500 MPa,  ρ ~ 2800 kg/m³

materials = {
    'Acer':    {'rho': 7800, 'sigma_max': 400e6,  'color': 'steelblue'},
    'Alumini': {'rho': 2800, 'sigma_max': 500e6,  'color': 'orange'},
    'CFRP':    {'rho': 1600, 'sigma_max': 1500e6, 'color': 'green'},
}

def omega_max_material(R, mat):
    """ω màxim per no superar tensió centrífuga"""
    return np.sqrt(mat['sigma_max'] / (mat['rho'] * R**2))

def energia_max_material(m, R, mat):
    """Energia màxima limitada per material"""
    om_max = omega_max_material(R, mat)
    I = m * R**2
    return 0.5 * I * om_max**2

# ── Estudi invers ────────────────────────────────────────────
objectius_J = {
    '1 kWh':   1e3 * 3600,
    '10 kWh':  10e3 * 3600,
    '100 kWh': 100e3 * 3600,
    '1 MWh':   1e6 * 3600,
}

print("=" * 80)
print("ESTUDI INVERS: de l'objectiu energètic als paràmetres físics")
print("E_k = ½·m·R²·ω²  |  σ_centr = ρ·ω²·R²  ≤  σ_max")
print("=" * 80)

# Per cada objectiu: donat R i material, quina m cal? quina ω?
Rs_estudi = [0.5, 1.0, 1.5, 2.0, 3.0]   # m

for nom_obj, E_obj in objectius_J.items():
    print(f"\n{'─'*80}")
    print(f"OBJECTIU: {nom_obj}  ({E_obj:.2e} J)")
    print(f"{'Material':>10} | {'R(m)':>5} | {'ω_max(rad/s)':>13} | "
          f"{'rpm':>7} | {'m_min(kg)':>10} | {'σ%(%)':>7} | {'viable':>7}")
    print("-" * 80)

    for mat_nom, mat in materials.items():
        for R in Rs_estudi:
            om_max = omega_max_material(R, mat)
            # Massa mínima per assolir E_obj a ω_max
            I_need = 2 * E_obj / om_max**2
            m_need = I_need / R**2   # I = m·R²
            rpm    = om_max * 60 / (2*np.pi)
            # Tensió real (% del màxim)
            sigma_pct = 100.0  # per definició som a om_max

            viable = "✓" if m_need < 50000 and rpm < 10000 else "✗ pesat/ràpid"
            if m_need < 1000 and rpm < 3000:
                viable = "★ òptim"

            print(f"{mat_nom:>10} | {R:>5.1f} | {om_max:>13.1f} | "
                  f"{rpm:>7.0f} | {m_need:>10.0f} | {sigma_pct:>7.0f} | {viable:>7}")

# ── Cas específic: disseny recomanat ─────────────────────────
print(f"\n{'='*80}")
print("DISSENY RECOMANAT per a 1 kWh (objectiu assolible)")
print(f"{'='*80}")

E_obj = 1e3 * 3600  # 1 kWh en J
mat = materials['Acer']

# Barrida R vs ω per veure la frontera viable
Rs   = np.linspace(0.2, 3.0, 50)
oms  = np.linspace(10, 600, 50)   # rad/s (fins ~5700 rpm)
RR, OO = np.meshgrid(Rs, oms)

# Massa necessària
II_need = 2 * E_obj / OO**2
MM_need = II_need / RR**2

# Tensió centrífuga
sigma   = mat['rho'] * OO**2 * RR**2 / 1e6   # MPa
sigma_lim = mat['sigma_max'] / 1e6

# ── Comparació amb volants industrials ───────────────────────
print("""
Volants industrials de referència (Beacon Power, Amber Kinetics, etc.):
  Beacon Power 20 kWh:   ~2300 kg,  R~0.6m,  ~9000 rpm  (CFRP)
  Amber Kinetics 32 kWh: ~2700 kg,  R~0.5m,  ~3000 rpm  (acer)
  Punts forts: resposta <1s, >100.000 cicles, eficiència >85%
  Punts febles: cost elevat, contention de seguretat

Kilòmetre (avantatge potencial):
  - Sinfín: control de parell fi sense electrònica complexa
  - Massa desplaçable: modifica I dinàmicament (Quijote integrat)
  - Submersió: refredament natural, menys soroll, seguretat passiva
  - ZypyZape: 44 unitats sincronitzades → bateria distribuïda
""")

# ── Niche del kilòmetre ──────────────────────────────────────
print("=" * 80)
print("NICHE COMPETITIU DEL KILÒMETRE")
print("=" * 80)

nichos = [
    ("Densitat energètica (Wh/kg)",
     "Li-ion: 150-250", "Volant CFRP: 20-100", "Kilòmetre: 5-30",
     "❌ No competeix en densitat"),
    ("Cicles de vida",
     "Li-ion: 500-2000", "Volant: >100.000", "Kilòmetre: >100.000",
     "✓ Competeix en longevitat"),
    ("Temps de resposta",
     "Li-ion: 100ms-1s", "Volant: <10ms", "Kilòmetre: <100ms",
     "✓ Ràpid (sinfín embriagable)"),
    ("Control de parell fi",
     "Li-ion: via inversor", "Volant: difícil", "Kilòmetre: sinfín directe",
     "★ Avantatge únic"),
    ("Regulació de freqüència (FFR)",
     "Li-ion: possible", "Volant: excel·lent", "Kilòmetre+ZypyZape: natiu",
     "★ Avantatge únic"),
    ("Cost de materials",
     "Li-ion: alt (Li, Co)", "Volant CFRP: molt alt", "Kilòmetre acer: baix",
     "✓ Competeix en cost"),
    ("Seguretat (fallada)",
     "Li-ion: foc/explosió", "Volant: projectil", "Kilòmetre submergit: passiva",
     "★ Avantatge únic submergit"),
]

for niche, ref1, ref2, km, conclusio in nichos:
    print(f"\n  {niche}:")
    print(f"    {ref1}  |  {ref2}  |  {km}")
    print(f"    → {conclusio}")

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 12))
fig.suptitle("KILÒMETRE — Estudi invers: paràmetres per objectiu energètic\n"
             "Comparació amb volants industrials",
             fontsize=12, fontweight='bold')

gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

# 1. Massa necessària vs R per 1 kWh (diverses ω, acer)
ax1 = fig.add_subplot(gs[0, 0])
E_1kWh = 1e3 * 3600
mat_acer = materials['Acer']
for om_plot in [50, 100, 200, 400]:
    m_plot = 2 * E_1kWh / (om_plot**2 * Rs**2)
    # Filtrar on σ > σ_max
    sigma_v = mat_acer['rho'] * om_plot**2 * Rs**2 / 1e6
    m_plot[sigma_v > mat_acer['sigma_max']/1e6] = np.nan
    ax1.plot(Rs, m_plot, lw=2, label=f'ω={om_plot} rad/s ({om_plot*60/(2*np.pi):.0f}rpm)')
ax1.axhline(2700, color='gray', lw=1, ls='--', label='Amber 32kWh (ref)')
ax1.set_xlabel('Radi R [m]'); ax1.set_ylabel('Massa [kg]')
ax1.set_title('Massa per 1 kWh (acer)\n(NaN = tensió excedida)')
ax1.set_ylim(0, 20000); ax1.legend(fontsize=7); ax1.grid(True, alpha=0.3)

# 2. ω_max vs R per material
ax2 = fig.add_subplot(gs[0, 1])
for mat_nom, mat in materials.items():
    om_max_v = np.sqrt(mat['sigma_max'] / (mat['rho'] * Rs**2))
    rpm_max  = om_max_v * 60 / (2*np.pi)
    ax2.plot(Rs, rpm_max, lw=2, color=mat['color'], label=mat_nom)
ax2.axhline(9000, color='gray', lw=1, ls='--', label='Beacon Power (ref)')
ax2.axhline(3000, color='gray', lw=1, ls=':',  label='Amber K. (ref)')
ax2.set_xlabel('Radi R [m]'); ax2.set_ylabel('rpm màximes')
ax2.set_title('Velocitat màxima per material\n(límit tensió centrífuga)')
ax2.set_yscale('log'); ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)

# 3. Energia màxima vs massa per R=1m
ax3 = fig.add_subplot(gs[0, 2])
masses_v = np.logspace(1, 5, 100)
R_fix = 1.0
for mat_nom, mat in materials.items():
    om_max = omega_max_material(R_fix, mat)
    E_max  = 0.5 * masses_v * R_fix**2 * om_max**2 / 3600  # Wh
    ax3.loglog(masses_v, E_max, lw=2, color=mat['color'], label=mat_nom)
# Referència Li-ion
ax3.loglog(masses_v, masses_v * 200, 'k--', lw=1, label='Li-ion 200Wh/kg')
ax3.loglog(masses_v, masses_v * 50,  'k:',  lw=1, label='Li-ion 50Wh/kg')
ax3.set_xlabel('Massa [kg]'); ax3.set_ylabel('Energia màxima [Wh]')
ax3.set_title(f'Energia màxima vs massa (R={R_fix}m)\nlimitat per tensió material')
ax3.legend(fontsize=8); ax3.grid(True, alpha=0.3)

# 4. Mapa m vs R per 10 kWh (acer, ω=100 rad/s)
ax4 = fig.add_subplot(gs[1, 0])
E_10kWh = 10e3 * 3600
om_fix  = 100  # rad/s ~ 955 rpm
mm_need = 2 * E_10kWh / (om_fix**2 * Rs**2)
sigma_v = mat_acer['rho'] * om_fix**2 * Rs**2 / 1e6
viable  = sigma_v <= mat_acer['sigma_max']/1e6
ax4.plot(Rs[viable], mm_need[viable], 'steelblue', lw=2.5, label='Massa necessària')
ax4.fill_between(Rs, 0, mm_need, where=viable, alpha=0.2, color='steelblue')
ax4.axvline(Rs[viable][0] if viable.any() else 0,
            color='red', lw=1, ls='--', label='R_min viable')
ax4.set_xlabel('Radi R [m]'); ax4.set_ylabel('Massa [kg]')
ax4.set_title(f'10 kWh, acer, ω={om_fix} rad/s ({om_fix*60/(2*np.pi):.0f} rpm)')
ax4.legend(fontsize=9); ax4.grid(True, alpha=0.3)

# 5. Comparació tecnologies: potència vs energia (Ragone plot)
ax5 = fig.add_subplot(gs[1, 1:])
# Tecnologies de referència (regions aproximades)
tecnologies = [
    ('Li-ion',        (0.1, 300),   (50, 300),   'blue',   0.3),
    ('Supercap',      (100, 10000), (1, 10),     'orange', 0.3),
    ('Volant CFRP',   (1000, 5000), (10, 100),   'green',  0.3),
    ('Volant acer',   (100, 1000),  (5, 30),     'gray',   0.3),
    ('Pb-àcid',       (0.05, 100),  (20, 50),    'brown',  0.3),
]

for nom, P_range, E_range, color, alpha in tecnologies:
    ax5.fill_between(E_range, P_range[0], P_range[1],
                     alpha=alpha, color=color, label=nom)

# Kilòmetre: posició estimada
# E: 1-100 Wh/kg (depèn de R i ω)
# P: alta (sinfín ràpid)
ax5.fill_between([5, 50], [500, 500], [5000, 5000],
                 alpha=0.5, color='red', label='Kilòmetre (estimat)')
ax5.text(15, 1500, 'KILÒMETRE\n(estimat)', color='red',
         fontsize=10, fontweight='bold', ha='center')

ax5.set_xscale('log'); ax5.set_yscale('log')
ax5.set_xlabel('Densitat energètica [Wh/kg]')
ax5.set_ylabel('Densitat de potència [W/kg]')
ax5.set_title('Diagrama de Ragone — Posicionament tecnològic')
ax5.legend(fontsize=8, loc='upper right')
ax5.grid(True, alpha=0.3, which='both')

plt.savefig('/mnt/user-data/outputs/kilometre_estudi_invers.png',
            dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_estudi_invers.py  |  .png")
