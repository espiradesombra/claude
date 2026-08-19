"""
KILÒMETRE STIRLING — Flota de tubs xicotets
Principi: maximitzar superfície de captació tèrmica
          minimitzar cost per unitat
          moltes unitats × 120W = MW totals

Física clau:
  P_unitat ≈ Q_flux × A_superficie × η_real
  A_sup = 2π·R·L  → per tub xicot, cal molts tubs
  
  Per R=0.5m, L=2m: A = 6.28m²  → P ≈ 10W/m² × 6.28 × 0.15 ≈ 9W
  Per R=0.1m, L=0.4m: A = 0.25m² → P ≈ 10W/m² × 0.25 × 0.15 ≈ 0.4W
  
  Però: molts tubs xicotets → millor intercanvi tèrmic (Nusselt)
  Tubs xicotets: temps trànsit curt PERÒ ω gran → v_axial similar

  La densitat NO importa tant:
  W_Stirling = ΔFn × L = Δ(ρ_fluid·V_obj·g - m_obj·g) × L
             = (Δρ_fluid·V_obj + ρ_fluid·ΔV_obj) × g × L
  
  L'efecte dominant és ΔV_obj (expansió tèrmica objecte)
  → ρ_fluid importa per amplificar ΔV, però no cal mercuri
  → Salmuera/Fe+Oli és suficient
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

g = 9.81

# ── Paràmetres físics fixos ──────────────────────────────────
T_fred_K  = 278   # 5°C  aigüa profunda
T_cal_K   = 348   # 75°C font volcànica suau (Açores/Canàries)
DT        = T_cal_K - T_fred_K   # 70°C

eta_Carnot = 1 - T_fred_K/T_cal_K
eta_real   = eta_Carnot * 0.35   # 35% de Carnot (motor Stirling real)

# Fluids (ara xicotets → salmuera o Fe+Oli)
fluids = {
    'Salmuera':   {'rho': 1200, 'cost_m3': 0.5,  'color': 'steelblue'},
    'Fe+Oli 30%': {'rho': 2671, 'cost_m3': 150,  'color': 'peru'},
    'Fe+Oli 60%': {'rho': 5072, 'cost_m3': 300,  'color': 'saddlebrown'},
}

# Objecte 3 cambres (escalat per R)
# V_gas_fred = frac_gas × V_tub
# Expansió tèrmica gas: ΔV/V = ΔT/T_fred (gas ideal)
frac_gas   = 0.03   # 3% del volum tub → cambra gas
frac_liq   = 0.02   # 2% → líquid transferència
expansion  = DT / T_fred_K   # ΔV/V gas ideal = 0.252 (25.2%)

print("=" * 70)
print("KILÒMETRE STIRLING — Flota de tubs xicotets")
print(f"ΔT = {DT}°C  |  η_Carnot = {eta_Carnot*100:.1f}%  |  η_real = {eta_real*100:.1f}%")
print(f"Expansió gas ideal: {expansion*100:.1f}%")
print("=" * 70)

# ── Funció: anàlisi d'un tub Stirling ────────────────────────
def tub_stirling(R, L_mult, rho_f, cost_m3_f,
                 Q_flux=10.0,   # W/m² flux tèrmic volcà
                 frac_gas=0.03, frac_liq=0.02,
                 eta=None):
    """
    Anàlisi complet d'un tub Stirling vertical.
    Retorna totes les mètriques per unitat.
    """
    if eta is None:
        eta = eta_real

    L         = R * L_mult
    V_tub     = np.pi * R**2 * L
    V_gas_ref = frac_gas * V_tub   # volum gas a T_fred
    V_liq     = frac_liq * V_tub
    V_obj_fred = V_gas_ref + V_liq  # volum total objecte fred
    V_obj_cal  = V_gas_ref*(1+expansion) + V_liq  # calent

    # Massa objecte (carcassa lleugera + líquid)
    rho_carcassa = 100   # kg/m³ equivalent carcassa lleugera
    m_obj = rho_carcassa * V_obj_fred + 800*V_liq  # oli+carcassa

    # Força neta en cada zona
    Fn_fred = rho_f * V_obj_fred * g - m_obj * g
    Fn_cal  = rho_f * V_obj_cal  * g - m_obj * g
    delta_Fn = Fn_cal - Fn_fred   # sempre positiu

    # Treball per cicle (tub vertical: objecte puja fred, baixa calent)
    W_cicle = delta_Fn * L   # J

    # Velocitat rotació (límit centrífug acer)
    om_max = np.sqrt(250e6 / (7800 * R**2))
    rpm    = om_max * 60 / (2*np.pi)
    freq   = om_max / (2*np.pi)   # Hz

    # Temps de trànsit (cal que sigui > temps d'intercanvi tèrmic)
    pas_sinfin = R * 0.1   # pas proporcional a R
    v_ax       = pas_sinfin * om_max
    t_transit  = (L/2) / v_ax   # s (mig recorregut)

    # Potència mecànica bruta
    P_mec = W_cicle * freq   # W

    # Superfície d'intercanvi tèrmic
    A_sup = 2 * np.pi * R * L   # m²

    # Calor disponible del volcà
    Q_disp = Q_flux * A_sup   # W (total disponible)

    # Potència útil: mínim de mecànica i calor disponible × η
    P_util = min(P_mec, Q_disp * eta)   # W

    # Cost per unitat
    rho_acer = 7800
    t_paret  = max(R/100, 0.002)  # tub xicotet → paret més fina
    m_cont   = rho_acer * 2*np.pi*R * L * t_paret
    m_fluid  = rho_f * V_tub
    cost_cont  = m_cont * 3
    cost_fluid = cost_m3_f * V_tub
    cost_obj   = m_obj * 20   # objecte 3 cambres, fabricació
    cost_unit  = cost_cont + cost_fluid + cost_obj + 500  # +500€ instal·lació

    # Watts per euro
    W_per_euro = P_util / cost_unit if cost_unit > 0 else 0

    return {
        'R': R, 'L': L, 'V_tub': V_tub,
        'delta_Fn': delta_Fn, 'W_cicle': W_cicle,
        'om_max': om_max, 'rpm': rpm, 'freq': freq,
        't_transit': t_transit,
        'P_mec': P_mec, 'Q_disp': Q_disp, 'P_util': P_util,
        'A_sup': A_sup,
        'cost_unit': cost_unit, 'W_per_euro': W_per_euro,
        'm_cont': m_cont, 'm_fluid': m_fluid,
        'V_gas_ref': V_gas_ref, 'V_obj_fred': V_obj_fred,
        'V_obj_cal': V_obj_cal, 'Fn_fred': Fn_fred, 'Fn_cal': Fn_cal,
    }

# ── Escombrat de mides (R=0.05m → 2m) ───────────────────────
print("\n=== ESCOMBRAT DE MIDES: R=5cm → 2m ===")
print(f"Fluid: Salmuera  |  Q_flux = 10 W/m²  |  η = {eta_real*100:.1f}%")
print(f"\n{'R(m)':>6} {'L(m)':>5} {'ΔFn(N)':>8} {'W_cicle(J)':>11} "
      f"{'rpm':>7} {'t_trans(s)':>11} {'P_util(W)':>10} "
      f"{'cost(€)':>8} {'W/€':>7}")
print("-" * 80)

Rs_sweep = [0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.0, 1.5, 2.0]
fl_sal   = fluids['Salmuera']
res_sweep = []
for R in Rs_sweep:
    r = tub_stirling(R, 4, fl_sal['rho'], fl_sal['cost_m3'])
    res_sweep.append(r)
    print(f"{R:>6.2f} {r['L']:>5.1f} {r['delta_Fn']:>8.3f} "
          f"{r['W_cicle']:>11.4f} {r['rpm']:>7.0f} "
          f"{r['t_transit']:>11.2f} {r['P_util']:>10.3f} "
          f"{r['cost_unit']:>8.0f} {r['W_per_euro']*1000:>7.3f}m")

# Troba l'òptim
idx_opt = np.argmax([r['W_per_euro'] for r in res_sweep])
r_opt   = res_sweep[idx_opt]
print(f"\n★ ÒPTIM: R={r_opt['R']}m  P={r_opt['P_util']:.3f}W  "
      f"W/€={r_opt['W_per_euro']*1000:.3f}mW/€")

# ── Quants tubs per MW? ──────────────────────────────────────
print(f"\n=== FLOTA: quants tubs per assolir 1MW? ===")
print(f"Q_flux = 10 W/m² (flanc volcànic actiu)")
print(f"\n{'R(m)':>6} {'P_unit(W)':>10} {'N per 1MW':>10} "
      f"{'A_total(m²)':>12} {'cost_total(k€)':>15} {'€/W':>7}")
print("-" * 65)

for R in [0.05, 0.10, 0.20, 0.50, 1.0, 2.0]:
    r = tub_stirling(R, 4, fl_sal['rho'], fl_sal['cost_m3'])
    if r['P_util'] > 0:
        N_1MW    = int(np.ceil(1e6 / r['P_util']))
        A_total  = N_1MW * r['A_sup']
        cost_tot = N_1MW * r['cost_unit'] / 1000
        euro_W   = N_1MW * r['cost_unit'] / 1e6
        print(f"{R:>6.2f} {r['P_util']:>10.3f} {N_1MW:>10,} "
              f"{A_total:>12,.0f} {cost_tot:>15,.0f} {euro_W:>7.2f}")

# ── Efecte Q_flux (densitat del volcà) ──────────────────────
print(f"\n=== EFECTE Q_flux: quin volcà necessitem? ===")
print(f"R=0.20m (òptim econòmic)")
print(f"\n{'Q_flux(W/m²)':>13} {'Font':>25} {'P_unit(W)':>10} "
      f"{'N per 1MW':>10} {'A_km²':>8}")
print("-" * 72)

fonts_calor = [
    (0.1,  'OTEC difús'),
    (1.0,  'Font hidrotermal difusa'),
    (5.0,  'Flanc volcànic suau'),
    (10.0, 'Flanc volcànic actiu ★'),
    (50.0, 'Font blanca (white smoker)'),
    (100.0,'Font negra (edge)'),
    (500.0,'Font negra (core)'),
]

R_fix = 0.20
for Q_f, nom_f in fonts_calor:
    r = tub_stirling(R_fix, 4, fl_sal['rho'], fl_sal['cost_m3'], Q_flux=Q_f)
    if r['P_util'] > 0.001:
        N_1MW   = int(np.ceil(1e6 / r['P_util']))
        A_km2   = N_1MW * r['A_sup'] / 1e6
        mark    = " ★" if Q_f == 10.0 else ""
        print(f"{Q_f:>13.1f} {nom_f:>25} {r['P_util']:>10.4f} "
              f"{N_1MW:>10,} {A_km2:>8.3f}{mark}")

# ── Comparació fluids (R=0.20m, Q=10W/m²) ───────────────────
print(f"\n=== COMPARACIÓ FLUIDS (R=0.20m, Q=10W/m²) ===")
print(f"{'Fluid':>14} {'P_util(W)':>10} {'cost(€)':>9} "
      f"{'W/€(mW/€)':>11} {'N per 1MW':>11} {'cost_1MW(k€)':>13}")
print("-" * 72)

for f_nom, fl in fluids.items():
    r = tub_stirling(0.20, 4, fl['rho'], fl['cost_m3'], Q_flux=10.0)
    N = int(np.ceil(1e6/r['P_util'])) if r['P_util']>0 else 999999
    cost_1MW = N * r['cost_unit'] / 1000
    print(f"{f_nom:>14} {r['P_util']:>10.4f} {r['cost_unit']:>9.0f} "
          f"{r['W_per_euro']*1000:>11.4f} {N:>11,} {cost_1MW:>13,.0f}")

# ── La clau: 120W per unitat ─────────────────────────────────
print(f"\n=== LA CLAU: ~120W PER UNITAT (R=0.50m, Q=50W/m²) ===")
r_120 = tub_stirling(0.50, 4, fl_sal['rho'], fl_sal['cost_m3'], Q_flux=50.0)
print(f"""
  R=0.50m, L=2m, Salmuera, Q=50W/m² (font blanca suau):
  P_util   = {r_120['P_util']:.1f} W
  A_sup    = {r_120['A_sup']:.2f} m²
  cost_unit= {r_120['cost_unit']:.0f} €
  W/€      = {r_120['W_per_euro']*1000:.2f} mW/€

  Per 1 MW:
    N = {int(np.ceil(1e6/r_120['P_util'])):,} tubs
    Cost = {int(np.ceil(1e6/r_120['P_util']))*r_120['cost_unit']/1e6:.1f} M€
    A_total = {int(np.ceil(1e6/r_120['P_util']))*r_120['A_sup']/1e4:.0f} hectàrees

  Per 100 MW (parc gran):
    N = {int(np.ceil(100e6/r_120['P_util'])):,} tubs
    Cost = {int(np.ceil(100e6/r_120['P_util']))*r_120['cost_unit']/1e6:.1f} M€

  Comparació amb eòlica offshore:
    100 MW eòlica: ~200M€  ({200e6/(100e6):.1f} €/W)
    100 MW Stirling: {int(np.ceil(100e6/r_120['P_util']))*r_120['cost_unit']/1e6:.0f}M€  ({int(np.ceil(100e6/r_120['P_util']))*r_120['cost_unit']/(100e6):.1f} €/W)
""")

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
fig.suptitle("KILÒMETRE STIRLING — Flota de tubs xicotets\n"
             "Màxim W/€: molts tubs, calor del volcà, fluids econòmics",
             fontsize=13, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.38)

# 1. P_util vs R
ax = fig.add_subplot(gs[0, 0])
Rs_p = np.logspace(-2, 1, 100)
for f_nom, fl in fluids.items():
    Ps = [tub_stirling(R, 4, fl['rho'], fl['cost_m3'],
                       Q_flux=10)['P_util'] for R in Rs_p]
    ax.loglog(Rs_p, Ps, lw=2, color=fl['color'], label=f_nom)
ax.axhline(120, color='red', lw=2, ls='--', label='120W referència')
ax.axvline(0.50, color='gray', lw=1, ls=':')
ax.set_xlabel('R [m]'); ax.set_ylabel('P_útil per tub [W]')
ax.set_title('Potència per tub vs R\n(Q_flux=10W/m²)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 2. W/€ vs R (la mètrica clau)
ax = fig.add_subplot(gs[0, 1])
for f_nom, fl in fluids.items():
    WE = [tub_stirling(R, 4, fl['rho'], fl['cost_m3'],
                       Q_flux=10)['W_per_euro']*1000 for R in Rs_p]
    ax.semilogx(Rs_p, WE, lw=2, color=fl['color'], label=f_nom)
ax.axvline(r_opt['R'], color='red', lw=2, ls='--',
           label=f'R_òptim={r_opt["R"]}m')
ax.set_xlabel('R [m]'); ax.set_ylabel('mW per € [mW/€]')
ax.set_title('W/€ vs R — MÈTRICA CLAU\n(on és millor el retorn?)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# 3. N tubs per 1MW vs Q_flux
ax = fig.add_subplot(gs[0, 2])
Qs_p = np.logspace(-1, 3, 100)
for R_fix2 in [0.10, 0.20, 0.50, 1.0]:
    Ns = []
    for Q in Qs_p:
        r = tub_stirling(R_fix2, 4, fl_sal['rho'],
                         fl_sal['cost_m3'], Q_flux=Q)
        Ns.append(int(np.ceil(1e6/r['P_util'])) if r['P_util']>0.001 else 1e7)
    ax.loglog(Qs_p, Ns, lw=2, label=f'R={R_fix2}m')
ax.axvline(10,  color='orange', lw=1, ls='--', label='Flanc volcànic 10W/m²')
ax.axvline(50,  color='red',    lw=1, ls='--', label='Font blanca 50W/m²')
ax.axhline(1000,color='gray',   lw=1, ls=':', label='1000 tubs')
ax.axhline(10000,color='gray',  lw=1, ls='--',label='10000 tubs')
ax.set_xlabel('Q_flux [W/m²]'); ax.set_ylabel('N tubs per 1 MW')
ax.set_title('Tubs necessaris per 1MW\nvs flux volcànic')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

# 4. Cost per W vs Q_flux
ax = fig.add_subplot(gs[1, 0])
for R_fix2 in [0.10, 0.20, 0.50, 1.0]:
    euros_W = []
    for Q in Qs_p:
        r = tub_stirling(R_fix2, 4, fl_sal['rho'],
                         fl_sal['cost_m3'], Q_flux=Q)
        eW = r['cost_unit'] / r['P_util'] if r['P_util']>0.001 else 1e6
        euros_W.append(eW)
    ax.loglog(Qs_p, euros_W, lw=2, label=f'R={R_fix2}m')
# Referències
ax.axhline(1.0, color='blue',  lw=2, ls='--', label='Li-ion 1€/W')
ax.axhline(1.5, color='green', lw=1, ls='--', label='Eòlica 1.5€/W')
ax.axhline(3.0, color='orange',lw=1, ls='--', label='Solar 3€/W')
ax.axvline(50,  color='red',   lw=1, ls=':')
ax.set_xlabel('Q_flux [W/m²]'); ax.set_ylabel('€/W')
ax.set_title('Cost per Watt vs flux volcànic\n(competitiu si < 1-3 €/W)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

# 5. Diagrama escalabilitat: P vs cost (flota)
ax = fig.add_subplot(gs[1, 1])
Ns_flota = np.logspace(1, 6, 100)
for R_fix2, Q_fix, ls in [(0.20, 10, '-'), (0.50, 50, '--'),
                            (0.20, 50, ':')]:
    r = tub_stirling(R_fix2, 4, fl_sal['rho'],
                     fl_sal['cost_m3'], Q_flux=Q_fix)
    Ps_fl   = Ns_flota * r['P_util'] / 1e6   # MW
    costs_fl= Ns_flota * r['cost_unit'] / 1e6 # M€
    ax.loglog(costs_fl, Ps_fl, lw=2, ls=ls,
              label=f'R={R_fix2}m Q={Q_fix}W/m²')
# Refs eòlica offshore
ax.loglog([200, 2000], [100, 1000], 'k--', lw=2,
          label='Eòlica offshore (ref)')
ax.set_xlabel('Cost [M€]'); ax.set_ylabel('Potència [MW]')
ax.set_title('Escalabilitat flota\n(dalt-esquerra = millor)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 6. Resum visual: el principi bàsic
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
ax.set_title('El principi dels tubs xicotets', fontsize=11,
             fontweight='bold')

text = f"""
PRINCIPI FONAMENTAL:

  P_total = N × P_unit
           = N × (Q_flux × A_sup × η)
           = Q_flux × (N×A_sup) × η
           = Q_flux × A_TOTAL × η

A_TOTAL és el que importa.
Molts tubs xicotets = gran A_TOTAL.

AVANTATGES dels tubs xicotets:
  ✓ Millor intercanvi tèrmic (Nusselt↑)
  ✓ Temps fabricació/instal·lació menor
  ✓ Redundància (si falla un, res passa)
  ✓ Menys materials per unitat
  ✓ Pressió mecànica menor
  ✓ Manteniment modular

FLUID ÒPTIM:  Salmuera 25%
  ρ = 1200 kg/m³ → Δρ suficient
  Cost ≈ 0.5 €/m³  (quasi gratis)
  No tòxic, no sedimenta, abundant

PER 1 MW (Q=10W/m², R=0.2m):
  N = {int(np.ceil(1e6/tub_stirling(0.2,4,1200,0.5,Q_flux=10)['P_util'])):,} tubs
  A_total = {int(np.ceil(1e6/tub_stirling(0.2,4,1200,0.5,Q_flux=10)['P_util']))*tub_stirling(0.2,4,1200,0.5,Q_flux=10)['A_sup']:.0f} m²
  Cost ~ {int(np.ceil(1e6/tub_stirling(0.2,4,1200,0.5,Q_flux=10)['P_util']))*tub_stirling(0.2,4,1200,0.5,Q_flux=10)['cost_unit']/1e6:.1f} M€

La 120W és una conseqüència del flux,
no un límit intrínsec del disseny.
"""
ax.text(0.05, 0.95, text, transform=ax.transAxes,
        fontsize=8.5, va='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.savefig('/mnt/user-data/outputs/kilometre_stirling_flota_xicotets.png',
            dpi=150, bbox_inches='tight')
print("Fitxers guardats: kilometre_stirling_flota_xicotets.py  |  .png")
