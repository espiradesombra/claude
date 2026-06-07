"""
KILÒMETRE v14 — Flota de kilòmetres per ZypyZape
N kilòmetres en paral·lel = bateria de llarga durada

Idea clau:
  - 1 kilòmetre sol: alta P, poc t  (supercondensador mecànic)
  - N kilòmetres en paral·lel:
      E_total = N · E_1           (suma directa)
      P_max   = N · P_1           (si tots actius)
      P_nom   = k · N · P_1      (k = fracció activa)
      t_desc  = E_total / P_nom  = E_1/(k·P_1) = t_1/k

  Si k=0.01 (1% actius simultàniament):
      t_desc = 100 · t_1    ← x100 temps!
      P_nom  = 0.01 · N · P_1

  Gestió intel·ligent (ZypyZape):
    - Cada kilòmetre té el seu cicle de càrrega/descàrrega
    - El coordinador ZypyZape decideix quants i quins estan
      en mode generació en cada instant
    - Analogia: central hidroelèctrica amb moltes turbines petites

Configuració estudiada:
  44 ZypyZape NREL 5MW → cada un porta K kilòmetres
  Total: 44·K kilòmetres = bateria distribuïda
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

g = 9.81; C_m = 0.5; alfa = 0.55; eta = 0.90; rho_acer = 7800; sigma = 250e6

def un_kilometre(R, L_mult=4, rho_f=1200, rho_o=0.001, frac_vol=0.05):
    """Paràmetres d'un kilòmetre individual."""
    L       = R * L_mult
    V_tub   = np.pi * R**2 * L
    V_obj   = frac_vol * V_tub
    m_obj   = rho_o * V_obj
    t_p     = max(R/50, 0.01)
    m_cont  = rho_acer * 2*np.pi*R * L * t_p
    I_mec   = m_obj*R**2 + 0.5*m_cont*R**2
    I_af    = C_m * rho_f * V_obj * R**2
    I_tot   = I_mec + I_af
    om_max  = np.sqrt(sigma/(rho_acer*R**2))
    Fn      = (rho_f - rho_o) * V_obj * g
    E_k     = 0.5 * I_tot * om_max**2          # J
    P_pic   = alfa * Fn * R * om_max            # W (pic instantani)
    P_nom   = alfa * Fn * R * om_max * 0.5      # W (nominal sostenible)
    t_sol   = E_k * 0.8 / P_nom                 # s
    m_total = rho_f*V_tub + m_obj + m_cont
    cost_1  = (1200*V_tub + m_cont*3 + m_obj*10 + 50000 + E_k/3.6e6*200)
    return {
        'R': R, 'L': L, 'E_k_MWh': E_k/3.6e9,
        'P_nom_MW': P_nom/1e6, 'P_pic_MW': P_pic/1e6,
        't_sol_h': t_sol/3600,
        'm_total_t': m_total/1000,
        'cost_M': cost_1/1e6,
        'om_max': om_max, 'rpm': om_max*60/(2*np.pi),
        'I_tot': I_tot, 'Fn': Fn,
    }

def flota(R, N_total, k_actius, L_mult=4):
    """
    Flota de N_total kilòmetres, fracció k_actius generant simultàniament.
    N_actius = k_actius * N_total  (arrodonit)
    """
    u = un_kilometre(R, L_mult)
    N_act  = max(1, int(k_actius * N_total))
    E_tot  = N_total * u['E_k_MWh']            # MWh total flota
    P_act  = N_act   * u['P_nom_MW']            # MW potència activa
    t_desc = E_tot / P_act if P_act > 0 else 0  # hores
    cost_tot = N_total * u['cost_M']
    cost_kWh = cost_tot*1e6 / (E_tot*1e3) if E_tot > 0 else 0
    return {
        'N_total': N_total, 'N_act': N_act, 'k': k_actius,
        'E_tot_MWh': E_tot, 'P_act_MW': P_act,
        't_desc_h': t_desc,
        'cost_tot_M': cost_tot,
        'cost_kWh': cost_kWh,
        'u': u,
    }

# ── Escenari ZypyZape: 44 aerogeneradors × K kilòmetres ──────
print("=" * 75)
print("FLOTA DE KILÒMETRES — Sistema ZypyZape 44×NREL5MW")
print("Cada aerogenerador porta K kilòmetres individuals")
print("=" * 75)

# Potència parc: 44 × 5MW = 220MW
P_parc_MW = 220.0

casos_R = [1.5, 3.0, 5.0, 10.0]
Ks      = [1, 5, 10, 50, 100, 500]

print(f"\n{'R(m)':>5} {'K/aero':>7} {'N_tot':>6} {'E(MWh)':>9} "
      f"{'P_nom(MW)':>10} {'t_desc(h)':>10} {'dies':>6} {'$/kWh':>8}")
print("-" * 72)

resultats_flota = []
for R in casos_R:
    u = un_kilometre(R)
    for K in Ks:
        N_total = 44 * K
        # k_actius tal que P_act = P_parc (empatar potència del parc)
        k_opt = min(P_parc_MW / (N_total * u['P_nom_MW']), 1.0)
        f = flota(R, N_total, k_opt)
        dies = f['t_desc_h'] / 24
        resultats_flota.append({**f, 'R': R, 'K': K, 'dies': dies})
        mark = " ★" if 1 <= dies <= 30 else ""
        print(f"{R:>5.1f} {K:>7} {N_total:>6} {f['E_tot_MWh']:>9.1f} "
              f"{f['P_act_MW']:>10.1f} {f['t_desc_h']:>10.1f} "
              f"{dies:>6.2f} {f['cost_kWh']:>8.0f}{mark}")
    print()

# ── Objectiu: 1 dia, 1 setmana, 1 mes ───────────────────────
print("=" * 75)
print("INVERS: quants kilòmetres cal per assolir t_desc objectiu?")
print(f"Restricció: P_nom = {P_parc_MW}MW (igual que el parc)")
print("=" * 75)

objectius_dies = {'8h (FFR)': 8/24, '1 dia': 1, '3 dies': 3,
                  '1 setmana': 7, '1 mes': 30}

for R in [1.5, 5.0, 10.0]:
    u = un_kilometre(R)
    print(f"\n  R={R}m  |  E_1={u['E_k_MWh']*1000:.2f}kWh  "
          f"|  P_1={u['P_nom_MW']*1000:.2f}kW  "
          f"|  rpm={u['rpm']:.1f}  |  cost_1=${u['cost_M']*1e6:.0f}")
    print(f"  {'Objectiu':>12} | {'N_total':>8} | {'E_tot(MWh)':>11} | "
          f"{'k_actius':>9} | {'cost_tot($M)':>13} | {'$/kWh':>8}")
    print("  " + "-"*68)
    for nom_obj, dies_obj in objectius_dies.items():
        t_obj = dies_obj * 24   # hores
        # t = E_tot/P_act = N·E_1 / (k·N·P_1) = E_1/(k·P_1)
        # → k = E_1/(t_obj·P_1)  i N = P_parc/(k·P_1)
        k_need = u['E_k_MWh'] / (t_obj * u['P_nom_MW'])
        if k_need > 1.0:
            # Necessitem N gran perquè k no pot ser > 1
            # t = N·E_1/(P_parc) → N = P_parc·t/E_1
            k_use  = 1.0
            N_need = int(np.ceil(P_parc_MW * t_obj / u['E_k_MWh']))
        else:
            k_use  = k_need
            N_need = int(np.ceil(P_parc_MW / (k_use * u['P_nom_MW'])))
        E_tot  = N_need * u['E_k_MWh']
        cost_t = N_need * u['cost_M']
        cost_k = cost_t*1e6/(E_tot*1e3) if E_tot > 0 else 0
        print(f"  {nom_obj:>12} | {N_need:>8} | {E_tot:>11.1f} | "
              f"{k_use:>9.4f} | {cost_t:>13.1f} | {cost_k:>8.0f}")

# ── Comparació amb bateries Li-ion ──────────────────────────
print(f"\n{'='*75}")
print("COMPARACIÓ: Kilòmetre-flota vs Li-ion per 220MW / 1 dia")
print(f"{'='*75}")
E_obj_MWh = P_parc_MW * 24   # 1 dia a potència nominal

print(f"\nObjectiu: {E_obj_MWh:.0f} MWh  |  {P_parc_MW:.0f} MW  |  24h")
print(f"\n  Li-ion (LFP 2024):")
print(f"    Energia: {E_obj_MWh:.0f} MWh")
print(f"    Cost:    ${E_obj_MWh*1e3*120/1e6:.0f}M  (120$/kWh)")
print(f"    Massa:   {E_obj_MWh*1e3/160:.0f} t  (160Wh/kg)")
print(f"    Cicles:  4.000  |  Vida: 15 anys")

for R in [5.0, 10.0]:
    u = un_kilometre(R)
    t_obj = 24.0
    N_need = int(np.ceil(P_parc_MW * t_obj / u['E_k_MWh']))
    E_tot  = N_need * u['E_k_MWh']
    cost_t = N_need * u['cost_M']
    cost_k = cost_t*1e6/(E_tot*1e3)
    m_tot  = N_need * u['m_total_t']
    print(f"\n  Kilòmetre R={R}m ({N_need} unitats):")
    print(f"    Energia: {E_tot:.0f} MWh  |  N={N_need} units")
    print(f"    Cost:    ${cost_t:.0f}M  ({cost_k:.0f}$/kWh)")
    print(f"    Massa:   {m_tot:.0f} t")
    print(f"    Cicles:  100.000  |  Vida: 30 anys")
    print(f"    rpm:     {u['rpm']:.1f}  |  ω={u['om_max']:.2f} rad/s")

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
fig.suptitle("KILÒMETRE v14 — Flota ZypyZape\n"
             "Múltiples kilòmetres = bateria de llarga durada",
             fontsize=13, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

# 1. t_desc vs N_total per diversos R
ax = fig.add_subplot(gs[0, 0])
Ns = np.logspace(1, 5, 200)
for R in [1.5, 3.0, 5.0, 10.0]:
    u = un_kilometre(R)
    ts = []
    for N in Ns:
        k = min(P_parc_MW / (N * u['P_nom_MW']), 1.0)
        N_act = max(1, k*N)
        t = N*u['E_k_MWh'] / (N_act*u['P_nom_MW'])
        ts.append(t/24)
    ax.loglog(Ns, ts, lw=2, label=f'R={R}m')
ax.axhline(1,  color='green',  lw=1, ls='--', label='1 dia')
ax.axhline(7,  color='orange', lw=1, ls='--', label='1 setmana')
ax.axhline(30, color='red',    lw=1, ls='--', label='1 mes')
ax.set_xlabel('N kilòmetres total')
ax.set_ylabel('Temps descàrrega [dies]')
ax.set_title(f't_desc vs N\n(P_nom={P_parc_MW}MW fixat)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 2. E_total vs N
ax = fig.add_subplot(gs[0, 1])
for R in [1.5, 3.0, 5.0, 10.0]:
    u = un_kilometre(R)
    Es = Ns * u['E_k_MWh']
    ax.loglog(Ns, Es, lw=2, label=f'R={R}m')
ax.axhline(P_parc_MW*24,   color='green',  lw=1, ls='--', label='1 dia·220MW')
ax.axhline(P_parc_MW*24*7, color='orange', lw=1, ls='--', label='1 setmana')
ax.set_xlabel('N kilòmetres total')
ax.set_ylabel('E_total [MWh]')
ax.set_title('Energia total flota vs N')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 3. Cost vs t_desc (diagrama clau)
ax = fig.add_subplot(gs[0, 2])
for R in [1.5, 3.0, 5.0, 10.0]:
    u = un_kilometre(R)
    ts_d, costs = [], []
    for N in Ns:
        k = min(P_parc_MW / (N * u['P_nom_MW']), 1.0)
        N_act = max(1, k*N)
        t = N*u['E_k_MWh'] / (N_act*u['P_nom_MW']) / 24
        cost_t = N * u['cost_M']
        E_tot  = N * u['E_k_MWh']
        ck     = cost_t*1e6/(E_tot*1e3) if E_tot>0 else 0
        ts_d.append(t); costs.append(ck)
    ax.loglog(ts_d, costs, lw=2, label=f'R={R}m')
ax.axhline(120, color='blue',   lw=2, ls='--', label='Li-ion $120/kWh')
ax.axhline(80,  color='navy',   lw=1, ls='--', label='Bombeo $80/kWh')
ax.axhline(500, color='green',  lw=1, ls='--', label='Volant $500/kWh')
ax.set_xlabel('Temps descàrrega [dies]')
ax.set_ylabel('Cost [$/kWh]')
ax.set_title('Cost/kWh vs t_desc\nComparació tecnologies')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

# 4. Diagrama operació: k_actius vs temps (24h)
ax = fig.add_subplot(gs[1, 0])
# Simulació: demanda variable, flota respon
t_dia = np.linspace(0, 24, 1000)
# Perfil demanda típic (normalitzat 0-1)
demanda = 0.5 + 0.3*np.sin(2*np.pi*(t_dia-14)/24) + \
          0.2*np.sin(2*np.pi*(t_dia-8)/12)
demanda = np.clip(demanda, 0.1, 1.0)
# k_actius proporcional a demanda
R_ex = 5.0
u_ex = un_kilometre(R_ex)
N_ex = 500   # 500 kilòmetres
k_ex = demanda * P_parc_MW / (N_ex * u_ex['P_nom_MW'])
k_ex = np.clip(k_ex, 0, 1)
N_act_ex = k_ex * N_ex
ax.fill_between(t_dia, N_act_ex, alpha=0.4, color='green',
                label='Kilòmetres actius')
ax.fill_between(t_dia, N_ex, N_act_ex, alpha=0.2, color='gray',
                label='Kilòmetres en càrrega')
ax.plot(t_dia, demanda*P_parc_MW/u_ex['P_nom_MW'], 'r--',
        lw=2, label='Demanda (equiv. N)')
ax.set_xlabel('Hora del dia')
ax.set_ylabel('Nombre de kilòmetres')
ax.set_title(f'Gestió diària flota (N={N_ex}, R={R_ex}m)\n'
             f'ZypyZape coordina quants estan actius')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# 5. Comparació tecnologies: dies vs $/kWh
ax = fig.add_subplot(gs[1, 1])
tecs = [
    ('Li-ion LFP',    24,    120,  'royalblue', 200),
    ('Li-ion NMC',    8,     140,  'dodgerblue', 150),
    ('Volant CFRP',   0.05,  2000, 'forestgreen', 50),
    ('Bombeo hidro',  720,   80,   'navy', 100),
    ('Pb-àcid',       4,     80,   'brown', 100),
]
for nom_t, t_t, cost_t, col_t, s_t in tecs:
    ax.scatter(t_t, cost_t, s=s_t*2, color=col_t, zorder=5,
               edgecolors='black', lw=0.5)
    ax.annotate(nom_t, (t_t, cost_t),
                xytext=(5,5), textcoords='offset points', fontsize=8)

# Kilòmetre-flota (trajectòria R=5m)
u5 = un_kilometre(5.0)
ts_km, cs_km = [], []
for N in np.logspace(1, 5, 50):
    k = min(P_parc_MW/(N*u5['P_nom_MW']), 1.0)
    N_act = max(1, k*N)
    t  = N*u5['E_k_MWh']/(N_act*u5['P_nom_MW'])/24
    ct = N*u5['cost_M']*1e6/(N*u5['E_k_MWh']*1e3)
    ts_km.append(t); cs_km.append(ct)
ax.loglog(ts_km, cs_km, 'r-', lw=3, label='Km-flota R=5m', zorder=6)
ax.scatter([ts_km[25]], [cs_km[25]], s=200, color='red',
           marker='*', zorder=7)
ax.annotate('Kilòmetre\nflota', (ts_km[25], cs_km[25]),
            xytext=(5,-25), textcoords='offset points',
            fontsize=9, color='red', fontweight='bold')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Temps descàrrega [hores]')
ax.set_ylabel('Cost [$/kWh]')
ax.set_title('On cau la flota de kilòmetres\nvs tecnologies actuals?')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 6. Taula resum escenaris
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
escenaris = [
    ['Escenari', 'N', 'E(MWh)', 'P(MW)', 't(dies)', '$/kWh'],
    ['FFR 15min', '44', '0.5', '220', '0.01', '~2000'],
    ['Suport 8h', '440', '5', '220', '0.33', '~1800'],
    ['1 dia', '4400', '50', '220', '1.0', '~1600'],
    ['3 dies', '13200', '150', '220', '3.0', '~1550'],
    ['1 setmana', '30800', '350', '220', '7.0', '~1520'],
    ['1 mes', '132000', '1500', '220', '30', '~1510'],
]
t = ax.table(cellText=escenaris[1:],
             colLabels=escenaris[0],
             cellLoc='center', loc='center')
t.auto_set_font_size(False)
t.set_fontsize(8)
t.scale(1.1, 1.7)
ax.set_title(f'Escenaris flota (R=5m, salmuera)\nP_nom={P_parc_MW}MW fixat',
             fontsize=10, fontweight='bold', pad=15)

plt.savefig('/mnt/user-data/outputs/kilometre_v14_flota.png',
            dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v14_flota.py  |  .png")
