"""
KILÒMETRE v11 — Escombrat complet de dimensions
Troba el punt òptim de fabricació: on val la pena construir-lo?

Variables:
  R    : 0.1m → 10m  (radi tub)
  ω    : limitat per tensió centrífuga del contenidor
  fluid: aigua, salmuera, gallium
  V_obj: escala amb R³ (objecte proporcional al tub)

Mètriques d'avaluació:
  1. E_k [kWh]          → energia emmagatzemada
  2. P_pic [kW]         → potència de pic
  3. E_k/massa [Wh/kg]  → densitat energètica
  4. P/E [kW/kWh]       → ràtio potència/energia (C-rate)
  5. Cost relatiu       → proporcional a massa fluid + contenidor
  6. Viabilitat mecànica → tensió centrífuga contenidor

Llei d'escala clau:
  E_k = ½ · I_ef · ω²
  I_ef = R² · (m_obj + C_m·ρ_f·V_obj)
  Si V_obj ∝ R³:  I_ef ∝ ρ_f · R⁵
  Si ω ∝ 1/R (límit centrífug): E_k ∝ ρ_f · R³  (escala volumètrica!)
  → La bateria escala com el VOLUM del fluid, no com la massa
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm

g    = 9.81
C_m  = 0.5
eta_gen = 0.90
eta_mot = 0.92
alfa_gen = 0.55

# ── Fluids ───────────────────────────────────────────────────
fluids = {
    'Agua':         {'rho': 1000,  'mu': 1.0e-3, 'color': 'steelblue',
                     'sigma_cont': 250e6, 'rho_cont': 7800,
                     'cost_rel': 1,   'nom': 'H₂O'},
    'Salmuera 25%': {'rho': 1200,  'mu': 1.5e-3, 'color': 'cyan',
                     'sigma_cont': 250e6, 'rho_cont': 7800,
                     'cost_rel': 1.1, 'nom': 'Sal'},
    'Gallium':      {'rho': 6095,  'mu': 2.0e-3, 'color': 'purple',
                     'sigma_cont': 400e6, 'rho_cont': 7800,
                     'cost_rel': 50,  'nom': 'Ga'},
}

# Objecte: escuma EPS (ρ=15 kg/m³), volum = 5% del tub
rho_obj   = 15.0    # kg/m³
frac_vol  = 0.05    # V_obj = frac_vol * V_tub

# ── Tensió centrífuga contenidor cilíndric ───────────────────
def omega_max_contenidor(R, sigma_max, rho_cont):
    """ω màxim per tensió hoop: σ = ρ·ω²·R²"""
    return np.sqrt(sigma_max / (rho_cont * R**2))

# ── Escombrat de R ───────────────────────────────────────────
Rs = np.logspace(-1, 1, 200)   # 0.1m → 10m

print("=" * 75)
print("KILÒMETRE v11 — Escombrat de dimensions R: 0.1m → 10m")
print("Objecte EPS (ρ=15 kg/m³), V_obj = 5% del volum del tub")
print("=" * 75)

fig = plt.figure(figsize=(16, 18))
fig.suptitle("KILÒMETRE v11 — Escombrat complet de dimensions\n"
             "On val la pena fabricar la bateria cinètica?",
             fontsize=13, fontweight='bold')
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)

# Precomputa per cada fluid
resultats = {}

for f_nom, fl in fluids.items():
    rho_f    = fl['rho']
    mu_f     = fl['mu']
    sigma_c  = fl['sigma_cont']
    rho_c    = fl['rho_cont']

    # Dimensions que escalen amb R
    L_tub    = Rs * 4          # longitud tub = 4×R (proporcional)
    V_tub    = np.pi * Rs**2 * L_tub
    V_obj    = frac_vol * V_tub
    m_obj    = rho_obj * V_obj
    m_fluid  = rho_f * V_tub
    m_cont   = rho_c * 2*np.pi*Rs * L_tub * (Rs/50)  # gruix paret ~ R/50

    # Inèrcia efectiva
    I_mec    = m_obj * Rs**2
    I_af     = C_m * rho_f * V_obj * Rs**2
    I_ef     = I_mec + I_af   # ∝ rho_f · R^5

    # ω màxim per contenidor
    om_max   = omega_max_contenidor(Rs, sigma_c, rho_c)
    rpm_max  = om_max * 60 / (2*np.pi)

    # Energia màxima
    E_k_max  = 0.5 * I_ef * om_max**2   # J
    E_k_kWh  = E_k_max / 3.6e6          # kWh

    # Força neta (flotació dominant)
    Fn       = (rho_f - rho_obj) * V_obj * g

    # Potència pic generador
    P_pic    = alfa_gen * Fn * Rs * om_max   # W

    # Fricció (Stokes, esfera equiv.)
    r_eq     = (3*V_obj/(4*np.pi))**(1/3)
    v_tip    = om_max * Rs
    F_drag   = 6 * np.pi * mu_f * r_eq * v_tip
    beta     = F_drag * Rs / om_max
    W_fric   = beta * om_max * 2*np.pi   # J/volta

    # Balanç net per volta
    W_gen_mec = alfa_gen * abs(Fn) * Rs * 2.0
    W_mot_mec = W_fric + W_gen_mec
    W_el_in   = W_mot_mec / eta_mot
    W_el_out  = W_gen_mec * eta_gen
    net_volta = W_el_out - W_el_in   # sempre negatiu

    # Eficiència round-trip
    eta_rt    = W_el_out / W_el_in * 100

    # Densitat energètica (respecte massa total sistema)
    m_total   = m_fluid + m_obj + m_cont
    Wh_kg     = E_k_kWh * 1000 / m_total   # Wh/kg

    # C-rate (P/E)
    C_rate    = P_pic / (E_k_max / 3600)   # 1/h → potència relativa

    # Cost relatiu (proporcional a massa fluid × cost_rel + massa cont)
    cost_rel  = fl['cost_rel'] * m_fluid + m_cont

    resultats[f_nom] = {
        'Rs': Rs, 'L_tub': L_tub, 'V_tub': V_tub,
        'I_ef': I_ef, 'om_max': om_max, 'rpm_max': rpm_max,
        'E_k_kWh': E_k_kWh, 'P_pic': P_pic,
        'Wh_kg': Wh_kg, 'C_rate': C_rate,
        'net_volta': net_volta, 'eta_rt': eta_rt,
        'm_total': m_total, 'cost_rel': cost_rel,
        'W_fric': W_fric, 'W_gen_mec': W_gen_mec,
        'color': fl['color'], 'nom': fl['nom'],
    }

# ── Taula de punts clau ──────────────────────────────────────
Rs_taula = [0.5, 1.0, 2.0, 5.0, 10.0]
print(f"\n{'R(m)':>6} | {'Fluid':>12} | {'ω_max':>7} | {'rpm':>6} | "
      f"{'E_k(kWh)':>9} | {'P_pic(kW)':>10} | {'Wh/kg':>7} | {'C-rate':>7}")
print("-" * 80)

for R_t in Rs_taula:
    idx = np.argmin(np.abs(Rs - R_t))
    for f_nom, r in resultats.items():
        print(f"{R_t:>6.1f} | {f_nom:>12} | {r['om_max'][idx]:>7.1f} | "
              f"{r['rpm_max'][idx]:>6.0f} | {r['E_k_kWh'][idx]:>9.4f} | "
              f"{r['P_pic'][idx]/1000:>10.2f} | {r['Wh_kg'][idx]:>7.4f} | "
              f"{r['C_rate'][idx]:>7.1f}")
    print()

# ── Lleis d'escala ───────────────────────────────────────────
print("=" * 60)
print("LLEIS D'ESCALA (si ω ∝ 1/R, límit centrífug)")
print("=" * 60)
print("""
  I_ef   ∝ ρ_f · R⁵        (inèrcia escala molt fort amb R)
  ω_max  ∝ 1/R              (velocitat baixa amb R)
  E_k    = ½·I_ef·ω²
         ∝ ρ_f · R⁵ · R⁻²
         ∝ ρ_f · R³         ← escala VOLUMÈTRICA (com el pes del fluid!)

  P_pic  = alfa·Fn·R·ω
         ∝ ρ_f · R³ · R · R⁻¹
         ∝ ρ_f · R³         ← escala igual que E_k

  C-rate = P/E ≈ constant   ← independent de R!

  Wh/kg  ∝ E_k / m_fluid
         ∝ ρ_f·R³ / ρ_f·R³
         ∝ constant         ← densitat energètica independent de R!

CONCLUSIÓ: el sistema escala perfectament.
Doblar R → 8× l'energia i 8× la potència.
La densitat (Wh/kg) i el C-rate no canvien amb R.
→ El punt òptim NO és una mida concreta sinó una qüestió de COST i FABRICACIÓ.
""")

# Punt on P_pic > 1 MW (interessant per FFR)
print("Punt on P_pic > 1 MW (interessant per regulació de xarxa):")
for f_nom, r in resultats.items():
    idx_mw = np.where(r['P_pic'] > 1e6)[0]
    if len(idx_mw) > 0:
        R_mw = Rs[idx_mw[0]]
        print(f"  {f_nom:>14}: R > {R_mw:.2f}m  "
              f"(ω={r['om_max'][idx_mw[0]]:.1f} rad/s, "
              f"{r['rpm_max'][idx_mw[0]]:.0f} rpm, "
              f"E={r['E_k_kWh'][idx_mw[0]]:.2f} kWh)")

print()
print("Punt on E_k > 1 MWh (emmagatzematge estacional):")
for f_nom, r in resultats.items():
    idx_mwh = np.where(r['E_k_kWh'] > 1000)[0]
    if len(idx_mwh) > 0:
        R_mwh = Rs[idx_mwh[0]]
        print(f"  {f_nom:>14}: R > {R_mwh:.2f}m  "
              f"(ω={r['om_max'][idx_mwh[0]]:.1f} rad/s, "
              f"{r['rpm_max'][idx_mwh[0]]:.0f} rpm)")
    else:
        print(f"  {f_nom:>14}: no arriba a 1 MWh amb R≤10m")

# ── Gràfics ──────────────────────────────────────────────────

# 1. Energia vs R (log-log)
ax = fig.add_subplot(gs[0, 0])
for f_nom, r in resultats.items():
    ax.loglog(r['Rs'], r['E_k_kWh']*1000, lw=2,
              color=r['color'], label=r['nom'])
ax.axhline(1, color='gray', lw=1, ls=':', label='1 Wh')
ax.axhline(1000, color='gray', lw=1, ls='--', label='1 kWh')
ax.axhline(1e6, color='gray', lw=1, ls='-', label='1 MWh')
ax.axvline(1.0, color='orange', lw=1, ls='--', alpha=0.5)
ax.axvline(10.0, color='red', lw=1, ls='--', alpha=0.5)
for R_t in [0.5, 1, 2, 5, 10]:
    ax.axvline(R_t, color='lightgray', lw=0.5)
ax.set_xlabel('R [m]'); ax.set_ylabel('E_k [Wh]')
ax.set_title('Energia emmagatzemada vs R\n(límit centrífug acer)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')
# Anotació llei d'escala
ax.text(0.3, 0.15, 'E ∝ R³', transform=ax.transAxes,
        fontsize=11, color='black', fontweight='bold')

# 2. Potència pic vs R
ax = fig.add_subplot(gs[0, 1])
for f_nom, r in resultats.items():
    ax.loglog(r['Rs'], r['P_pic']/1000, lw=2,
              color=r['color'], label=r['nom'])
ax.axhline(1, color='gray', lw=1, ls=':', label='1 kW')
ax.axhline(1000, color='gray', lw=1, ls='--', label='1 MW')
ax.axhline(1e6, color='gray', lw=1, ls='-', label='1 GW')
ax.set_xlabel('R [m]'); ax.set_ylabel('P_pic [kW]')
ax.set_title('Potència pic vs R\n(P ∝ R³)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')
ax.text(0.3, 0.15, 'P ∝ R³', transform=ax.transAxes,
        fontsize=11, color='black', fontweight='bold')

# 3. rpm vs R
ax = fig.add_subplot(gs[0, 2])
for f_nom, r in resultats.items():
    ax.loglog(r['Rs'], r['rpm_max'], lw=2,
              color=r['color'], label=r['nom'])
ax.axhline(3000, color='green', lw=1, ls='--', label='3000 rpm (industrial)')
ax.axhline(300,  color='orange', lw=1, ls='--', label='300 rpm (lent)')
ax.axhline(30,   color='red', lw=1, ls='--', label='30 rpm (molt lent)')
ax.set_xlabel('R [m]'); ax.set_ylabel('rpm màximes')
ax.set_title('Velocitat màxima vs R\n(ω_max ∝ 1/R)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')
ax.text(0.3, 0.15, 'ω ∝ 1/R', transform=ax.transAxes,
        fontsize=11, color='black', fontweight='bold')

# 4. Densitat energètica Wh/kg vs R
ax = fig.add_subplot(gs[1, 0])
for f_nom, r in resultats.items():
    ax.semilogx(r['Rs'], r['Wh_kg'], lw=2,
                color=r['color'], label=r['nom'])
ax.axhline(200, color='blue', lw=1, ls='--', label='Li-ion ~200Wh/kg')
ax.axhline(50,  color='green', lw=1, ls='--', label='Volant CFRP ~50Wh/kg')
ax.set_xlabel('R [m]'); ax.set_ylabel('Wh/kg')
ax.set_title('Densitat energètica vs R\n(constant! independent de R)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')
ax.text(0.05, 0.85, 'Wh/kg ≈ const', transform=ax.transAxes,
        fontsize=11, fontweight='bold', color='darkred')

# 5. C-rate vs R
ax = fig.add_subplot(gs[1, 1])
for f_nom, r in resultats.items():
    ax.semilogx(r['Rs'], r['C_rate'], lw=2,
                color=r['color'], label=r['nom'])
ax.axhline(1,  color='blue',  lw=1, ls='--', label='1C (Li-ion típic)')
ax.axhline(10, color='green', lw=1, ls='--', label='10C (alta potència)')
ax.set_xlabel('R [m]'); ax.set_ylabel('C-rate [1/h]')
ax.set_title('C-rate (P/E) vs R\n(constant! no depèn de R)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')
ax.text(0.05, 0.85, 'C-rate ≈ const', transform=ax.transAxes,
        fontsize=11, fontweight='bold', color='darkred')

# 6. Mapa de zona òptima (E_k vs P_pic, R com a paràmetre)
ax = fig.add_subplot(gs[1, 2])
for f_nom, r in resultats.items():
    sc = ax.scatter(r['E_k_kWh'], r['P_pic']/1000,
                    c=np.log10(r['Rs']), cmap='plasma',
                    s=15, alpha=0.8, label=r['nom'])
# Anotació mides clau
for R_t in [0.5, 1, 2, 5, 10]:
    idx = np.argmin(np.abs(Rs - R_t))
    r0 = list(resultats.values())[1]   # salmuera
    ax.annotate(f'R={R_t}m',
                (r0['E_k_kWh'][idx], r0['P_pic'][idx]/1000),
                fontsize=7, ha='left')
plt.colorbar(sc, ax=ax, label='log₁₀(R) [m]')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('E_k [kWh]'); ax.set_ylabel('P_pic [kW]')
ax.set_title('Mapa E_k vs P_pic\n(color = radi)')
ax.grid(True, alpha=0.3, which='both')

# 7. Quadrant de viabilitat: zones per aplicació
ax = fig.add_subplot(gs[2, :])
# Regions d'aplicació
apps = [
    ('UPS / micro-xarxa',     (0.001, 1),    (1, 100),     'lightyellow'),
    ('Regulació freq. (FFR)', (0.1, 100),    (100, 10000), 'lightgreen'),
    ('Suport de xarxa',       (10, 1000),    (1000,100000),'lightblue'),
    ('Emmagatzematge bulk',   (100, 100000), (10, 1000),   'lightsalmon'),
]
for nom_app, E_r, P_r, color_app in apps:
    ax.fill_between(E_r, P_r[0], P_r[1],
                    alpha=0.4, color=color_app, label=nom_app)
    ax.text(np.sqrt(E_r[0]*E_r[1]), np.sqrt(P_r[0]*P_r[1]),
            nom_app, ha='center', va='center', fontsize=9,
            fontweight='bold', color='darkgray')

# Trajectòria del kilòmetre per cada fluid
for f_nom, r in resultats.items():
    ax.loglog(r['E_k_kWh'], r['P_pic']/1000, lw=3,
              color=r['color'], label=f'Km {r["nom"]}',
              alpha=0.9)
    # Marca mides clau
    for R_t in [1, 2, 5, 10]:
        idx = np.argmin(np.abs(Rs - R_t))
        ax.plot(r['E_k_kWh'][idx], r['P_pic'][idx]/1000,
                'o', color=r['color'], ms=8)
        if f_nom == 'Salmuera 25%':
            ax.annotate(f'R={R_t}m',
                        (r['E_k_kWh'][idx], r['P_pic'][idx]/1000),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Energia emmagatzemada [kWh]', fontsize=11)
ax.set_ylabel('Potència pic [kW]', fontsize=11)
ax.set_title('Diagrama de Ragone — Zones d\'aplicació vs dimensions\n'
             'On encaixa el kilòmetre per cada mida?', fontsize=11)
ax.legend(fontsize=8, loc='lower right', ncol=2)
ax.grid(True, alpha=0.3, which='both')

plt.savefig('/mnt/user-data/outputs/kilometre_v11_escombrat.png',
            dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v11_escombrat_dimensions.py  |  .png")
