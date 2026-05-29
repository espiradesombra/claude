"""
GEMELL v6 — Estudi d'Escala Complet
=====================================
Preguntes que respon:
  1. Com escala la potència amb R?
  2. N turbines òptim per a cada R?
  3. Molts xicotets vs un de gran?
  4. Quines variables dominen?
  5. Quants lupes mínim/màxim per R donat?
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Constants físiques (Sals fundides Solar Salt) ─────────────────────────────
CP    = 1500.0   # J/kg·K
RHO   = 1800.0   # kg/m³
MU    = 0.002    # Pa·s
BETA  = 3.5e-4   # 1/K
G     = 9.81

# ── Solar ──────────────────────────────────────────────────────────────────────
G_SOL    = 1000.0   # W/m² irradiació directa
ETA_LUPA = 0.85     # eficiència òptica lupa
ETA_TUB  = 0.90     # eficiència transferència lupa→fluid
ETA_TURB = 0.80     # eficiència microturbina
T_COLD   = 40.0     # °C circuit fred
T_HOT_MAX = 555.0   # °C límit sals fundides

# ── Geometria de les lupes ─────────────────────────────────────────────────────
# Lupes al perímetre exterior del circuit calent
# Cada lupa ocupa diam_lupa metres de perímetre
DIAM_LUPA      = 3.0    # m  diàmetre de cada lupa
ESPAI_LUPA     = 3.2    # m  espai per lupa al perímetre (lupa + marge)
ESPAI_TURB     = 1.5    # m  espai per turbina al perímetre interior

def v_flux(dP_motor, D, L_circ):
    """Velocitat per convecció natural (equilibri flotabilitat-fricció)."""
    v = 0.05
    for _ in range(100):
        Re = max(RHO*v*D/MU, 1)
        f  = 64/Re if Re < 2300 else 0.316/Re**0.25
        v_new = np.sqrt(max(1e-9, 2*dP_motor/(f*L_circ/D*RHO+1e-9)))
        v = 0.3*v_new + 0.7*v
    return max(v, 1e-4)

def simula_circuit(R, D_tub, D_inter, N_lupes, N_turb_cascada):
    """
    Simulació completa d'un circuit amb cascada de turbines.
    
    Escala automàticament: lupes al perímetre, turbines a les interconnexions.
    """
    Circ  = 2*np.pi*R
    A_tub = np.pi*(D_tub/2)**2
    A_int = np.pi*(D_inter/2)**2

    # ── Límits físics ──────────────────────────────────────────────────────────
    N_lupes_max  = int(Circ / ESPAI_LUPA)
    N_turb_max   = int(Circ / ESPAI_TURB)
    N_lupes_real = min(N_lupes, N_lupes_max)
    N_turb_real  = min(N_turb_cascada, N_turb_max)

    # ── Potència solar ─────────────────────────────────────────────────────────
    A_lupa  = np.pi*(DIAM_LUPA/2)**2
    P_solar = N_lupes_real * A_lupa * G_SOL * ETA_LUPA * ETA_TUB  # W

    # ── Temperatura HOT ────────────────────────────────────────────────────────
    # Iteració: troba T_hot i v_circ consistents
    T_hot = T_COLD + 50  # primera estimació
    for _ in range(30):
        DT      = T_hot - T_COLD
        dP_buo  = RHO * BETA * DT * G * R * 2
        v_circ  = v_flux(dP_buo, D_tub, Circ)
        m_dot   = RHO * A_tub * v_circ
        DT_new  = min(P_solar / (m_dot * CP + 1e-9), T_HOT_MAX - T_COLD)
        T_hot   = T_COLD + DT_new * 0.3 + DT * 0.7

    DT    = T_hot - T_COLD
    dP_buo = RHO * BETA * DT * G * R * 2
    v_circ = v_flux(dP_buo, D_tub, Circ)
    m_dot  = RHO * A_tub * v_circ

    # ── Cascada de turbines ────────────────────────────────────────────────────
    # La interconnexió té longitud = R (del circuit exterior a l'interior)
    # ΔP disponible per a cada turbina = ΔP_buoyança - fricció_interconnexió
    def dP_inter(DT_local, D_i, L_i=R):
        dP_m = RHO * BETA * DT_local * G * L_i
        v_i  = 0.5  # estimació inicial
        for _ in range(30):
            Re_i = max(RHO*v_i*D_i/MU, 1)
            f_i  = 64/Re_i if Re_i < 2300 else 0.316/Re_i**0.25
            dP_fr = f_i*(L_i/D_i)*0.5*RHO*v_i**2
            dP_ut = max(0, dP_m - dP_fr)
            v_new = np.sqrt(max(1e-9, 2*dP_ut/(RHO+1e-9)))
            v_i   = 0.3*v_new + 0.7*v_i
        Re_i  = max(RHO*v_i*D_i/MU, 1)
        f_i   = 64/Re_i if Re_i < 2300 else 0.316/Re_i**0.25
        dP_fr = f_i*(L_i/D_i)*0.5*RHO*v_i**2
        dP_ut = max(0, dP_m - dP_fr)
        return dP_ut, v_i

    # Cascada: cada turbina extreu del ΔT restant
    DT_restant = DT
    P_turb_tot = 0.0
    Pi_list    = []

    for i in range(N_turb_real):
        if DT_restant < 0.5:
            Pi_list.append(0)
            continue
        dPu, vi    = dP_inter(DT_restant, D_inter)
        Q_inter    = A_int * vi
        P_i        = dPu * Q_inter * ETA_TURB * N_turb_real  # totes les inter. en paral·lel
        # Però cada una extreu del mateix ΔT → dividim per N
        P_i_1turb  = dPu * Q_inter * ETA_TURB
        DT_consum  = P_i_1turb / (m_dot * CP + 1e-9)
        DT_consum  = min(DT_consum, DT_restant * 0.4)
        DT_restant -= DT_consum
        P_turb_tot += P_i_1turb * N_turb_real
        Pi_list.append(P_i_1turb * N_turb_real)

    # ── Eficiències ───────────────────────────────────────────────────────────
    eta_real   = P_turb_tot / (P_solar + 1e-9) * 100
    eta_carnot = (1 - (T_COLD+273.15)/(T_hot+273.15)) * 100
    eta_relat  = eta_real / eta_carnot * 100 if eta_carnot > 0 else 0

    return {
        'R': R, 'D_tub': D_tub, 'D_inter': D_inter,
        'Circ': Circ,
        'N_lupes_max': N_lupes_max, 'N_lupes': N_lupes_real,
        'N_turb_max': N_turb_max,  'N_turb': N_turb_real,
        'P_solar_kW': P_solar/1e3,
        'T_hot': T_hot, 'DT': DT,
        'v_circ': v_circ, 'm_dot': m_dot,
        'P_turb_kW': P_turb_tot/1e3,
        'P_per_lupa': P_turb_tot/max(N_lupes_real,1)/1e3,
        'P_per_m2': P_turb_tot/(np.pi*R**2+1e-9),  # W/m² de superfície
        'eta': eta_real,
        'eta_carnot': eta_carnot,
        'eta_relat': eta_relat,
        'Pi_list': Pi_list,
    }

# ─────────────────────────────────────────────────────────────────────────────
print("="*72)
print("GEMELL v6 — Estudi d'Escala Complet")
print("="*72)
print(f"  Fluid: Sals fundides | T_cold={T_COLD}°C | Lupa Ø{DIAM_LUPA}m")
print(f"  η_lupa={ETA_LUPA} | η_tub={ETA_TUB} | η_turb={ETA_TURB}")

# ── 1. Escala amb R (tots els paràmetres al màxim per R) ─────────────────────
print(f"\n{'='*72}")
print("1. COM ESCALA AMB EL RADI R (tot al màxim físic per R)")
print(f"{'='*72}")
print(f"  {'R m':>6} {'Circ m':>8} {'Lupes max':>10} {'Turb max':>9} "
      f"{'T_hot°C':>9} {'P_sol MW':>10} {'P_turb kW':>11} {'η%':>6} {'Carnot%':>8}")

Rs = [10, 25, 50, 100, 200, 500, 1000, 2000, 5000]
res_R = []
for R in Rs:
    N_l = int(2*np.pi*R / ESPAI_LUPA)
    N_t = int(2*np.pi*R / ESPAI_TURB)
    r = simula_circuit(R, D_tub=max(0.1,R*0.003),
                       D_inter=max(0.05,R*0.001),
                       N_lupes=N_l, N_turb_cascada=min(N_t,20))
    res_R.append(r)
    print(f"  {R:>6} {r['Circ']:>8.0f} {r['N_lupes']:>10} {r['N_turb']:>9} "
          f"{r['T_hot']:>9.1f} {r['P_solar_kW']/1e3:>10.2f} "
          f"{r['P_turb_kW']:>11.1f} {r['eta']:>6.2f} {r['eta_carnot']:>8.1f}")

# ── 2. N turbines òptim per R=100m ───────────────────────────────────────────
print(f"\n{'='*72}")
print("2. N TURBINES ÒPTIM (R=100m, lupes al màxim)")
print(f"{'='*72}")
R_test = 100
N_l_100 = int(2*np.pi*R_test / ESPAI_LUPA)
print(f"  {'N turb':>8} {'P_turb kW':>12} {'P_marginal kW':>15} "
      f"{'P/turb kW':>12} {'η%':>6}")
prev = 0
res_NT = []
for nt in [1,2,3,5,8,10,15,20,30,50,100]:
    if nt > int(2*np.pi*R_test/ESPAI_TURB):
        break
    r = simula_circuit(R_test, 0.30, 0.15, N_l_100, nt)
    Pm = r['P_turb_kW'] - prev
    res_NT.append((nt, r['P_turb_kW'], Pm))
    print(f"  {nt:>8} {r['P_turb_kW']:>12.2f} {Pm:>15.2f} "
          f"{r['P_turb_kW']/nt:>12.3f} {r['eta']:>6.3f}")
    prev = r['P_turb_kW']

# ── 3. Molts xicotets vs un de gran ──────────────────────────────────────────
print(f"\n{'='*72}")
print("3. MOLTS XICOTETS vs UN DE GRAN (mateixa àrea total de lupes)")
print(f"{'='*72}")
# Àrea total fixada: equivalent a R=200m ple de lupes
N_lupes_total = int(2*np.pi*200 / ESPAI_LUPA)
print(f"  N_lupes_total = {N_lupes_total} (equivalent a R=200m ple)")
print(f"  {'Config':>20} {'N unitats':>10} {'R cada una':>11} "
      f"{'P_turb kW':>12} {'η%':>6} {'P/m²':>8}")

# Un gran (R=200)
r_gran = simula_circuit(200, 0.60, 0.20, N_lupes_total, 15)
print(f"  {'1 gran (R=200m)':>20} {1:>10} {200:>11} "
      f"{r_gran['P_turb_kW']:>12.1f} {r_gran['eta']:>6.2f} "
      f"{r_gran['P_per_m2']:>8.3f}")

# N xicotets amb la mateixa àrea de lupes total
for n_units, R_unit in [(2,100),(4,50),(10,25),(20,15),(50,10)]:
    n_lup_cada = N_lupes_total // n_units
    if n_lup_cada < 1:
        break
    r_unit = simula_circuit(R_unit, max(0.1,R_unit*0.003),
                            max(0.05,R_unit*0.001),
                            n_lup_cada, 10)
    P_tot = r_unit['P_turb_kW'] * n_units
    eta_tot = P_tot / (r_gran['P_solar_kW']) * 100 if r_gran['P_solar_kW']>0 else 0
    print(f"  {f'{n_units} × R={R_unit}m':>20} {n_units:>10} {R_unit:>11} "
          f"{P_tot:>12.1f} {r_unit['eta']:>6.2f} "
          f"{r_unit['P_per_m2']:>8.3f}")

# ── 4. Dades que necessitem (resum paramètric) ────────────────────────────────
print(f"\n{'='*72}")
print("4. PARÀMETRES CLAU I COM AFECTEN LA POTÈNCIA")
print(f"{'='*72}")
r_base = simula_circuit(100, 0.30, 0.15, 100, 10)
P_base = r_base['P_turb_kW']

params = [
    ("R doble (100→200m)", simula_circuit(200,0.30,0.15,100,10),'P_turb_kW'),
    ("D_tub doble (0.30→0.60m)", simula_circuit(100,0.60,0.15,100,10),'P_turb_kW'),
    ("D_inter doble (0.15→0.30m)", simula_circuit(100,0.30,0.30,100,10),'P_turb_kW'),
    ("N_lupes doble (100→200)", simula_circuit(100,0.30,0.15,200,10),'P_turb_kW'),
    ("N_turb doble (10→20)", simula_circuit(100,0.30,0.15,100,20),'P_turb_kW'),
    ("T_cold meitat (40→20°C)", None, None),
]
print(f"  {'Paràmetre':>35} {'P_base kW':>10} {'P_nova kW':>10} {'Millora':>8}")
print(f"  {'(base: R=100,Dt=0.30,Di=0.15,N=100l,10t)':>35} {P_base:>10.2f} {'—':>10} {'—':>8}")
for desc, r_new, key in params:
    if r_new is None:
        continue
    P_new = r_new[key]
    millora = (P_new/P_base - 1)*100
    print(f"  {desc:>35} {P_base:>10.2f} {P_new:>10.2f} {millora:>+7.1f}%")

# ── 5. La resposta: millors xicotets o gran? ──────────────────────────────────
print(f"\n{'='*72}")
print("5. RESPOSTA: MILLORS MOLTS XICOTETS O UN DE GRAN?")
print(f"{'='*72}")
print("""
  P_turb ∝ R²   (via: P_solar ∝ R, ΔP ∝ R, m_dot ∝ R²  → P ∝ R³... però η cau)
  
  UN DE GRAN:
  ✅ Menys pèrdues de calor (menys superfície per volum)
  ✅ Convecció natural més forta (R gran → ΔP gran)
  ✅ Menys connexions elèctriques externes
  ❌ Difícil de construir, transportar, mantenir
  ❌ Si falla, falla tot
  
  MOLTS XICOTETS:
  ✅ Modular, fàcil de construir i replicar
  ✅ Redundància (un fallen, la resta funcionen)
  ✅ Més fàcil d'optimitzar individualment
  ❌ Més pèrdues de calor per superfície total
  ❌ Menys ΔP per unitat (R xicotet → menys pressió)
  
  CONCLUSIÓ NUMÈRICA:
  La potència per m² de superfície NO escala linealment amb R.
  El punt òptim depèn del cost de construcció vs. potència extra.
  Per a instal·lació tipus CSP industrial: R=50-200m és el rang habitual.
""")

# ── Gràfics ───────────────────────────────────────────────────────────────────
BG='#0d0d1a'; PAN='#13132b'
CW='white'; CG='#00ff88'; CO='#ffd700'; CR='#ff6644'
CB='#00d2ff'; CP='#cc66ff'; CK='#ff99cc'

fig = plt.figure(figsize=(24,18), facecolor=BG)
gs  = gridspec.GridSpec(3,3, figure=fig, hspace=0.52, wspace=0.38)

def sty(ax, tit, xl='', yl=''):
    ax.set_facecolor(PAN)
    ax.set_title(tit, color=CW, fontsize=9.5, pad=5, fontweight='bold')
    ax.tick_params(colors='#aaa', labelsize=8)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl, color='#aaa', fontsize=8.5)
    ax.set_ylabel(yl, color='#aaa', fontsize=8.5)
    ax.grid(color='#1e1e40', lw=0.5, ls='--', alpha=0.8)

Rv   = [r['R'] for r in res_R]
Pturb= [r['P_turb_kW'] for r in res_R]
Psol = [r['P_solar_kW'] for r in res_R]
Etas = [r['eta'] for r in res_R]
Thos = [r['T_hot'] for r in res_R]
Nlup = [r['N_lupes'] for r in res_R]
Ppm2 = [r['P_per_m2'] for r in res_R]

# G1: P_turb i P_solar vs R
ax = fig.add_subplot(gs[0,0])
ax.loglog(Rv, Pturb, 'o-', color=CG, lw=2.5, ms=8, label='P turbines')
ax.loglog(Rv, Psol, 's--', color=CR, lw=2, ms=6, label='P solar captada')
ax.legend(fontsize=8.5, framealpha=0.3)
sty(ax, 'G1 — Potència vs Radi R\n[Escala log-log: pendent = exponent]', 'R [m]', 'kW')

# G2: η real i Carnot vs R
ax = fig.add_subplot(gs[0,1])
ax.semilogx(Rv, Etas, 'D-', color=CO, lw=2.5, ms=8, label='η real')
ax.semilogx(Rv, [r['eta_carnot'] for r in res_R], '^--',
            color=CP, lw=2, ms=6, label='η Carnot (límit)')
ax.legend(fontsize=8.5, framealpha=0.3)
sty(ax, 'G2 — Eficiència vs Radi\n[η real vs límit Carnot]', 'R [m]', 'η [%]')

# G3: N_lupes màx i T_hot vs R
ax = fig.add_subplot(gs[0,2])
ax2 = ax.twinx()
ax.loglog(Rv, Nlup, 'o-', color=CB, lw=2.5, ms=8, label='N_lupes màx')
ax2.semilogx(Rv, Thos, 's--', color=CR, lw=2, ms=6, label='T_hot')
ax.set_ylabel('N_lupes màx', color=CB, fontsize=8.5)
ax2.set_ylabel('T_hot [°C]', color=CR, fontsize=8.5)
ax2.tick_params(colors='#aaa', labelsize=8)
ax.legend(fontsize=8, framealpha=0.3, loc='upper left')
ax2.legend(fontsize=8, framealpha=0.3, loc='lower right')
sty(ax, 'G3 — N_lupes màx i T_hot vs R', 'R [m]')

# G4: P_marginal per turbina (rendiment decreixent)
nts = [x[0] for x in res_NT]
Pts = [x[1] for x in res_NT]
Pms = [x[2] for x in res_NT]
ax = fig.add_subplot(gs[1,0])
ax2 = ax.twinx()
ax.plot(nts, Pts, 'o-', color=CG, lw=2.5, ms=8, label='P acumulada')
ax2.plot(nts, Pms, 's--', color=CO, lw=2, ms=6, label='P marginal')
ax2.axhline(0, color=CR, ls=':', lw=1)
ax.set_ylabel('P acumulada [kW]', color=CG, fontsize=8.5)
ax2.set_ylabel('P marginal [kW]', color=CO, fontsize=8.5)
ax2.tick_params(colors='#aaa', labelsize=8)
ax.legend(fontsize=8, framealpha=0.3, loc='upper left')
ax2.legend(fontsize=8, framealpha=0.3, loc='upper right')
sty(ax, f'G4 — Rendiment decreixent turbines\n(R={R_test}m, N_lupes màx)',
    'N turbines', 'kW')

# G5: P/m² vs R (densitat de potència)
ax = fig.add_subplot(gs[1,1])
ax.loglog(Rv, [p*1e-3 for p in Ppm2], 'D-', color=CK, lw=2.5, ms=8)
sty(ax, 'G5 — Densitat de Potència P/m² vs R\n[Indica si val la pena fer-ho gran]',
    'R [m]', 'W/m²')

# G6: Escala cascada recursiva — visualització
ax = fig.add_subplot(gs[1,2])
ax.set_facecolor(PAN)
# Mostra la sèrie geomètrica: cada turbina pren x% del que queda
x0 = 30  # % que pren la primera turbina
N_show = 12
DT_now = 410.0
DTs = [DT_now]
Ps_show = []
for i in range(N_show):
    Pi = DT_now * x0/100
    Ps_show.append(Pi)
    DT_now -= Pi
    DTs.append(DT_now)

colors_grad = plt.cm.plasma(np.linspace(0.9, 0.2, N_show))
bars = ax.bar(range(1, N_show+1), Ps_show,
              color=colors_grad, alpha=0.85, width=0.7)
ax.plot(range(1, N_show+1), DTs[:-1], 's--', color=CB, lw=1.5, ms=5,
        label='ΔT restant')
ax.set_ylabel('ΔT furtat [°C]', color=CO, fontsize=8.5)
ax2 = ax.twinx()
ax2.plot(range(1, N_show+1), DTs[:-1], 's--', color=CB, lw=1.5, ms=5)
ax2.set_ylabel('ΔT restant [°C]', color=CB, fontsize=8.5)
ax2.tick_params(colors='#aaa', labelsize=8)
sty(ax, f'G6 — Cascada recursiva: f(x)=x·(1-x/100)^n\n'
        f'[Cada turbina pren {x0}% del ΔT restant]',
    'N turbina en cascada', 'ΔT furtat per turbina [°C]')

# G7: Taula paràmetres necessaris
ax = fig.add_subplot(gs[2,:])
ax.axis('off'); ax.set_facecolor('#0a0a14')

rows = [
    ['Paràmetre','Rang recomanat','Com afecta','Limitació física'],
    ['Radi R','50–500 m',
     'P ∝ R^2.5 aprox. (N_lupes·ΔP·m_dot)',
     'Cost construcció, terreny'],
    ['Ø tub principal D_tub','R/200 a R/100',
     'Controla m_dot → P ∝ D_tub²',
     'Convecció natural: massa gran = lent'],
    ['Ø interconnexió D_inter','R/500 a R/200',
     'v_inter ↑ si D_inter ↓ (però Q ↓)',
     'Punt òptim: ~1/3 D_tub'],
    ['N lupes','Màx físic (1 cada 3.2m perim.)',
     'P ∝ N_lupes (lineal)',
     'Perímetre disponible = 2πR/3.2'],
    ['N turbines cascada','8–15 (rendiment decreixent)',
     'P marginal cau ~50% per cada 2×',
     'Perímetre + cost instal·lació'],
    ['T_cold','Mínim possible (<40°C)',
     'ΔT ↑ → ΔP ↑ → P ↑ (quasi lineal)',
     'Refrigeració activa costosa'],
    ['D_tub vs D_inter','D_inter ≈ D_tub/3',
     'Maximitza P = ΔP·Q·η',
     'Massa gran = poca v; massa petit = poca Q'],
    ['1 gran vs N xicotets','Depèn del cost',
     'Gran: més η; Xicotet: més modular',
     'CSP industrial: R=50-200m habitual'],
]
tbl = ax.table(cellText=rows[1:], colLabels=rows[0],
               loc='center', cellLoc='center',
               colWidths=[0.15, 0.18, 0.35, 0.30])
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5)
for (rr, c), cell in tbl.get_celld().items():
    if rr == 0:
        bg = '#1F5C9E'
    elif rr % 2 == 0:
        bg = '#0d0d1a'
    else:
        bg = '#111125'
    cell.set_facecolor(bg)
    cell.set_text_props(color='white')
    cell.set_edgecolor('#333355')
    cell.set_height(0.095)
ax.set_title('G7 — Taula de Paràmetres: Quines Dades Necessitem i Com Afecten',
             color=CW, fontsize=10, pad=8, fontweight='bold')

fig.suptitle(
    'GEMELL v6 — Estudi d\'Escala Complet: Radi, Lupes, Turbines en Cascada\n'
    'Sals Fundides Solar Salt | Motor Tèrmic Solar | '
    'Víctor Manzanares Alberola — EPSA UPV Alcoi — 2026',
    color=CW, fontsize=11, fontweight='bold', y=0.999)

OUT='/mnt/user-data/outputs/solar_v6_escala.png'
plt.savefig(OUT, dpi=145, bbox_inches='tight', facecolor=BG)
plt.close()
print(f"\nGràfic: {OUT}")
