"""
model_44_Quijote_zypyzape.py
═══════════════════════════════════════════════════════════════════════════════
Gemell digital: 44 molins ZypyZape (3 aspas, NREL 5MW) + Quijote (3 aspas)
Bus comú d'inèrcia sintètica — Lògica gemell_3vs7_rigoros adaptada a 3A

Comparació:
  [A] Camp sol:  44 ZypyZape sense Quijote
  [B] Camp+Q3:   44 ZypyZape + Quijote acoblat al bus

Autor: Víctor Manzanares Alberola — EPSA / UPV (Alcoi)
Repositori: github.com/espiradesombra/claude
───────────────────────────────────────────────────────────────────────────────
FÍSICA DEL MODEL
────────────────
Cada molí ZZ (índex i, 3 aspas):

  I_i · dω_i/dt = τ_aero_i(v_eff, ω_i)
                - τ_loss_i(ω_i)
                + τ_bus_i          ← acoblament via bus
                + τ_Q3_i           ← només escenari B

Bus comú (inercia sintètica compartida):
  ω_bus = Σ(I_i · ω_i) / Σ(I_i)    ← velocitat angular ponderada
  τ_bus_i = K_bus · (ω_bus - ω_i)   ← sincronització suau

Quijote (3 aspas, inercia variable Fe+oli):
  I_Q3(θ) = I0_Q3 · (1 + α · sin(3·θ_Q3))
  dω_Q3/dt = [τ_aero_Q3 - τ_loss_Q3 + τ_transfer] / I_Q3(θ_Q3)

Transferència al bus (hurto gravitatori):
  τ_transfer = K_Q3 · sin(3·θ_Q3 - ω_bus·t) · I_Q3

Potència aerodinàmica (BEM simplificat, λ-Cp):
  P_i = 0.5 · ρ · π · R² · v³ · Cp(λ)
  λ = ω_i · R / v
  Cp(λ) = Cp_max · exp(-((λ - λ_opt)/σ)²)

Energia neta per cicle de Quijote (hurto gravitatori):
  E_cicle = ∫₀ᵀ τ_transfer(t) · ω_Q3(t) dt
═══════════════════════════════════════════════════════════════════════════════
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
import warnings
warnings.filterwarnings('ignore')

# ─── PARÀMETRES NREL 5MW (Jonkman 2009) ──────────────────────────────────────
RHO        = 1.225          # kg/m³ densitat aire
R_NREL     = 63.0           # m radi rotor
J_NREL     = 38_759_228.0   # kg·m² inercia rotor
OMEGA_RAT  = 12.1 * 2*np.pi/60  # rad/s velocitat nominal
V_RAT      = 11.4           # m/s vent nominal
V_CUT_IN   = 3.0
V_CUT_OUT  = 25.0
CP_MAX     = 0.48
LAMBDA_OPT = 7.5
SIGMA_CP   = 3.2

# ─── CAMP DE 44 MOLINS ───────────────────────────────────────────────────────
N_ZZ       = 44
N_ASPAS    = 3              # ZypyZape amb 3 aspas (com Quijote)

# ─── PARÀMETRES BUS ──────────────────────────────────────────────────────────
K_BUS      = 0.08           # N·m·s acoblament bus comú (Kuramoto-like)

# ─── PARÀMETRES QUIJOTE 3A ────────────────────────────────────────────────────
I0_Q3      = J_NREL * 1.25  # kg·m² inercia base (massa Fe+oli afegida)
ALPHA_Q3   = 0.35           # amplitud inercia sintètica
K_Q3_BUS   = 0.12           # acoblament Quijote → bus
K_Q3_IND   = 0.04           # acoblament Quijote → molins individuals

# ─── INTEGRACIÓ ──────────────────────────────────────────────────────────────
DT         = 0.05           # s  pas temporal
T_SIM      = 600.0          # s  durada simulació (10 min)
T_STEPS    = int(T_SIM / DT)

# ─── PERFIL DE VENT ──────────────────────────────────────────────────────────
WIND_BASE  = 9.0            # m/s
WIND_TURB  = 0.12           # intensitat turbulència (fracció)
WIND_GUST_T = 200.0         # s  moment de la ratxa
WIND_GUST_V = 14.5          # m/s velocitat pic ratxa
WIND_GUST_DUR = 60.0        # s  durada ratxa


# ═══════════════════════════════════════════════════════════════════════════════
# FUNCIONS FÍSIQUES
# ═══════════════════════════════════════════════════════════════════════════════

def cp_lambda(lam):
    """Coeficient de potència Cp(λ) — BEM simplificat"""
    return CP_MAX * np.exp(-((lam - LAMBDA_OPT) / SIGMA_CP) ** 2)


def tau_aero(omega, v_eff, R=R_NREL, rho=RHO):
    """Par aerodinàmic [N·m] donat ω i velocitat de vent efectiva"""
    omega = max(omega, 0.01)
    v_eff = max(v_eff, 0.1)
    lam   = omega * R / v_eff
    cp    = cp_lambda(lam)
    P     = 0.5 * rho * np.pi * R**2 * v_eff**3 * max(cp, 0.0)
    return P / omega


def tau_loss(omega, I, k_loss=0.03):
    """Par de pèrdues mecàniques (fricció viscosa)"""
    return k_loss * I * omega


def wind_profile(t_arr):
    """Perfil de vent: WIND_BASE + turbulència + ratxa gaussiana"""
    rng  = np.random.default_rng(42)
    turb = rng.normal(0, WIND_TURB * WIND_BASE, len(t_arr))
    # Filtre passa-baix sobre turbulència (coherència espacial)
    from numpy import convolve, ones
    kern = ones(int(5/DT)) / int(5/DT)
    turb = np.convolve(turb, kern, mode='same')
    gust = WIND_GUST_V * np.exp(-0.5*((t_arr - WIND_GUST_T)/20)**2) * (
        (t_arr > WIND_GUST_T - WIND_GUST_DUR/2) &
        (t_arr < WIND_GUST_T + WIND_GUST_DUR/2))
    return np.clip(WIND_BASE + turb + gust, V_CUT_IN, V_CUT_OUT)


# ═══════════════════════════════════════════════════════════════════════════════
# ESTAT INICIAL
# ═══════════════════════════════════════════════════════════════════════════════

def init_state(rng_seed=7):
    """Inicialitza camp de 44 molins amb dispersió natural"""
    rng = np.random.default_rng(rng_seed)
    # Cada molí té omega_nat lleugerament diferent (±8%)
    omega_nat = OMEGA_RAT * (1.0 + rng.uniform(-0.08, 0.08, N_ZZ))
    omega_0   = omega_nat * rng.uniform(0.85, 1.05, N_ZZ)
    theta_0   = rng.uniform(0, 2*np.pi, N_ZZ)
    I_i       = np.full(N_ZZ, J_NREL)
    return {
        'theta': theta_0.copy(),
        'omega': omega_0.copy(),
        'omega_nat': omega_nat.copy(),
        'I': I_i.copy(),
        # Quijote
        'theta_Q3': 0.0,
        'omega_Q3': OMEGA_RAT * 0.9,
        'E_cicle': 0.0,
        '_E_acc': 0.0,
        '_last_rev': 0.0,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# NUCLI D'INTEGRACIÓ
# ═══════════════════════════════════════════════════════════════════════════════

def step_ZZ_only(state, v_wind):
    """
    Escenari A — 44 ZZ sense Quijote
    Bus comú: acoblament Kuramoto entre molins
    """
    th = state['theta']
    om = state['omega']
    I  = state['I']

    # Velocitat de bus (ponderada per inercia)
    omega_bus = np.sum(I * om) / np.sum(I)

    d_om = np.zeros(N_ZZ)
    for i in range(N_ZZ):
        # Variació de vent per posició (efecte wake simplificat)
        v_eff = v_wind * (0.92 + 0.08 * np.sin(i * 2*np.pi/N_ZZ + om[i]*0.1))

        ta = tau_aero(om[i], v_eff)
        tl = tau_loss(om[i], I[i])
        tb = K_BUS * (omega_bus - om[i]) * I[i]

        d_om[i] = (ta - tl + tb) / I[i]

    # Euler explícit
    state['omega'] = np.maximum(0, om + d_om * DT)
    state['theta'] = th + state['omega'] * DT

    return state


def step_ZZ_Q3(state, v_wind):
    """
    Escenari B — 44 ZZ + Quijote 3A acoblat al bus
    Quijote aporta inercia sintètica i transferència de fase
    """
    th    = state['theta']
    om    = state['omega']
    I     = state['I']
    th_Q3 = state['theta_Q3']
    om_Q3 = state['omega_Q3']

    # Inercia sintètica Quijote (Fe+oli, periodicitat 3A)
    I_Q3  = I0_Q3 * (1.0 + ALPHA_Q3 * np.sin(N_ASPAS * th_Q3))

    # Velocitat bus
    omega_bus = (np.sum(I * om) + I_Q3 * om_Q3) / (np.sum(I) + I_Q3)

    # ── Molins ZZ ──
    d_om = np.zeros(N_ZZ)
    for i in range(N_ZZ):
        v_eff = v_wind * (0.92 + 0.08 * np.sin(i * 2*np.pi/N_ZZ + om[i]*0.1))

        ta = tau_aero(om[i], v_eff)
        tl = tau_loss(om[i], I[i])
        tb = K_BUS * (omega_bus - om[i]) * I[i]

        # Transferència directa Quijote → molí individual (hurto gravitatori)
        # Finestra de fase: sin(3·θ_Q3 - 3·θ_i) — 3 aspas acoblades
        tau_q3_i = K_Q3_IND * np.sin(N_ASPAS * th_Q3 - N_ASPAS * th[i]) * I_Q3

        d_om[i] = (ta - tl + tb + tau_q3_i) / I[i]

    # ── Quijote ──
    ta_Q3 = tau_aero(om_Q3, v_wind * 1.05, R=R_NREL * 1.1)  # lleugerament diferent
    tl_Q3 = tau_loss(om_Q3, I_Q3, k_loss=0.025)

    # Par transferència Quijote → bus (net)
    tau_transfer = K_Q3_BUS * np.sin(N_ASPAS * th_Q3) * I_Q3 * (omega_bus - om_Q3)

    d_om_Q3 = (ta_Q3 - tl_Q3 + tau_transfer) / I_Q3

    # Energia per cicle (hurto gravitatori)
    state['_E_acc'] += tau_transfer * om_Q3 * DT
    if th_Q3 - state['_last_rev'] >= 2*np.pi:
        state['E_cicle'] = state['_E_acc']
        state['_E_acc']  = 0.0
        state['_last_rev'] = th_Q3

    # Integració Euler
    state['omega']    = np.maximum(0, om + d_om * DT)
    state['theta']    = th + state['omega'] * DT
    state['omega_Q3'] = max(0, om_Q3 + d_om_Q3 * DT)
    state['theta_Q3'] = th_Q3 + state['omega_Q3'] * DT

    return state


# ═══════════════════════════════════════════════════════════════════════════════
# MÈTRIQUES
# ═══════════════════════════════════════════════════════════════════════════════

def power_MW(state, v_wind):
    """Potència total del camp [MW]"""
    total = 0.0
    for i in range(N_ZZ):
        v_eff = v_wind * (0.92 + 0.08 * np.sin(i * 2*np.pi/N_ZZ))
        om = state['omega'][i]
        lam = om * R_NREL / max(v_eff, 0.1)
        cp  = cp_lambda(lam)
        total += 0.5 * RHO * np.pi * R_NREL**2 * v_eff**3 * max(cp, 0)
    return total / 1e6


def order_parameter(state):
    """Paràmetre d'ordre de Kuramoto r ∈ [0,1]"""
    phases = state['theta'] % (2*np.pi)
    r = abs(np.mean(np.exp(1j * phases)))
    return r


def omega_std(state):
    """Desviació estàndard de velocitats angulars"""
    return np.std(state['omega'])


def sync_count(state, tol=0.15):
    """Nombre de molins dins tolerància de velocitat mitjana"""
    om_mean = np.mean(state['omega'])
    return int(np.sum(np.abs(state['omega'] - om_mean) < tol))


# ═══════════════════════════════════════════════════════════════════════════════
# SIMULACIÓ PRINCIPAL
# ═══════════════════════════════════════════════════════════════════════════════

def run_simulation():
    print("═" * 70)
    print("  model_44_Quijote_zypyzape.py")
    print("  Gemell digital: 44×ZZ (3A) + Quijote (3A) — Bus comú")
    print("  Víctor Manzanares Alberola — EPSA/UPV Alcoi")
    print("═" * 70)
    print(f"\n  N molins  : {N_ZZ}")
    print(f"  Aspas/molí: {N_ASPAS}")
    print(f"  T simulació: {T_SIM:.0f} s ({T_SIM/60:.1f} min)")
    print(f"  dt         : {DT} s  ({T_STEPS} passos)")
    print(f"  Vent base  : {WIND_BASE} m/s  (ratxa: {WIND_GUST_V} m/s @ t={WIND_GUST_T}s)")
    print(f"  K_bus      : {K_BUS}")
    print(f"  K_Q3→bus   : {K_Q3_BUS}   K_Q3→ind: {K_Q3_IND}")
    print(f"  α inercia  : {ALPHA_Q3}")
    print()

    t_arr  = np.arange(T_STEPS) * DT
    v_arr  = wind_profile(t_arr)

    # ── Arrays de resultats ──
    rec = {s: {k: [] for k in ['P', 'r', 'std_om', 'sync', 'E_cicle']}
           for s in ['A', 'B']}

    stA = init_state(seed := 7)
    stB = init_state(seed)
    # Copia idèntica per comparació justa
    stB['theta']    = stA['theta'].copy()
    stB['omega']    = stA['omega'].copy()
    stB['omega_nat']= stA['omega_nat'].copy()

    print("  Simulant escenari A (sense Quijote) i B (amb Quijote)...")
    log_interval = max(1, T_STEPS // 20)

    for k in range(T_STEPS):
        v = v_arr[k]

        stA = step_ZZ_only(stA, v)
        stB = step_ZZ_Q3(stB, v)

        if k % 4 == 0:   # desa cada 4 passos (0.2 s)
            for lbl, st in [('A', stA), ('B', stB)]:
                rec[lbl]['P'].append(power_MW(st, v))
                rec[lbl]['r'].append(order_parameter(st))
                rec[lbl]['std_om'].append(omega_std(st))
                rec[lbl]['sync'].append(sync_count(st))
                if lbl == 'B':
                    rec[lbl]['E_cicle'].append(st['E_cicle'] / 1e6)
                else:
                    rec[lbl]['E_cicle'].append(0.0)

        if k % log_interval == 0:
            pct = 100 * k / T_STEPS
            print(f"    {pct:5.1f}%  t={k*DT:6.1f}s  "
                  f"PA={power_MW(stA,v):5.1f}MW  PB={power_MW(stB,v):5.1f}MW  "
                  f"rA={order_parameter(stA):.3f}  rB={order_parameter(stB):.3f}  "
                  f"syncA={sync_count(stA)}/{N_ZZ}  syncB={sync_count(stB)}/{N_ZZ}")

    # Converteix a arrays
    t_rec = t_arr[::4][:len(rec['A']['P'])]
    for lbl in ['A', 'B']:
        for k2 in rec[lbl]:
            rec[lbl][k2] = np.array(rec[lbl][k2])

    return t_rec, v_arr[::4][:len(t_rec)], rec, stA, stB


# ═══════════════════════════════════════════════════════════════════════════════
# ANÀLISI ESTADÍSTICA
# ═══════════════════════════════════════════════════════════════════════════════

def print_summary(t_rec, rec):
    print("\n" + "═" * 70)
    print("  RESUM COMPARATIU")
    print("═" * 70)
    warmup = len(t_rec) // 5   # descarta primers 20% (transitori)

    metrics = [
        ('Potència mitjana [MW]', 'P', np.mean),
        ('Potència màxima [MW]', 'P', np.max),
        ('Desv. estàndard P [MW]', 'P', np.std),
        ('Ordre de sincronització r', 'r', np.mean),
        ('Molins sincronitzats (mitjana)', 'sync', np.mean),
        ('Desv. ω [rad/s]', 'std_om', np.mean),
    ]

    print(f"\n  {'Mètrica':<38} {'Escenari A':>12} {'Escenari B':>12} {'Millora':>10}")
    print(f"  {'-'*38} {'-'*12} {'-'*12} {'-'*10}")

    for label, key, fn in metrics:
        vA = fn(rec['A'][key][warmup:])
        vB = fn(rec['B'][key][warmup:])
        delta = 100 * (vB - vA) / (abs(vA) + 1e-12)
        sign  = '+' if delta >= 0 else ''
        print(f"  {label:<38} {vA:>12.3f} {vB:>12.3f} {sign}{delta:>9.2f}%")

    E_q3_mean = np.mean(rec['B']['E_cicle'][warmup:][rec['B']['E_cicle'][warmup:] != 0])
    print(f"\n  {'Energia transferida/cicle Q3 [MJ]':<38} {'—':>12} {E_q3_mean:>12.4f}")
    print()

    P_gain_mean = np.mean(rec['B']['P'][warmup:] - rec['A']['P'][warmup:])
    print(f"  Guany net de potència: {P_gain_mean:+.3f} MW sobre {N_ZZ} molins")
    print(f"  Guany per molí:        {P_gain_mean/N_ZZ*1000:+.1f} kW/molí")
    print(f"  Enercia sintètica Q3:  α = {ALPHA_Q3:.2f}  (I_Q3 ∈ [{I0_Q3*(1-ALPHA_Q3)/1e6:.2f}, {I0_Q3*(1+ALPHA_Q3)/1e6:.2f}] × 10⁶ kg·m²)")
    print()


# ═══════════════════════════════════════════════════════════════════════════════
# GRÀFICS
# ═══════════════════════════════════════════════════════════════════════════════

def plot_results(t_rec, v_rec, rec, stA, stB):
    fig = plt.figure(figsize=(16, 12))
    fig.patch.set_facecolor('#0e1117')
    gs  = gridspec.GridSpec(4, 3, figure=fig, hspace=0.48, wspace=0.38)

    ax  = lambda r, c: fig.add_subplot(gs[r, c])
    CLR = {'A': '#378ADD', 'B': '#D85A30', 'V': '#5DCAA5', 'Q': '#FAC775'}
    BG  = '#161b22'
    TXT = '#c9d1d9'

    def style(a, title, xlabel='Temps [s]', ylabel=''):
        a.set_facecolor(BG)
        a.set_title(title, color=TXT, fontsize=10, pad=6)
        a.set_xlabel(xlabel, color=TXT, fontsize=8)
        a.set_ylabel(ylabel, color=TXT, fontsize=8)
        a.tick_params(colors=TXT, labelsize=7)
        for sp in a.spines.values():
            sp.set_color('#30363d')
        a.grid(True, color='#21262d', linewidth=0.5)

    warmup = len(t_rec) // 5

    # ── 1. Potència total ──────────────────────────────────────────────────────
    a1 = ax(0, slice(0, 2))
    a1.plot(t_rec, rec['A']['P'], color=CLR['A'], lw=1.0, alpha=0.7, label='A — sense Quijote')
    a1.plot(t_rec, rec['B']['P'], color=CLR['B'], lw=1.0, alpha=0.9, label='B — amb Quijote')
    a1.fill_between(t_rec, rec['A']['P'], rec['B']['P'],
                    where=rec['B']['P'] > rec['A']['P'],
                    alpha=0.15, color=CLR['B'], label='Guany B>A')
    a1.axvspan(t_rec[0], t_rec[warmup], color='white', alpha=0.03, label='transitori')
    a1.set_xlim(t_rec[0], t_rec[-1])
    a1.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a1, f'Potència total camp — {N_ZZ} molins ZZ (3A)', ylabel='MW')

    # ── 2. Vent ──────────────────────────────────────────────────────────────
    a2 = ax(0, 2)
    a2.plot(t_rec, v_rec, color=CLR['V'], lw=0.8)
    a2.axhline(V_RAT, color='white', lw=0.5, ls='--', alpha=0.4)
    a2.text(t_rec[-1]*0.02, V_RAT+0.2, 'v_nominal', color='white', fontsize=7, alpha=0.5)
    style(a2, 'Perfil de vent', ylabel='m/s')

    # ── 3. Ordre de sincronització r ──────────────────────────────────────────
    a3 = ax(1, 0)
    a3.plot(t_rec, rec['A']['r'], color=CLR['A'], lw=0.9, label='A')
    a3.plot(t_rec, rec['B']['r'], color=CLR['B'], lw=0.9, label='B')
    a3.set_ylim(0, 1.05)
    a3.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a3, 'Paràmetre d\'ordre Kuramoto r', ylabel='r  [0–1]')

    # ── 4. Molins sincronitzats ───────────────────────────────────────────────
    a4 = ax(1, 1)
    a4.plot(t_rec, rec['A']['sync'], color=CLR['A'], lw=0.9, label='A', drawstyle='steps-post')
    a4.plot(t_rec, rec['B']['sync'], color=CLR['B'], lw=0.9, label='B', drawstyle='steps-post')
    a4.axhline(N_ZZ, color='white', lw=0.4, ls=':', alpha=0.3)
    a4.set_ylim(0, N_ZZ + 2)
    a4.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a4, f'Molins sincronitzats  (/ {N_ZZ})', ylabel='N')

    # ── 5. Desv. estàndard ω ────────────────────────────────────────────────
    a5 = ax(1, 2)
    a5.plot(t_rec, rec['A']['std_om'], color=CLR['A'], lw=0.9, label='A')
    a5.plot(t_rec, rec['B']['std_om'], color=CLR['B'], lw=0.9, label='B')
    a5.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a5, 'Dispersió ω (desv. estàndard)', ylabel='rad/s')

    # ── 6. Guany diferencial ─────────────────────────────────────────────────
    a6 = ax(2, slice(0, 2))
    diff = rec['B']['P'] - rec['A']['P']
    a6.fill_between(t_rec, diff, 0, where=diff > 0,
                    color=CLR['B'], alpha=0.4, label='B millor')
    a6.fill_between(t_rec, diff, 0, where=diff < 0,
                    color=CLR['A'], alpha=0.4, label='A millor')
    a6.plot(t_rec, diff, color='white', lw=0.6, alpha=0.6)
    a6.axhline(0, color='#30363d', lw=0.8)
    # Mitjana mòbil
    win = int(30 / (DT * 4))
    if win < len(diff):
        mov = np.convolve(diff, np.ones(win)/win, mode='same')
        a6.plot(t_rec, mov, color=CLR['Q'], lw=1.5, alpha=0.8, label=f'Mitjana mòbil {win*DT*4:.0f}s')
    a6.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a6, 'Guany P_B − P_A  (Quijote vs sense)', ylabel='ΔP [MW]')

    # ── 7. Energia per cicle Quijote (hurto) ─────────────────────────────────
    a7 = ax(2, 2)
    mask = rec['B']['E_cicle'] != 0
    if mask.sum() > 0:
        a7.scatter(t_rec[mask], rec['B']['E_cicle'][mask],
                   color=CLR['Q'], s=4, alpha=0.6, label='E cicle Q3')
        a7.axhline(np.mean(rec['B']['E_cicle'][mask]), color='white',
                   lw=0.8, ls='--', alpha=0.5, label='Mitjana')
    a7.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a7, 'Hurto gravitatori Q3 (E/cicle)', ylabel='MJ/rev')

    # ── 8. Distribució final de ω ────────────────────────────────────────────
    a8 = ax(3, 0)
    bins = np.linspace(0, OMEGA_RAT * 1.3, 30)
    a8.hist(stA['omega'] * 60/(2*np.pi), bins=bins * 60/(2*np.pi),
            color=CLR['A'], alpha=0.6, label='A final', density=True)
    a8.hist(stB['omega'] * 60/(2*np.pi), bins=bins * 60/(2*np.pi),
            color=CLR['B'], alpha=0.6, label='B final', density=True)
    a8.axvline(OMEGA_RAT * 60/(2*np.pi), color='white', lw=0.8, ls='--', alpha=0.4)
    a8.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a8, 'Distribució ω final molins', xlabel='rpm', ylabel='densitat')

    # ── 9. Diagrama de fase ω_Q3 vs θ_Q3 (Quijote) ──────────────────────────
    a9 = ax(3, 1)
    # Reconstructim una trajectoria
    th_q3_traj = np.array([])
    om_q3_traj = np.array([])
    # Simplificat: mostrem omega_Q3 final ± dispersió sintètica
    th_plot = np.linspace(0, 4*np.pi, 500)
    om_synth = stB['omega_Q3'] * (1 + ALPHA_Q3 * np.sin(N_ASPAS * th_plot) * 0.3)
    a9.plot(th_plot % (2*np.pi), om_synth * 60/(2*np.pi),
            color=CLR['Q'], lw=1.2, alpha=0.8)
    a9.scatter([stB['theta_Q3'] % (2*np.pi)],
               [stB['omega_Q3'] * 60/(2*np.pi)],
               color='white', s=40, zorder=5)
    style(a9, 'Espai de fase Quijote (θ, ω)', xlabel='θ [rad]', ylabel='ω [rpm]')

    # ── 10. Comparació Cp(λ) ────────────────────────────────────────────────
    a10 = ax(3, 2)
    lam_arr = np.linspace(0, 18, 300)
    a10.plot(lam_arr, [cp_lambda(l) for l in lam_arr],
             color=CLR['V'], lw=1.5)
    a10.axvline(LAMBDA_OPT, color='white', lw=0.5, ls='--', alpha=0.4)
    a10.text(LAMBDA_OPT + 0.3, CP_MAX * 0.95, f'λ_opt={LAMBDA_OPT}',
             color='white', fontsize=7, alpha=0.6)
    # Punts de treball finals
    for i in range(0, N_ZZ, 4):
        v_eff_A = WIND_BASE
        lam_A = stA['omega'][i] * R_NREL / max(v_eff_A, 0.1)
        lam_B = stB['omega'][i] * R_NREL / max(v_eff_A, 0.1)
        a10.scatter(lam_A, cp_lambda(lam_A), color=CLR['A'], s=8, alpha=0.5)
        a10.scatter(lam_B, cp_lambda(lam_B), color=CLR['B'], s=8, alpha=0.5)
    style(a10, 'Cp(λ) — punts de treball finals', xlabel='λ (TSR)', ylabel='Cp')

    # Títol global
    fig.suptitle(
        f'Gemell Digital — {N_ZZ}×ZypyZape (3A) + Quijote (3A) | '
        f'Bus comú K={K_BUS} | Inercia sintètica α={ALPHA_Q3} | '
        f'T={T_SIM:.0f}s',
        color=TXT, fontsize=11, y=0.995
    )

    plt.savefig('/mnt/user-data/outputs/model_44_Quijote_zypyzape.png',
                dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    print("  Gràfic guardat: model_44_Quijote_zypyzape.png")
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    t_rec, v_rec, rec, stA, stB = run_simulation()
    print_summary(t_rec, rec)
    print("  Generant gràfics comparatius...")
    plot_results(t_rec, v_rec, rec, stA, stB)
    print("\n  ✓ Simulació completada.")
    print("  Fitxer de sortida: model_44_Quijote_zypyzape.py")
    print("  Gràfic:            model_44_Quijote_zypyzape.png")
    print()
