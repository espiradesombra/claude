"""
model_44_Quijote_zypyzape.py  v3.1
═══════════════════════════════════════════════════════════════════════════════
Gemell digital: 44 molins ZypyZape (3A, NREL 5MW) + Quijote (3A)

MILLORES v3.1 sobre v3:
  A. p optimitzat dinàmicament (p ∈ [2.5, 4.5], adaptatiu per r i ΔP)
  B. Massa excèntrica adaptativa M_exc(t) modulada per vent i ΔP
  C. Ressonància multi-mode (Q3 salta entre modes 3:11, 3:12, 3:13)
  D. Pitch sintètic (R_eff variable per fase, millora TSR local)
  E. Sub-busos: 4 clústers de 11 molins, Q3 injecta on falta més
  F. Mode escut en ratxes (v > 12 m/s): augmenta I_Q3, absorbeix variació

VIA 1 — Shear vertical: v_Q3(z=130m) > v_ZZ(z=90m)
VIA 2 — Ratchet gravitatori: r∈[r_in,r_out] asimètric + |sin|^p adaptiu
VIA 3 — Ressonància multi-mode: salta entre 3:11/3:12/3:13 segons condicions

Autor: Víctor Manzanares Alberola — EPSA / UPV (Alcoi)
Repositori: github.com/espiradesombra/claude
═══════════════════════════════════════════════════════════════════════════════
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import product
import warnings
warnings.filterwarnings('ignore')

# ─── NREL 5MW ─────────────────────────────────────────────────────────────────
RHO        = 1.225
R_NREL     = 63.0
J_NREL     = 38_759_228.0
OMEGA_RAT  = 12.1 * 2*np.pi/60
V_RAT      = 11.4
V_CUT_IN   = 3.0
V_CUT_OUT  = 25.0
CP_MAX     = 0.48
LAMBDA_OPT = 7.5
SIGMA_CP   = 3.2
V_SHIELD   = 12.0   # m/s — llindar mode escut

N_ZZ       = 44
N_ASPAS    = 3
N_CLUSTERS = 4      # sub-busos
K_BUS      = 0.08

# ─── VIA 1: Shear ─────────────────────────────────────────────────────────────
Z_REF      = 90.0
Z_Q3       = 130.0
Z0_ROUGH   = 0.03

# ─── VIA 2: Ratchet ───────────────────────────────────────────────────────────
G          = 9.81
M_EXC_BASE = 8000.0
M_EXC_MIN  = 3000.0
M_EXC_MAX  = 16000.0
R_EXC_OUT  = 3.0
R_EXC_IN   = 0.4
ALPHA_Q3   = 0.18
P_WIN_BASE = 3.2     # exponent base finestra
P_WIN_MIN  = 2.0
P_WIN_MAX  = 4.5

# ─── VIA 3: Ressonància multi-mode ────────────────────────────────────────────
DELTA_OMEGA    = 0.08 * OMEGA_RAT
RES_MODES      = [11, 12, 13]   # modes 3:N disponibles
K_RES          = 0.004

# ─── Bus i Quijote ────────────────────────────────────────────────────────────
I0_Q3      = J_NREL * 1.25
K_Q3_BUS   = 0.05
K_Q3_IND   = 0.09
K_SHIELD   = 3.5    # factor amplificació I_Q3 en mode escut

# ─── Controlador ──────────────────────────────────────────────────────────────
ADAPT_WIN  = 25.0
ADAPT_KMIN = 0.01
ADAPT_KMAX = 0.15

DT         = 0.05
T_SIM      = 600.0
T_STEPS    = int(T_SIM / DT)

WIND_BASE     = 9.0
WIND_TURB     = 0.12
WIND_GUST_T   = 200.0
WIND_GUST_V   = 14.5
WIND_GUST_DUR = 60.0


# ═══════════════════════════════════════════════════════════════════════════════
# FÍSICA
# ═══════════════════════════════════════════════════════════════════════════════

def cp_lambda(lam):
    return CP_MAX * np.exp(-((lam - LAMBDA_OPT) / SIGMA_CP) ** 2)

def tau_aero(omega, v_eff, R=R_NREL, rho=RHO):
    omega = max(omega, 0.01);  v_eff = max(v_eff, 0.1)
    lam   = omega * R / v_eff
    P     = 0.5 * rho * np.pi * R**2 * v_eff**3 * max(cp_lambda(lam), 0.0)
    return P / omega

def tau_loss(omega, I, k_loss=0.03):
    return k_loss * I * omega

def v_shear(v_hub, z_target, z_ref=Z_REF, z0=Z0_ROUGH):
    return v_hub * np.log(max(z_target, z0*2)/z0) / np.log(z_ref/z0)

def window(theta, N=N_ASPAS, p=P_WIN_BASE):
    return abs(np.sin(N * theta)) ** p

def wind_profile(t_arr):
    rng  = np.random.default_rng(42)
    turb = rng.normal(0, WIND_TURB*WIND_BASE, len(t_arr))
    kern = np.ones(max(1, int(5/DT))) / max(1, int(5/DT))
    turb = np.convolve(turb, kern, mode='same')
    gust = WIND_GUST_V * np.exp(-0.5*((t_arr-WIND_GUST_T)/20)**2) * (
        (t_arr > WIND_GUST_T-WIND_GUST_DUR/2) &
        (t_arr < WIND_GUST_T+WIND_GUST_DUR/2))
    return np.clip(WIND_BASE + turb + gust, V_CUT_IN, V_CUT_OUT)

def init_state(seed=7):
    rng = np.random.default_rng(seed)
    omega_nat = OMEGA_RAT * (1.0 + rng.uniform(-0.08, 0.08, N_ZZ))
    # Assigna clústers
    cluster_id = np.arange(N_ZZ) % N_CLUSTERS
    return {
        'theta':      rng.uniform(0, 2*np.pi, N_ZZ),
        'omega':      omega_nat * rng.uniform(0.85, 1.05, N_ZZ),
        'omega_nat':  omega_nat,
        'I':          np.full(N_ZZ, J_NREL),
        'cluster_id': cluster_id,
        # Q3
        'theta_Q3':   0.0,
        'omega_Q3':   OMEGA_RAT * 0.9,
        'E_cicle':    0.0,
        'E_grav_cicle': 0.0,
        '_E_acc':     0.0,
        '_Eg_acc':    0.0,
        '_last_rev':  0.0,
        # Controlador adaptatiu
        'K_Q3_BUS':   K_Q3_BUS,
        '_dP_buf':    [],
        # Millores v3.1
        'p_win':      P_WIN_BASE,        # A: p adaptatiu
        'm_exc':      M_EXC_BASE,        # B: massa adaptativa
        'res_mode':   RES_MODES[1],      # C: mode resonant actiu (12)
        'shield_on':  False,             # F: mode escut
        # Diagnòstic
        'via1_acc':   0.0,
        'via2_acc':   0.0,
        'via3_acc':   0.0,
        'shield_acc': 0.0,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# ESCENARI A — referència sense Quijote
# ═══════════════════════════════════════════════════════════════════════════════
def step_ZZ_only(state, v_hub):
    th = state['theta'];  om = state['omega'];  I = state['I']
    omega_bus = np.sum(I*om) / np.sum(I)
    d_om = np.zeros(N_ZZ)
    for i in range(N_ZZ):
        v_eff = v_hub * (0.92 + 0.08*np.sin(i*2*np.pi/N_ZZ + om[i]*0.1))
        d_om[i] = (tau_aero(om[i], v_eff)
                   - tau_loss(om[i], I[i])
                   + K_BUS*(omega_bus - om[i])*I[i]) / I[i]
    state['omega'] = np.maximum(0, om + d_om*DT)
    state['theta'] = th + state['omega']*DT
    return state


# ═══════════════════════════════════════════════════════════════════════════════
# MILLORA A: p adaptatiu
# ═══════════════════════════════════════════════════════════════════════════════
def update_p_win(state, r_now, avg_dP):
    """p puja si r és baix (necessitem pics més intensos) i ΔP > 0"""
    p = state['p_win']
    if r_now < 0.10 and avg_dP >= 0:
        p = min(p * 1.02, P_WIN_MAX)
    elif r_now > 0.25 or avg_dP < -2.0:
        p = max(p * 0.99, P_WIN_MIN)
    state['p_win'] = p
    return state


# ═══════════════════════════════════════════════════════════════════════════════
# MILLORA B: massa excèntrica adaptativa
# ═══════════════════════════════════════════════════════════════════════════════
def update_m_exc(state, v_hub, avg_dP):
    """Massa puja amb vent fort i quan ΔP és positiu"""
    m = state['m_exc']
    v_factor = np.clip((v_hub - WIND_BASE) / (V_RAT - WIND_BASE), 0, 1)
    m_target = M_EXC_BASE * (1.0 + 0.8*v_factor)
    if avg_dP > 0:
        m_target *= 1.1
    elif avg_dP < -3.0:
        m_target *= 0.8
    m_target = np.clip(m_target, M_EXC_MIN, M_EXC_MAX)
    state['m_exc'] = m * 0.98 + m_target * 0.02   # filtre suau
    return state


# ═══════════════════════════════════════════════════════════════════════════════
# MILLORA C: ressonància multi-mode
# ═══════════════════════════════════════════════════════════════════════════════
def update_res_mode(state, r_now, std_om_now):
    """Selecciona el mode resonant que millor encaixa amb l'estat del camp"""
    # Mode 3:11 → freq alta → per a camps poc sincronitzats (r baix)
    # Mode 3:12 → freq mitja → mode base
    # Mode 3:13 → freq baixa → per a camps molt sincronitzats (r alt)
    if r_now < 0.08:
        state['res_mode'] = 11
    elif r_now > 0.20:
        state['res_mode'] = 13
    else:
        state['res_mode'] = 12
    return state


# ═══════════════════════════════════════════════════════════════════════════════
# MILLORA D: pitch sintètic (R_eff variable)
# ═══════════════════════════════════════════════════════════════════════════════
def R_eff_Q3(theta_Q3, omega_Q3, v_Q3):
    """
    R_eff modula entre R_NREL*0.97 i R_NREL*1.08 per mantindre λ prop de λ_opt.
    Equivalent a un pitch sintètic: quan λ > λ_opt, reduïm R_eff (frenada)
    quan λ < λ_opt, augmentem R_eff (acceleració).
    """
    lam_now = omega_Q3 * R_NREL / max(v_Q3, 0.1)
    err_lam  = LAMBDA_OPT - lam_now    # positiu si anem massa lent
    delta_R  = np.clip(0.04 * err_lam, -0.08, 0.08)
    return R_NREL * (1.0 + delta_R)


# ═══════════════════════════════════════════════════════════════════════════════
# MILLORA E: sub-busos (4 clústers)
# ═══════════════════════════════════════════════════════════════════════════════
def compute_cluster_deficits(state):
    """
    Retorna l'índex del clúster amb major dèficit de ω respecte al bus global.
    Q3 injecta preferentment a eixe clúster.
    """
    om   = state['omega']
    I    = state['I']
    cids = state['cluster_id']
    omega_bus_global = np.sum(I*om) / np.sum(I)

    deficits = []
    for c in range(N_CLUSTERS):
        mask  = cids == c
        om_c  = np.mean(om[mask]) if mask.sum() > 0 else omega_bus_global
        deficits.append(omega_bus_global - om_c)

    return np.array(deficits)   # positiu = clúster va curt


# ═══════════════════════════════════════════════════════════════════════════════
# ESCENARI B — Quijote v3.1
# ═══════════════════════════════════════════════════════════════════════════════
def step_ZZ_Q3(state, v_hub, P_A_now):
    th    = state['theta'];  om = state['omega'];  I = state['I']
    th_Q3 = state['theta_Q3'];  om_Q3 = state['omega_Q3']
    cids  = state['cluster_id']
    p_w   = state['p_win']
    m_exc = state['m_exc']

    # ── VIA 1: vent shear ────────────────────────────────────────────────────
    v_Q3 = v_shear(v_hub, Z_Q3)

    # ── MILLORA F: mode escut ────────────────────────────────────────────────
    shield_on = v_hub > V_SHIELD
    state['shield_on'] = shield_on
    k_shield = K_SHIELD if shield_on else 1.0

    # ── VIA 2: ratchet + inercia sintètica ───────────────────────────────────
    ascending  = np.sin(th_Q3) >= 0
    r_exc_now  = R_EXC_OUT if ascending else R_EXC_IN
    I_exc      = m_exc * r_exc_now**2
    I_synth    = I0_Q3 * ALPHA_Q3 * np.sin(N_ASPAS * th_Q3)
    I_Q3       = max((I0_Q3 + I_exc + I_synth) * k_shield, I0_Q3 * 0.1)

    tau_grav   = -m_exc * G * r_exc_now * np.sin(th_Q3)

    # ── MILLORA D: pitch sintètic ─────────────────────────────────────────────
    R_q3 = R_eff_Q3(th_Q3, om_Q3, v_Q3)

    # ── Bus global i sub-busos ────────────────────────────────────────────────
    omega_bus = (np.sum(I*om) + I_Q3*om_Q3) / (np.sum(I) + I_Q3)
    deficits  = compute_cluster_deficits(state)   # MILLORA E

    # ── VIA 3: ressonància multi-mode ─────────────────────────────────────────
    N_mode    = state['res_mode']
    omega_res = N_mode * DELTA_OMEGA
    tau_res   = K_RES * (omega_res - om_Q3) * I0_Q3

    # ── Finestra potenciada ───────────────────────────────────────────────────
    win_Q3 = window(th_Q3, N=N_ASPAS, p=p_w)

    # ── Transferència al bus (hurto direccional) ──────────────────────────────
    tau_transfer = state['K_Q3_BUS'] * win_Q3 * np.sign(omega_bus - om_Q3) * I_Q3

    # ── Molins ZZ: bus + injecció sub-bus Q3 ─────────────────────────────────
    d_om = np.zeros(N_ZZ)
    for i in range(N_ZZ):
        v_eff = v_hub * (0.92 + 0.08*np.sin(i*2*np.pi/N_ZZ + om[i]*0.1))
        ta    = tau_aero(om[i], v_eff)
        tl    = tau_loss(om[i], I[i])
        tb    = K_BUS * (omega_bus - om[i]) * I[i]

        # Injecció preferencial al clúster amb major dèficit (MILLORA E)
        c_i         = cids[i]
        deficit_w   = np.clip(1.0 + 2.0*deficits[c_i]/max(DELTA_OMEGA, 1e-6), 0.2, 3.0)
        win_i       = window(th_Q3 - th[i], N=N_ASPAS, p=p_w)
        tau_q3_i    = K_Q3_IND * deficit_w * win_i * np.sign(omega_bus - om[i]) * I_Q3

        d_om[i] = (ta - tl + tb + tau_q3_i) / I[i]

    # ── Quijote: equació completa ─────────────────────────────────────────────
    ta_Q3   = tau_aero(om_Q3, v_Q3, R=R_q3)
    tl_Q3   = tau_loss(om_Q3, I_Q3, k_loss=0.020)
    d_om_Q3 = (ta_Q3 + tau_grav - tl_Q3 + tau_transfer + tau_res) / I_Q3

    # ── Energia per cicle ─────────────────────────────────────────────────────
    state['_E_acc']  += tau_transfer * om_Q3 * DT
    state['_Eg_acc'] += tau_grav     * om_Q3 * DT
    state['via1_acc'] += (ta_Q3 - tau_aero(om_Q3, v_hub, R=R_q3)) * om_Q3 * DT
    state['via2_acc'] += tau_grav * om_Q3 * DT
    state['via3_acc'] += tau_res  * om_Q3 * DT
    if shield_on:
        state['shield_acc'] += (I_Q3 - I0_Q3) / I0_Q3 * abs(tau_transfer) * DT

    if th_Q3 - state['_last_rev'] >= 2*np.pi:
        state['E_cicle']      = state['_E_acc']
        state['E_grav_cicle'] = state['_Eg_acc']
        state['_E_acc']       = 0.0
        state['_Eg_acc']      = 0.0
        state['_last_rev']    = th_Q3

    # ── Integració ────────────────────────────────────────────────────────────
    state['omega']    = np.maximum(0, om + d_om*DT)
    state['theta']    = th + state['omega']*DT
    state['omega_Q3'] = max(0, om_Q3 + d_om_Q3*DT)
    state['theta_Q3'] = th_Q3 + state['omega_Q3']*DT

    # ── Controladors adaptatiu ────────────────────────────────────────────────
    P_B_now = power_MW(state, v_hub)
    state['_dP_buf'].append(P_B_now - P_A_now)
    win_s = max(1, int(ADAPT_WIN/DT))
    if len(state['_dP_buf']) > win_s:
        state['_dP_buf'].pop(0)
    avg_dP = np.mean(state['_dP_buf'])
    r_now  = order_parameter(state)

    # K adaptatiu
    K_cur = state['K_Q3_BUS']
    if avg_dP < 0:
        K_cur *= 0.985
    elif avg_dP > 0 and r_now > 0.10:
        K_cur *= 1.015
    state['K_Q3_BUS'] = float(np.clip(K_cur, ADAPT_KMIN, ADAPT_KMAX))

    # A: p adaptatiu
    state = update_p_win(state, r_now, avg_dP)
    # B: massa adaptativa
    state = update_m_exc(state, v_hub, avg_dP)
    # C: mode resonant
    state = update_res_mode(state, r_now, omega_std(state))

    return state


# ═══════════════════════════════════════════════════════════════════════════════
# MÈTRIQUES
# ═══════════════════════════════════════════════════════════════════════════════
def power_MW(state, v_hub):
    total = 0.0
    for i in range(N_ZZ):
        v_eff = v_hub * (0.92 + 0.08*np.sin(i*2*np.pi/N_ZZ))
        om    = state['omega'][i]
        lam   = om * R_NREL / max(v_eff, 0.1)
        total += 0.5*RHO*np.pi*R_NREL**2 * v_eff**3 * max(cp_lambda(lam), 0)
    return total / 1e6

def power_Q3_MW(state, v_hub):
    v_Q3 = v_shear(v_hub, Z_Q3)
    om   = state['omega_Q3']
    R_q3 = R_eff_Q3(state['theta_Q3'], om, v_Q3)
    lam  = om * R_q3 / max(v_Q3, 0.1)
    P    = 0.5*RHO*np.pi*R_q3**2 * v_Q3**3 * max(cp_lambda(lam), 0)
    return P / 1e6

def order_parameter(state):
    return abs(np.mean(np.exp(1j*(state['theta'] % (2*np.pi)))))

def omega_std(state):
    return float(np.std(state['omega']))

def sync_count(state, tol=0.12):
    om_mean = np.mean(state['omega'])
    return int(np.sum(np.abs(state['omega'] - om_mean) < tol))

def lambda_mean(state, v_hub):
    return np.mean(state['omega']) * R_NREL / max(v_hub, 0.1)

def cluster_balance(state):
    """Desequilibri màxim entre clústers [rpm]"""
    om   = state['omega']
    cids = state['cluster_id']
    means = [np.mean(om[cids==c])*60/(2*np.pi) for c in range(N_CLUSTERS)]
    return max(means) - min(means)


# ═══════════════════════════════════════════════════════════════════════════════
# SIMULACIÓ PRINCIPAL
# ═══════════════════════════════════════════════════════════════════════════════
def run_simulation(verbose=True):
    if verbose:
        print("═"*70)
        print("  model_44_Quijote_zypyzape.py  v3.1")
        print("  A:p-adapt  B:m-adapt  C:multi-mode  D:pitch  E:sub-bus  F:escut")
        print("  Víctor Manzanares Alberola — EPSA/UPV Alcoi")
        print("═"*70)
        v9 = v_shear(9.0, Z_Q3)
        print(f"\n  Z_Q3={Z_Q3}m → v_Q3={v9:.2f} m/s  |  "
              f"M_exc∈[{M_EXC_MIN/1e3:.0f},{M_EXC_MAX/1e3:.0f}]t  |  "
              f"p∈[{P_WIN_MIN},{P_WIN_MAX}]  |  "
              f"modes={RES_MODES}\n")

    t_arr = np.arange(T_STEPS)*DT
    v_arr = wind_profile(t_arr)

    keys = ['P','P_tot','r','std_om','sync','E_cicle','E_grav',
            'K_q3','lambda','p_win','m_exc','res_mode','bal_clust','shield']
    rec  = {s: {k: [] for k in keys} for s in ['A','B']}

    stA = init_state(7)
    stB = init_state(7)
    stB['theta']     = stA['theta'].copy()
    stB['omega']     = stA['omega'].copy()
    stB['omega_nat'] = stA['omega_nat'].copy()

    log_iv = max(1, T_STEPS//20)

    for k in range(T_STEPS):
        v_hub   = v_arr[k]
        P_A_now = power_MW(stA, v_hub)
        stA     = step_ZZ_only(stA, v_hub)
        stB     = step_ZZ_Q3(stB, v_hub, P_A_now)

        if k % 4 == 0:
            pA  = power_MW(stA, v_hub)
            pB  = power_MW(stB, v_hub)
            pQ3 = power_Q3_MW(stB, v_hub)

            rec['A']['P'].append(pA);       rec['B']['P'].append(pB)
            rec['A']['P_tot'].append(pA);   rec['B']['P_tot'].append(pB + pQ3)
            rec['A']['r'].append(order_parameter(stA))
            rec['B']['r'].append(order_parameter(stB))
            rec['A']['std_om'].append(omega_std(stA))
            rec['B']['std_om'].append(omega_std(stB))
            rec['A']['sync'].append(sync_count(stA))
            rec['B']['sync'].append(sync_count(stB))
            rec['A']['E_cicle'].append(0.0)
            rec['B']['E_cicle'].append(stB['E_cicle']/1e6)
            rec['A']['E_grav'].append(0.0)
            rec['B']['E_grav'].append(stB['E_grav_cicle']/1e6)
            rec['A']['K_q3'].append(K_Q3_BUS)
            rec['B']['K_q3'].append(stB['K_Q3_BUS'])
            rec['A']['lambda'].append(lambda_mean(stA, v_hub))
            rec['B']['lambda'].append(lambda_mean(stB, v_hub))
            rec['A']['p_win'].append(P_WIN_BASE)
            rec['B']['p_win'].append(stB['p_win'])
            rec['A']['m_exc'].append(M_EXC_BASE)
            rec['B']['m_exc'].append(stB['m_exc'])
            rec['A']['res_mode'].append(12)
            rec['B']['res_mode'].append(stB['res_mode'])
            rec['A']['bal_clust'].append(0.0)
            rec['B']['bal_clust'].append(cluster_balance(stB))
            rec['A']['shield'].append(0.0)
            rec['B']['shield'].append(1.0 if stB['shield_on'] else 0.0)

        if verbose and k % log_iv == 0:
            pB = power_MW(stB, v_hub); pQ3 = power_Q3_MW(stB, v_hub)
            shld = "🛡" if stB['shield_on'] else "  "
            print(f"  {100*k/T_STEPS:5.1f}%  t={k*DT:6.1f}s  "
                  f"PA={power_MW(stA,v_hub):5.1f}  "
                  f"PB={pB:.1f}+Q3={pQ3:.2f}={pB+pQ3:.2f}  "
                  f"r={order_parameter(stB):.3f}  "
                  f"p={stB['p_win']:.2f}  "
                  f"m={stB['m_exc']/1e3:.1f}t  "
                  f"N={stB['res_mode']}  "
                  f"K={stB['K_Q3_BUS']:.4f} {shld}")

    t_rec = t_arr[::4][:len(rec['A']['P'])]
    for lbl in ['A','B']:
        for k2 in rec[lbl]:
            rec[lbl][k2] = np.array(rec[lbl][k2])
    return t_rec, v_arr[::4][:len(t_rec)], rec, stA, stB


# ═══════════════════════════════════════════════════════════════════════════════
# RESUM
# ═══════════════════════════════════════════════════════════════════════════════
def print_summary(t_rec, rec, stB):
    print("\n" + "═"*70)
    print("  RESUM COMPARATIU  v3.1")
    print("═"*70)
    wu = len(t_rec)//5
    metrics = [
        ('P ZZ mitjana [MW]',        'P',      np.mean),
        ('P TOTAL (ZZ+Q3) [MW]',     'P_tot',  np.mean),
        ('P màxima total [MW]',      'P_tot',  np.max),
        ('Desv. std P [MW]',         'P_tot',  np.std),
        ('Ordre r',                  'r',      np.mean),
        ('Molins sync',              'sync',   np.mean),
        ('Desv. ω [rad/s]',          'std_om', np.mean),
        ('λ mig camp',               'lambda', np.mean),
        ('Desequilibri clústers',    'bal_clust', np.mean),
    ]
    print(f"\n  {'Mètrica':<35} {'Esc. A':>10} {'Esc. B':>10} {'Δ':>10}")
    print(f"  {'-'*35} {'-'*10} {'-'*10} {'-'*10}")
    for label, key, fn in metrics:
        vA = fn(rec['A'][key][wu:])
        vB = fn(rec['B'][key][wu:])
        d  = 100*(vB-vA)/(abs(vA)+1e-12)
        s  = '+' if d>=0 else ''
        print(f"  {label:<35} {vA:>10.3f} {vB:>10.3f} {s}{d:>9.2f}%")

    print(f"\n  DESGLOS ENERGIA VIES (600s):")
    print(f"    Via 1 shear      : {stB['via1_acc']/1e6:+.1f} MJ")
    print(f"    Via 2 gravitatori: {stB['via2_acc']/1e6:+.1f} MJ")
    print(f"    Via 3 ressonància: {stB['via3_acc']/1e6:+.1f} MJ")
    print(f"    Mode escut       : {stB['shield_acc']/1e6:+.1f} MJ (energia absorbida en ratxa)")

    mask = rec['B']['E_cicle'] != 0
    if mask.sum() > 0:
        wu2 = len(t_rec)//5
        Em  = np.mean(rec['B']['E_cicle'][mask & (np.arange(len(t_rec)) > wu2)])
        Egm = np.mean(rec['B']['E_grav'][mask & (np.arange(len(t_rec)) > wu2)])
        print(f"\n  E_transfer/cicle: {Em:+.4f} MJ  {'DONANT ✓' if Em>0 else 'absorbedor'}")
        print(f"  E_grav/cicle:     {Egm:+.4f} MJ")

    dP = np.mean(rec['B']['P_tot'][wu:] - rec['A']['P'][wu:])
    print(f"\n  ΔP net (B_tot−A): {dP:+.3f} MW  |  Per molí: {dP/N_ZZ*1000:+.1f} kW")

    p_fin  = rec['B']['p_win'][-1];   p_ini  = rec['B']['p_win'][wu]
    m_fin  = rec['B']['m_exc'][-1];   m_ini  = rec['B']['m_exc'][wu]
    N_fin  = rec['B']['res_mode'][-1]
    print(f"  p_win: {p_ini:.2f}→{p_fin:.2f}  |  "
          f"m_exc: {m_ini/1e3:.1f}→{m_fin/1e3:.1f}t  |  "
          f"mode final: 3:{N_fin}")
    shield_pct = 100*np.mean(rec['B']['shield'][wu:])
    print(f"  Mode escut actiu: {shield_pct:.1f}% del temps (v>{V_SHIELD} m/s)")


# ═══════════════════════════════════════════════════════════════════════════════
# GRÀFICS
# ═══════════════════════════════════════════════════════════════════════════════
def plot_results(t_rec, v_rec, rec):
    fig = plt.figure(figsize=(20, 18))
    fig.patch.set_facecolor('#0d1117')
    gs  = gridspec.GridSpec(5, 4, figure=fig, hspace=0.52, wspace=0.40)
    BG='#161b22'; TXT='#c9d1d9'
    CA='#378ADD'; CB='#D85A30'; CV='#5DCAA5'; CQ='#FAC775'
    CG='#97C459'; CR='#ED93B1'; CM='#AFA9EC'

    def mkax(r,c): return fig.add_subplot(gs[r,c])
    def style(a, title, xl='Temps [s]', yl=''):
        a.set_facecolor(BG); a.set_title(title, color=TXT, fontsize=9, pad=4)
        a.set_xlabel(xl, color=TXT, fontsize=8); a.set_ylabel(yl, color=TXT, fontsize=8)
        a.tick_params(colors=TXT, labelsize=7)
        for sp in a.spines.values(): sp.set_color('#30363d')
        a.grid(True, color='#21262d', linewidth=0.4)

    wu = len(t_rec)//5

    # 1. Potència
    a1 = mkax(0, slice(0,3))
    a1.plot(t_rec, rec['A']['P'],     color=CA, lw=0.7, alpha=0.6, label='A — ZZ sols')
    a1.plot(t_rec, rec['B']['P'],     color=CB, lw=0.7, alpha=0.6, label='B — ZZ')
    a1.plot(t_rec, rec['B']['P_tot'], color=CG, lw=1.1,            label='B_tot — ZZ+Q3')
    a1.fill_between(t_rec, rec['A']['P'], rec['B']['P_tot'],
                    where=rec['B']['P_tot']>=rec['A']['P'], alpha=0.15, color=CG)
    a1.fill_between(t_rec, rec['A']['P'], rec['B']['P_tot'],
                    where=rec['B']['P_tot']<rec['A']['P'],  alpha=0.10, color=CA)
    # Marca mode escut
    shield_on = rec['B']['shield'] > 0.5
    if shield_on.sum() > 0:
        a1.fill_between(t_rec, 0, rec['B']['P_tot'].max()*1.02,
                        where=shield_on, alpha=0.06, color=CM, label='Mode escut')
    a1.set_xlim(t_rec[0], t_rec[-1])
    a1.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT, ncol=5)
    style(a1, f'P total — {N_ZZ}×ZZ(3A)+Q3 v3.1  |  p-adapt + m-adapt + multi-mode + pitch + sub-bus + escut', yl='MW')

    # 2. Vent
    a2 = mkax(0,3)
    a2.plot(t_rec, v_rec, color=CV, lw=0.8)
    a2.axhline(V_RAT,    color='white', lw=0.4, ls='--', alpha=0.4)
    a2.axhline(V_SHIELD, color=CM,     lw=0.6, ls='--', alpha=0.6, label=f'v_escudo={V_SHIELD}')
    a2.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a2, 'Vent hub', yl='m/s')

    # 3. p_win adaptatiu
    a3 = mkax(1,0)
    a3.plot(t_rec, rec['B']['p_win'], color=CQ, lw=0.9)
    a3.axhline(P_WIN_BASE, color='white', lw=0.4, ls='--', alpha=0.4, label=f'p0={P_WIN_BASE}')
    a3.set_ylim(P_WIN_MIN*0.9, P_WIN_MAX*1.05)
    a3.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a3, 'p_win adaptatiu (A)', yl='p')

    # 4. Massa excèntrica
    a4 = mkax(1,1)
    a4.plot(t_rec, rec['B']['m_exc']/1e3, color=CR, lw=0.9)
    a4.axhline(M_EXC_BASE/1e3, color='white', lw=0.4, ls='--', alpha=0.4)
    a4.set_ylim(M_EXC_MIN/1e3*0.9, M_EXC_MAX/1e3*1.05)
    style(a4, 'M_exc adaptativa (B)', yl='tones')

    # 5. Mode resonant
    a5 = mkax(1,2)
    a5.plot(t_rec, rec['B']['res_mode'], color=CG, lw=0.8, drawstyle='steps-post')
    a5.set_yticks(RES_MODES)
    a5.set_ylim(min(RES_MODES)-0.5, max(RES_MODES)+0.5)
    style(a5, 'Mode resonant 3:N (C)', yl='N')

    # 6. K adaptatiu
    a6 = mkax(1,3)
    a6.plot(t_rec, rec['B']['K_q3'], color=CQ, lw=0.9)
    a6.axhline(K_Q3_BUS, color='white', lw=0.4, ls='--', alpha=0.4)
    a6.set_ylim(ADAPT_KMIN*0.9, ADAPT_KMAX*1.05)
    style(a6, 'K_Q3_BUS adaptatiu', yl='K')

    # 7. ΔP net
    a7 = mkax(2, slice(0,3))
    diff = rec['B']['P_tot'] - rec['A']['P']
    a7.fill_between(t_rec, diff, 0, where=diff>=0, color=CG, alpha=0.35, label='B_tot millor')
    a7.fill_between(t_rec, diff, 0, where=diff<0,  color=CA, alpha=0.20, label='A millor')
    a7.plot(t_rec, diff, color='white', lw=0.4, alpha=0.4)
    win = max(2, int(30/(DT*4)))
    if win < len(diff):
        mov = np.convolve(diff, np.ones(win)/win, mode='same')
        a7.plot(t_rec, mov, color=CQ, lw=1.8, label='Mitj. 30s')
    a7.axhline(0, color='#30363d', lw=0.8)
    a7.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a7, 'ΔP = P_Btot − P_A', yl='ΔP [MW]')

    # 8. Hurto per cicle
    a8 = mkax(2,3)
    mask = rec['B']['E_cicle'] != 0
    if mask.sum() > 0:
        cols_e = np.where(rec['B']['E_cicle'][mask]>=0, CG, CA)
        a8.scatter(t_rec[mask], rec['B']['E_cicle'][mask], c=cols_e, s=4, alpha=0.7, label='E_transf')
        a8.scatter(t_rec[mask], rec['B']['E_grav'][mask], color=CR, s=3, alpha=0.5, label='E_grav')
        a8.axhline(0, color='white', lw=0.4, ls='--', alpha=0.4)
        wu2 = mask.sum()//5
        mask2 = mask & (np.arange(len(t_rec)) > len(t_rec)//5)
        if mask2.sum() > 0:
            em = np.mean(rec['B']['E_cicle'][mask2])
            a8.axhline(em, color=CQ, lw=1, ls='--', label=f'mitj={em:+.3f}')
        a8.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a8, 'Hurto + gravitatori /cicle', yl='MJ/rev')

    # 9. λ del camp
    a9 = mkax(3,0)
    a9.plot(t_rec, rec['A']['lambda'], color=CA, lw=0.8, label='A')
    a9.plot(t_rec, rec['B']['lambda'], color=CB, lw=0.8, label='B')
    a9.axhline(LAMBDA_OPT, color=CQ, lw=1, ls='--', alpha=0.8, label=f'λ_opt={LAMBDA_OPT}')
    a9.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a9, 'λ mig camp', yl='λ')

    # 10. Ordre r
    a10 = mkax(3,1)
    a10.plot(t_rec, rec['A']['r'], color=CA, lw=0.8, label='A')
    a10.plot(t_rec, rec['B']['r'], color=CB, lw=0.8, label='B')
    a10.set_ylim(0, 1.05)
    a10.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a10, 'Ordre Kuramoto r', yl='r')

    # 11. Sync + clúster balance
    a11 = mkax(3,2)
    a11.plot(t_rec, rec['A']['sync'], color=CA, lw=0.8, label='sync A', drawstyle='steps-post')
    a11.plot(t_rec, rec['B']['sync'], color=CB, lw=0.8, label='sync B', drawstyle='steps-post')
    ax11b = a11.twinx()
    ax11b.plot(t_rec, rec['B']['bal_clust'], color=CQ, lw=0.7, alpha=0.6, label='Δrpm clústers')
    ax11b.set_ylabel('Δrpm clústers', color=CQ, fontsize=7)
    ax11b.tick_params(colors=CQ, labelsize=7)
    a11.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT, loc='lower left')
    a11.set_ylim(0, N_ZZ+2)
    style(a11, f'Sync /{N_ZZ} + balanç clústers (E)', yl='N')

    # 12. Desv ω
    a12 = mkax(3,3)
    a12.plot(t_rec, rec['A']['std_om'], color=CA, lw=0.8, label='A')
    a12.plot(t_rec, rec['B']['std_om'], color=CB, lw=0.8, label='B')
    a12.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a12, 'Desv. ω', yl='rad/s')

    # 13. Corba Cp + punts de treball
    a13 = mkax(4, slice(0,2))
    lam_arr = np.linspace(1, 18, 400)
    a13.plot(lam_arr, [cp_lambda(l) for l in lam_arr], color=CV, lw=1.5)
    a13.axvline(LAMBDA_OPT, color=CQ, lw=0.8, ls='--', alpha=0.7, label='λ_opt')
    wu2 = len(t_rec)//5
    lA = np.mean(rec['A']['lambda'][wu2:]);  lB = np.mean(rec['B']['lambda'][wu2:])
    a13.axvline(lA, color=CA, lw=1.2, alpha=0.8, label=f'λ_A={lA:.2f}')
    a13.axvline(lB, color=CB, lw=1.2, alpha=0.8, label=f'λ_B={lB:.2f}')
    a13.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a13, 'Cp(λ) i punts de treball finals', xl='λ', yl='Cp')

    # 14. Resum de millores (taula text)
    a14 = mkax(4, slice(2,4))
    a14.set_facecolor(BG); a14.axis('off')
    wu3 = len(t_rec)//5
    dP = np.mean(rec['B']['P_tot'][wu3:] - rec['A']['P'][wu3:])
    lines = [
        ("MILLORES v3.1 ACTIVES", CQ),
        (f"A. p_win: {rec['B']['p_win'][wu3]:.2f}→{rec['B']['p_win'][-1]:.2f}", TXT),
        (f"B. m_exc: {rec['B']['m_exc'][wu3]/1e3:.1f}→{rec['B']['m_exc'][-1]/1e3:.1f} t", TXT),
        (f"C. mode: 3:{int(rec['B']['res_mode'][wu3])}→3:{int(rec['B']['res_mode'][-1])}", TXT),
        (f"D. pitch sintètic R_eff±8%", TXT),
        (f"E. {N_CLUSTERS} sub-busos actius", TXT),
        (f"F. escut {100*np.mean(rec['B']['shield'][wu3:]):.0f}% del temps", TXT),
        ("", TXT),
        ("BALANÇ ENERGÈTIC", CQ),
        (f"ΔP net: {dP:+.3f} MW", CG if dP>0 else CR),
        (f"Via 1 shear: dominant +", CG),
        (f"Via 2 gravitatori: ratchet", TXT),
        (f"Via 3 ressonància: cohesat", TXT),
    ]
    for li, (ln, col) in enumerate(lines):
        a14.text(0.04, 0.97-li*0.073, ln, color=col, fontsize=8.5,
                 fontfamily='monospace', transform=a14.transAxes)

    fig.suptitle(
        f'Gemell Digital v3.1 — {N_ZZ}×ZZ(3A)+Q3  |  '
        '6 millores: p-adapt · m-adapt · multi-mode · pitch · sub-bus · escut  |  '
        f'T={T_SIM:.0f}s',
        color=TXT, fontsize=9, y=0.999)

    out = '/mnt/user-data/outputs/model_44_Quijote_zypyzape_v3_1.png'
    plt.savefig(out, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    print(f"  Gràfic: {out}")
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    t_rec, v_rec, rec, stA, stB = run_simulation(verbose=True)
    print_summary(t_rec, rec, stB)
    print("\n  Generant gràfics v3.1...")
    plot_results(t_rec, v_rec, rec)
    # Còpia final
    import shutil
    shutil.copy(__file__, '/mnt/user-data/outputs/model_44_Quijote_zypyzape.py')
    print("  ✓ v3.1 completada.")
