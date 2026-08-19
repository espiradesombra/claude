"""
model_44_Quijote_zypyzape.py  v4
═══════════════════════════════════════════════════════════════════════════════
Gemell digital: 44 molins ZypyZape (3A, NREL 5MW) + Quijote (3A)

MILLORES v4 sobre v3.1 (basades en revisió externa):
────────────────────────────────────────────────────
  1. POLÍTICA APRESA (Policy Gradient Finit)
     Vector estat x = [p_win, m_exc_norm, K_bus, N_mode_norm, K_ind]
     Funció cost J = -<ΔP> + λ1·std(P) + λ2·bal_cluster
     Gradient finit (pertorbació ±ε) cada T_POLICY s
     Actualització: x ← x + η·∇J  (gradient ascent sobre ΔP net)

  2. RESTRICCIONS ESTRUCTURALS / CÀRREGUES
     Penalització si τ_var = var(τ_aero_Q3) > τ_var_max
     Penalització si I_Q3 > I_max_struct
     Incorporat a la funció cost J per no sobreexcitar

  3. K_Q3_IND ADAPTATIU (lligat a cluster_balance)
     K_ind = K_ind_base · (1 + γ · max(deficits) / Δω_ref)
     Si un clúster va molt curt, K_ind puja per injectar-hi més

  4. ANÀLISI ESPECTRAL (FFT de ω)
     Cada T_FFT s, FFT de l'historial de omega_bus
     Detecta freq. dominant → escull mode resonant directament
     (no heurístic per r, sinó per espectre real)

  5. TSR KEEPER EXPLÍCIT
     Objectiu: λ_camp_B ∈ [λ_min, λ_max] = [7.2, 7.8]
     PID suau sobre ω_bus target per mantindre λ prop de λ_opt
     Actua via K_bus (augmenta/redueix sincronització) i pitch Q3

Arquitectura de control v4:
  ┌─────────────────────────────────────────────────────┐
  │  FAST LOOP (dt=0.05s):  física + injecció Q3        │
  │  MID LOOP (T=25s):      K_bus adapt + shield + sub  │
  │  SLOW LOOP (T=60s):     FFT spectre + TSR PID       │
  │  POLICY LOOP (T=120s):  gradient finit → nova x     │
  └─────────────────────────────────────────────────────┘

Autor: Víctor Manzanares Alberola — EPSA / UPV (Alcoi)
Repositori: github.com/espiradesombra/claude
═══════════════════════════════════════════════════════════════════════════════
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ─── NREL 5MW ──────────────────────────────────────────────────────────────
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
V_SHIELD   = 12.0

N_ZZ       = 44
N_ASPAS    = 3
N_CLUSTERS = 4
K_BUS_BASE = 0.08

# ─── VIA 1: Shear ──────────────────────────────────────────────────────────
Z_REF      = 90.0
Z_Q3       = 130.0
Z0_ROUGH   = 0.03

# ─── VIA 2: Ratchet ────────────────────────────────────────────────────────
G          = 9.81
M_EXC_BASE = 8000.0
M_EXC_MIN  = 2000.0
M_EXC_MAX  = 18000.0
R_EXC_OUT  = 3.0
R_EXC_IN   = 0.4
ALPHA_Q3   = 0.18
P_WIN_MIN  = 1.5
P_WIN_MAX  = 4.5

# ─── VIA 3: Ressonància ────────────────────────────────────────────────────
DELTA_OMEGA = 0.08 * OMEGA_RAT
RES_MODES   = [11, 12, 13]
K_RES       = 0.004

# ─── Bus / Q3 ──────────────────────────────────────────────────────────────
I0_Q3       = J_NREL * 1.25
I_MAX_STRUCT= I0_Q3 * 6.0      # límit estructural
K_Q3_BUS_0  = 0.05
K_Q3_IND_0  = 0.09
K_SHIELD    = 3.5

# ─── TSR KEEPER (millora 5) ────────────────────────────────────────────────
LAMBDA_MIN  = 7.2
LAMBDA_MAX  = 7.8
TSR_KP      = 0.008            # guany proporcional PID
TSR_KI      = 0.0005           # guany integral
TSR_KD      = 0.002            # guany derivatiu

# ─── POLÍTICA (millora 1) ──────────────────────────────────────────────────
T_POLICY    = 120.0            # s entre actualitzacions de política
EPS_GRAD    = 0.08             # pertorbació gradient finit
ETA_POLICY  = 0.04             # taxa aprenentatge
LAM1_COST   = 0.15             # pes std(P) en cost
LAM2_COST   = 0.05             # pes bal_cluster

# ─── RESTRICCIONS CÀRREGUES (millora 2) ───────────────────────────────────
TAU_VAR_MAX = 1.5e12           # N²m² var màxima de parell Q3
K_LOAD_PEN  = 0.02             # penalització càrregues a política

# ─── LOOPS TEMPORALS ───────────────────────────────────────────────────────
ADAPT_WIN   = 25.0
T_FFT       = 60.0             # s entre anàlisis FFT
T_TSR_PID   = 1.0              # s període PID TSR
ADAPT_KMIN  = 0.01
ADAPT_KMAX  = 0.18

DT          = 0.05
T_SIM       = 600.0
T_STEPS     = int(T_SIM / DT)

WIND_BASE      = 9.0
WIND_TURB      = 0.12
WIND_GUST_T    = 200.0
WIND_GUST_V    = 14.5
WIND_GUST_DUR  = 60.0


# ═══════════════════════════════════════════════════════════════════════════
# FÍSICA BASE
# ═══════════════════════════════════════════════════════════════════════════

def cp_lambda(lam):
    return CP_MAX * np.exp(-((lam - LAMBDA_OPT) / SIGMA_CP)**2)

def tau_aero(omega, v_eff, R=R_NREL, rho=RHO):
    omega = max(omega, 0.01);  v_eff = max(v_eff, 0.1)
    lam   = omega * R / v_eff
    P     = 0.5*rho*np.pi*R**2*v_eff**3 * max(cp_lambda(lam), 0.0)
    return P / omega

def tau_loss(omega, I, k_loss=0.03):
    return k_loss * I * omega

def v_shear(v_hub, z=Z_Q3, z_ref=Z_REF, z0=Z0_ROUGH):
    return v_hub * np.log(max(z, z0*2)/z0) / np.log(z_ref/z0)

def window(theta, N=N_ASPAS, p=3.2):
    return abs(np.sin(N * theta)) ** p

def R_eff_Q3(omega_Q3, v_Q3, theta_Q3):
    """Pitch sintètic: R_eff per mantindre λ prop de λ_opt"""
    lam_now = omega_Q3 * R_NREL / max(v_Q3, 0.1)
    err     = LAMBDA_OPT - lam_now
    delta_R = np.clip(0.035 * err, -0.08, 0.08)
    return R_NREL * (1.0 + delta_R)

def wind_profile(t_arr):
    rng  = np.random.default_rng(42)
    turb = rng.normal(0, WIND_TURB*WIND_BASE, len(t_arr))
    kern = np.ones(max(1, int(5/DT))) / max(1, int(5/DT))
    turb = np.convolve(turb, kern, mode='same')
    gust = WIND_GUST_V * np.exp(-0.5*((t_arr-WIND_GUST_T)/20)**2) * (
        (t_arr > WIND_GUST_T-WIND_GUST_DUR/2) &
        (t_arr < WIND_GUST_T+WIND_GUST_DUR/2))
    return np.clip(WIND_BASE + turb + gust, V_CUT_IN, V_CUT_OUT)


# ═══════════════════════════════════════════════════════════════════════════
# ESTAT INICIAL
# ═══════════════════════════════════════════════════════════════════════════

def init_state(seed=7):
    rng = np.random.default_rng(seed)
    omega_nat = OMEGA_RAT * (1.0 + rng.uniform(-0.08, 0.08, N_ZZ))
    return {
        # Camp
        'theta':       rng.uniform(0, 2*np.pi, N_ZZ),
        'omega':       omega_nat * rng.uniform(0.85, 1.05, N_ZZ),
        'omega_nat':   omega_nat,
        'I':           np.full(N_ZZ, J_NREL),
        'cluster_id':  np.arange(N_ZZ) % N_CLUSTERS,
        # Q3
        'theta_Q3':    0.0,
        'omega_Q3':    OMEGA_RAT * 0.9,
        'E_cicle':     0.0,
        'E_grav_cicle':0.0,
        '_E_acc':      0.0,
        '_Eg_acc':     0.0,
        '_last_rev':   0.0,
        # Vies
        'via1_acc':    0.0,
        'via2_acc':    0.0,
        'via3_acc':    0.0,
        'shield_acc':  0.0,
        # Controladors
        'K_Q3_BUS':    K_Q3_BUS_0,
        'K_Q3_IND':    K_Q3_IND_0,
        'p_win':       3.2,
        'm_exc':       M_EXC_BASE,
        'res_mode':    12,
        'shield_on':   False,
        '_dP_buf':     [],
        '_P_buf':      [],
        # TSR PID (millora 5)
        '_tsr_err_int': 0.0,
        '_tsr_err_prev':0.0,
        '_K_bus_tsr':   K_BUS_BASE,
        # Política (millora 1)
        '_policy_x':   np.array([3.2,
                                  M_EXC_BASE/M_EXC_MAX,
                                  K_Q3_BUS_0/ADAPT_KMAX,
                                  12.0/max(RES_MODES),
                                  K_Q3_IND_0/0.20]),
        '_policy_J':   0.0,
        '_policy_buf': [],
        '_policy_t':   0.0,
        '_policy_phase': 0,   # 0=mesura, 1=pertorba+, 2=pertorba-
        '_policy_perturb_idx': 0,
        '_policy_grad': np.zeros(5),
        '_policy_J_plus':  np.zeros(5),
        '_policy_J_minus': np.zeros(5),
        # Espectre FFT (millora 4)
        '_fft_buf':    [],
        '_fft_t':      0.0,
        # Càrregues (millora 2)
        '_tau_buf':    [],
        'tau_var':     0.0,
        # Diagnòstic
        'J_cost':      0.0,
        'lambda_track': 0.0,
    }


# ═══════════════════════════════════════════════════════════════════════════
# ESCENARI A
# ═══════════════════════════════════════════════════════════════════════════

def step_ZZ_only(state, v_hub):
    th = state['theta'];  om = state['omega'];  I = state['I']
    omega_bus = np.sum(I*om) / np.sum(I)
    d_om = np.zeros(N_ZZ)
    for i in range(N_ZZ):
        v_eff = v_hub * (0.92 + 0.08*np.sin(i*2*np.pi/N_ZZ + om[i]*0.1))
        d_om[i] = (tau_aero(om[i], v_eff)
                   - tau_loss(om[i], I[i])
                   + K_BUS_BASE*(omega_bus - om[i])*I[i]) / I[i]
    state['omega'] = np.maximum(0, om + d_om*DT)
    state['theta'] = th + state['omega']*DT
    return state


# ═══════════════════════════════════════════════════════════════════════════
# MILLORA 4: FFT espectral
# ═══════════════════════════════════════════════════════════════════════════

def update_res_mode_fft(state):
    """Detecta freq dominant de ω_bus via FFT i escull mode resonant"""
    buf = state['_fft_buf']
    if len(buf) < 64:
        return state
    signal = np.array(buf[-128:]) if len(buf) >= 128 else np.array(buf)
    signal = signal - np.mean(signal)
    fft    = np.abs(np.fft.rfft(signal))
    freqs  = np.fft.rfftfreq(len(signal), d=DT*4)  # guardat cada 4 passos
    if len(fft) < 2:
        return state
    # Pica dominant (exclou DC)
    idx_dom = np.argmax(fft[1:]) + 1
    f_dom   = freqs[idx_dom]
    # Distància a cada mode resonant
    dists = [abs(f_dom - n*DELTA_OMEGA/(2*np.pi)) for n in RES_MODES]
    best_mode = RES_MODES[int(np.argmin(dists))]
    state['res_mode'] = best_mode
    return state


# ═══════════════════════════════════════════════════════════════════════════
# MILLORA 5: TSR KEEPER PID
# ═══════════════════════════════════════════════════════════════════════════

def tsr_keeper_pid(state, v_hub):
    """
    PID que manté λ_camp ∈ [λ_min, λ_max].
    Actua sobre K_bus (ajusta sincronització → ajusta ω_bus).
    """
    om_mean = np.mean(state['omega'])
    lam_now = om_mean * R_NREL / max(v_hub, 0.1)
    state['lambda_track'] = lam_now

    # Error: volem λ → λ_opt
    lam_target = LAMBDA_OPT
    err = lam_target - lam_now   # positiu si λ < λ_opt (anem massa lents → pujar ω)

    # PID
    state['_tsr_err_int']  += err * DT * int(T_TSR_PID/DT)
    state['_tsr_err_int']   = np.clip(state['_tsr_err_int'], -5.0, 5.0)
    derr = (err - state['_tsr_err_prev']) / (DT * int(T_TSR_PID/DT))
    state['_tsr_err_prev']  = err

    # Correcció K_bus: si λ > λ_opt augmentem K_bus (fre ω), si λ < λ_opt el baixem
    dK = -(TSR_KP*err + TSR_KI*state['_tsr_err_int'] + TSR_KD*derr)
    K_new = state['_K_bus_tsr'] + dK
    state['_K_bus_tsr'] = float(np.clip(K_new, 0.02, 0.18))
    return state


# ═══════════════════════════════════════════════════════════════════════════
# MILLORA 3: K_Q3_IND adaptatiu per cluster_balance
# ═══════════════════════════════════════════════════════════════════════════

def compute_cluster_deficits(state):
    om   = state['omega'];  cids = state['cluster_id']
    om_bus = np.mean(om)
    return np.array([np.mean(om[cids==c]) - om_bus
                     if (cids==c).sum()>0 else 0.0
                     for c in range(N_CLUSTERS)])

def update_K_ind(state):
    """K_ind puja si algun clúster té gran dèficit"""
    deficits = compute_cluster_deficits(state)
    max_def  = max(abs(deficits))
    gamma    = 2.5
    K_ind    = K_Q3_IND_0 * (1.0 + gamma * max_def / max(DELTA_OMEGA, 1e-6))
    state['K_Q3_IND'] = float(np.clip(K_ind, 0.02, 0.25))
    return state, compute_cluster_deficits(state)


# ═══════════════════════════════════════════════════════════════════════════
# MILLORA 1: POLÍTICA GRADIENT FINIT
# ═══════════════════════════════════════════════════════════════════════════

def policy_cost(dP_mean, P_std, bal_cluster, tau_var):
    """J = -ΔP + λ1·std(P) + λ2·bal + λ3·load_penalty"""
    load_pen = K_LOAD_PEN * max(0, tau_var - TAU_VAR_MAX) / TAU_VAR_MAX
    return -dP_mean + LAM1_COST*P_std + LAM2_COST*bal_cluster + load_pen

def decode_policy(x):
    """Descodifica vector x ∈ [0,1]^5 en paràmetres físics"""
    p_win   = P_WIN_MIN + x[0] * (P_WIN_MAX - P_WIN_MIN)
    m_exc   = M_EXC_MIN + x[1] * (M_EXC_MAX - M_EXC_MIN)
    K_bus   = ADAPT_KMIN + x[2] * (ADAPT_KMAX - ADAPT_KMIN)
    N_mode  = RES_MODES[int(np.clip(round(x[3] * (len(RES_MODES)-1)), 0, len(RES_MODES)-1))]
    K_ind   = 0.02 + x[4] * 0.23
    return p_win, m_exc, K_bus, N_mode, K_ind

def policy_gradient_step(state, t_now, P_A_arr, P_B_arr, bal_arr, tau_v):
    """
    Cada T_POLICY s: estima gradient finit i actualitza x.
    Fase 0 (T_POLICY/3): mesura J(x) baseline
    Fase 1 (T_POLICY/3): pertorba x[i]+ε, mesura J+
    Fase 2 (T_POLICY/3): pertorba x[i]-ε, mesura J-
    Actualització: x[i] += -η · (J+ - J-) / (2ε)
    """
    if len(P_A_arr) < 10 or len(P_B_arr) < 10:
        return state

    n = min(len(P_B_arr), len(P_A_arr))
    dP   = np.mean(np.array(P_B_arr[-n:]) - np.array(P_A_arr[-n:]))
    Pstd = np.std(P_B_arr)
    bal  = np.mean(bal_arr) if len(bal_arr)>0 else 0.0
    J    = policy_cost(dP, Pstd, bal, tau_v)

    state['_policy_J']    = J
    state['J_cost']       = J

    phase = state['_policy_phase']
    idx   = state['_policy_perturb_idx']
    x     = state['_policy_x'].copy()

    if phase == 0:
        # Guardem J baseline per a cada component
        state['_policy_J_plus'][idx]  = J  # reutilitzem com a J_base
        state['_policy_phase'] = 1
    elif phase == 1:
        state['_policy_J_plus'][idx]  = J
        state['_policy_phase'] = 2
        # Aplica pertorbació negativa
        x[idx] = np.clip(x[idx] - EPS_GRAD, 0.01, 0.99)
        state['_policy_x'] = x
    elif phase == 2:
        state['_policy_J_minus'][idx] = J
        # Gradient estimat per a idx
        grad_i = (state['_policy_J_plus'][idx] - state['_policy_J_minus'][idx]) / (2*EPS_GRAD)
        state['_policy_grad'][idx] = grad_i
        # Restaura i actualitza
        x[idx] = np.clip(x[idx] + EPS_GRAD - ETA_POLICY*grad_i, 0.01, 0.99)
        state['_policy_x'] = x
        # Passa al següent component
        next_idx = (idx + 1) % len(x)
        state['_policy_perturb_idx'] = next_idx
        state['_policy_phase'] = 0
        # Aplica pertorbació positiva del proper component
        x2 = state['_policy_x'].copy()
        x2[next_idx] = np.clip(x2[next_idx] + EPS_GRAD, 0.01, 0.99)
        state['_policy_x'] = x2

    # Decodifica nova política
    p_w, m_e, K_b, N_m, K_i = decode_policy(state['_policy_x'])
    state['p_win']    = float(p_w)
    state['m_exc']    = float(m_e)
    state['K_Q3_BUS'] = float(np.clip(K_b, ADAPT_KMIN, ADAPT_KMAX))
    state['res_mode'] = int(N_m)
    # K_ind es gestiona per cluster_balance, però la política pot donar-li un floor
    state['K_Q3_IND'] = max(state['K_Q3_IND'], float(K_i))
    return state


# ═══════════════════════════════════════════════════════════════════════════
# ESCENARI B — QUIJOTE v4
# ═══════════════════════════════════════════════════════════════════════════

def step_ZZ_Q3(state, v_hub, P_A_now, t_now):
    th    = state['theta'];  om = state['omega'];  I = state['I']
    th_Q3 = state['theta_Q3'];  om_Q3 = state['omega_Q3']
    cids  = state['cluster_id']
    p_w   = state['p_win']
    m_exc = state['m_exc']

    # ── VIA 1: shear ──────────────────────────────────────────────────────
    v_Q3 = v_shear(v_hub, Z_Q3)

    # ── SHIELD ────────────────────────────────────────────────────────────
    shield_on    = v_hub > V_SHIELD
    state['shield_on'] = shield_on
    k_shield     = K_SHIELD if shield_on else 1.0

    # ── VIA 2: ratchet ────────────────────────────────────────────────────
    ascending = np.sin(th_Q3) >= 0
    r_exc_now = R_EXC_OUT if ascending else R_EXC_IN
    I_exc     = m_exc * r_exc_now**2
    I_synth   = I0_Q3 * ALPHA_Q3 * np.sin(N_ASPAS * th_Q3)
    I_Q3      = max((I0_Q3 + I_exc + I_synth) * k_shield, I0_Q3*0.1)
    # Restricció estructural (millora 2)
    I_Q3      = min(I_Q3, I_MAX_STRUCT)

    tau_grav  = -m_exc * G * r_exc_now * np.sin(th_Q3)

    # ── MILLORA 4: guardar ω_bus per a FFT ────────────────────────────────
    omega_bus = (np.sum(I*om) + I_Q3*om_Q3) / (np.sum(I) + I_Q3)
    state['_fft_buf'].append(omega_bus)
    if len(state['_fft_buf']) > 512:
        state['_fft_buf'].pop(0)

    # ── MILLORA 5: TSR PID (cada T_TSR_PID s) ─────────────────────────────
    if int(t_now / T_TSR_PID) > int((t_now - DT) / T_TSR_PID):
        state = tsr_keeper_pid(state, v_hub)
    K_bus_eff = state['_K_bus_tsr']

    # ── MILLORA 3: K_ind adaptatiu ────────────────────────────────────────
    if int(t_now / (ADAPT_WIN)) > int((t_now-DT) / (ADAPT_WIN)):
        state, deficits = update_K_ind(state)
    else:
        deficits = compute_cluster_deficits(state)

    # ── VIA 3: ressonància ────────────────────────────────────────────────
    N_mode   = state['res_mode']
    omega_res= N_mode * DELTA_OMEGA
    tau_res  = K_RES * (omega_res - om_Q3) * I0_Q3

    # ── Finestra i transferència ──────────────────────────────────────────
    win_Q3       = window(th_Q3, N=N_ASPAS, p=p_w)
    tau_transfer = state['K_Q3_BUS'] * win_Q3 * np.sign(omega_bus - om_Q3) * I_Q3

    # ── MILLORA 2: diagnòstic càrregues ──────────────────────────────────
    ta_Q3_now = tau_aero(om_Q3, v_Q3, R=R_eff_Q3(om_Q3, v_Q3, th_Q3))
    state['_tau_buf'].append(ta_Q3_now)
    if len(state['_tau_buf']) > int(ADAPT_WIN/DT):
        state['_tau_buf'].pop(0)
    state['tau_var'] = float(np.var(state['_tau_buf'])) if len(state['_tau_buf'])>1 else 0.0

    # ── Molins ZZ: bus TSR + sub-bus ──────────────────────────────────────
    d_om = np.zeros(N_ZZ)
    for i in range(N_ZZ):
        v_eff = v_hub * (0.92 + 0.08*np.sin(i*2*np.pi/N_ZZ + om[i]*0.1))
        ta    = tau_aero(om[i], v_eff)
        tl    = tau_loss(om[i], I[i])
        tb    = K_bus_eff * (omega_bus - om[i]) * I[i]
        # Sub-bus: pes per dèficit del clúster (millora 3)
        c_i       = cids[i]
        def_i     = deficits[c_i] if c_i < len(deficits) else 0.0
        def_w     = np.clip(1.0 + 2.5*(-def_i)/max(abs(deficits).max(), DELTA_OMEGA*0.01), 0.2, 4.0)
        win_i     = window(th_Q3 - th[i], N=N_ASPAS, p=p_w)
        tau_q3_i  = state['K_Q3_IND'] * def_w * win_i * np.sign(omega_bus - om[i]) * I_Q3
        d_om[i]   = (ta - tl + tb + tau_q3_i) / I[i]

    # ── Q3 equació de moviment ─────────────────────────────────────────────
    R_q3    = R_eff_Q3(om_Q3, v_Q3, th_Q3)
    tl_Q3   = tau_loss(om_Q3, I_Q3, k_loss=0.020)
    d_om_Q3 = (ta_Q3_now + tau_grav - tl_Q3 + tau_transfer + tau_res) / I_Q3

    # ── Energia ───────────────────────────────────────────────────────────
    state['_E_acc']   += tau_transfer * om_Q3 * DT
    state['_Eg_acc']  += tau_grav     * om_Q3 * DT
    state['via1_acc'] += (ta_Q3_now - tau_aero(om_Q3, v_hub, R=R_q3)) * om_Q3 * DT
    state['via2_acc'] += tau_grav * om_Q3 * DT
    state['via3_acc'] += tau_res  * om_Q3 * DT
    if shield_on:
        state['shield_acc'] += (k_shield - 1.0) * abs(tau_transfer) * DT

    if th_Q3 - state['_last_rev'] >= 2*np.pi:
        state['E_cicle']       = state['_E_acc']
        state['E_grav_cicle']  = state['_Eg_acc']
        state['_E_acc']        = 0.0;  state['_Eg_acc'] = 0.0
        state['_last_rev']     = th_Q3

    # ── Integració ────────────────────────────────────────────────────────
    state['omega']    = np.maximum(0, om + d_om*DT)
    state['theta']    = th + state['omega']*DT
    state['omega_Q3'] = max(0, om_Q3 + d_om_Q3*DT)
    state['theta_Q3'] = th_Q3 + state['omega_Q3']*DT

    # ── Loops mid/slow/policy ──────────────────────────────────────────────
    P_B_now = power_MW(state, v_hub)
    state['_dP_buf'].append(P_B_now - P_A_now)
    state['_P_buf'].append(P_B_now)
    if len(state['_dP_buf']) > int(ADAPT_WIN/DT):
        state['_dP_buf'].pop(0)
    if len(state['_P_buf']) > int(T_POLICY/DT):
        state['_P_buf'].pop(0)

    # FFT loop (millora 4)
    if int(t_now/T_FFT) > int((t_now-DT)/T_FFT):
        state = update_res_mode_fft(state)

    # Policy loop (millora 1)
    if int(t_now/T_POLICY) > int((t_now-DT)/T_POLICY) and t_now > T_POLICY:
        wu    = max(0, len(state['_dP_buf'])//5)
        PA_sl = [P_A_now] * max(1, len(state['_P_buf']))  # aproximació
        PB_sl = state['_P_buf'][wu:]
        bal   = cluster_balance(state)
        state['_policy_buf'].append(bal)
        if len(state['_policy_buf']) > 20: state['_policy_buf'].pop(0)
        state = policy_gradient_step(
            state, t_now, PA_sl, PB_sl,
            state['_policy_buf'], state['tau_var'])

    return state


# ═══════════════════════════════════════════════════════════════════════════
# MÈTRIQUES
# ═══════════════════════════════════════════════════════════════════════════

def power_MW(state, v_hub):
    total = 0.0
    for i in range(N_ZZ):
        v_eff = v_hub*(0.92+0.08*np.sin(i*2*np.pi/N_ZZ))
        om    = state['omega'][i]
        lam   = om*R_NREL/max(v_eff,0.1)
        total += 0.5*RHO*np.pi*R_NREL**2*v_eff**3*max(cp_lambda(lam),0)
    return total/1e6

def power_Q3_MW(state, v_hub):
    v_Q3 = v_shear(v_hub, Z_Q3)
    om   = state['omega_Q3']
    R_q3 = R_eff_Q3(om, v_Q3, state['theta_Q3'])
    lam  = om*R_q3/max(v_Q3,0.1)
    return 0.5*RHO*np.pi*R_q3**2*v_Q3**3*max(cp_lambda(lam),0)/1e6

def order_parameter(state):
    return abs(np.mean(np.exp(1j*(state['theta']%(2*np.pi)))))

def omega_std(state):
    return float(np.std(state['omega']))

def sync_count(state, tol=0.12):
    om_m = np.mean(state['omega'])
    return int(np.sum(np.abs(state['omega']-om_m)<tol))

def cluster_balance(state):
    om = state['omega'];  cids = state['cluster_id']
    means = [np.mean(om[cids==c])*60/(2*np.pi) for c in range(N_CLUSTERS)]
    return max(means)-min(means)

def lambda_mean(state, v_hub):
    return np.mean(state['omega'])*R_NREL/max(v_hub,0.1)


# ═══════════════════════════════════════════════════════════════════════════
# SIMULACIÓ
# ═══════════════════════════════════════════════════════════════════════════

def run_simulation(verbose=True):
    if verbose:
        print("═"*70)
        print("  model_44_Quijote_zypyzape.py  v4")
        print("  1:PG-finit  2:càrregues  3:K_ind-adapt  4:FFT  5:TSR-PID")
        print("  Víctor Manzanares Alberola — EPSA/UPV Alcoi")
        print("═"*70)
        print(f"  T={T_SIM:.0f}s  dt={DT}s  T_policy={T_POLICY:.0f}s  "
              f"T_fft={T_FFT:.0f}s  λ_target=[{LAMBDA_MIN},{LAMBDA_MAX}]\n")

    t_arr = np.arange(T_STEPS)*DT
    v_arr = wind_profile(t_arr)

    keys  = ['P','P_tot','r','std_om','sync','E_cicle','E_grav','K_q3',
             'lambda','p_win','m_exc','res_mode','bal_clust','shield',
             'J_cost','K_ind','tau_var','K_bus_tsr']
    rec   = {s: {k: [] for k in keys} for s in ['A','B']}

    stA = init_state(7);  stB = init_state(7)
    stB['theta']     = stA['theta'].copy()
    stB['omega']     = stA['omega'].copy()
    stB['omega_nat'] = stA['omega_nat'].copy()

    log_iv = max(1, T_STEPS//20)

    for k in range(T_STEPS):
        v_hub   = v_arr[k]
        t_now   = k * DT
        P_A_now = power_MW(stA, v_hub)
        stA     = step_ZZ_only(stA, v_hub)
        stB     = step_ZZ_Q3(stB, v_hub, P_A_now, t_now)

        if k % 4 == 0:
            pA  = power_MW(stA, v_hub)
            pB  = power_MW(stB, v_hub)
            pQ3 = power_Q3_MW(stB, v_hub)

            rec['A']['P'].append(pA);         rec['B']['P'].append(pB)
            rec['A']['P_tot'].append(pA);     rec['B']['P_tot'].append(pB+pQ3)
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
            for lbl,st,kbus in [('A',stA,K_Q3_BUS_0),('B',stB,stB['K_Q3_BUS'])]:
                rec[lbl]['K_q3'].append(kbus if lbl=='A' else stB['K_Q3_BUS'])
                rec[lbl]['lambda'].append(lambda_mean(st,v_hub))
                rec[lbl]['p_win'].append(st.get('p_win',3.2))
                rec[lbl]['m_exc'].append(st.get('m_exc',M_EXC_BASE))
                rec[lbl]['res_mode'].append(st.get('res_mode',12))
                rec[lbl]['bal_clust'].append(cluster_balance(st))
                rec[lbl]['shield'].append(1.0 if st.get('shield_on',False) else 0.0)
                rec[lbl]['J_cost'].append(st.get('J_cost',0.0))
                rec[lbl]['K_ind'].append(st.get('K_Q3_IND',K_Q3_IND_0))
                rec[lbl]['tau_var'].append(st.get('tau_var',0.0))
                rec[lbl]['K_bus_tsr'].append(st.get('_K_bus_tsr',K_BUS_BASE))

        if verbose and k % log_iv == 0:
            pB = power_MW(stB,v_hub); pQ3 = power_Q3_MW(stB,v_hub)
            shld = "🛡" if stB['shield_on'] else "  "
            lam  = stB['lambda_track']
            print(f"  {100*k/T_STEPS:5.1f}%  t={t_now:6.1f}s  "
                  f"PA={power_MW(stA,v_hub):5.1f}  "
                  f"PB_tot={pB+pQ3:.1f}  "
                  f"r={order_parameter(stB):.3f}  "
                  f"λ={lam:.2f}  "
                  f"p={stB['p_win']:.2f}  "
                  f"m={stB['m_exc']/1e3:.1f}t  "
                  f"N={stB['res_mode']}  "
                  f"J={stB['J_cost']:+.2f}  "
                  f"K_tsr={stB['_K_bus_tsr']:.4f} {shld}")

    t_rec = t_arr[::4][:len(rec['A']['P'])]
    for lbl in ['A','B']:
        for k2 in rec[lbl]:
            rec[lbl][k2] = np.array(rec[lbl][k2])
    return t_rec, v_arr[::4][:len(t_rec)], rec, stA, stB


# ═══════════════════════════════════════════════════════════════════════════
# RESUM
# ═══════════════════════════════════════════════════════════════════════════

def print_summary(t_rec, rec, stB):
    print("\n"+"═"*70)
    print("  RESUM COMPARATIU  v4")
    print("═"*70)
    wu = len(t_rec)//5
    metrics = [
        ('P ZZ mitjana [MW]',       'P',       np.mean),
        ('P TOTAL ZZ+Q3 [MW]',      'P_tot',   np.mean),
        ('P màxima total [MW]',     'P_tot',   np.max),
        ('Desv. std P [MW]',        'P_tot',   np.std),
        ('Ordre r',                 'r',       np.mean),
        ('Molins sync',             'sync',    np.mean),
        ('Desv. ω [rad/s]',         'std_om',  np.mean),
        ('λ mig camp',              'lambda',  np.mean),
        ('Balanç clústers [rpm]',   'bal_clust', np.mean),
        ('Cost J política',         'J_cost',  np.mean),
    ]
    print(f"\n  {'Mètrica':<35} {'Esc. A':>10} {'Esc. B':>10} {'Δ':>10}")
    print(f"  {'-'*35} {'-'*10} {'-'*10} {'-'*10}")
    for label, key, fn in metrics:
        vA = fn(rec['A'][key][wu:]);  vB = fn(rec['B'][key][wu:])
        d  = 100*(vB-vA)/(abs(vA)+1e-12)
        print(f"  {label:<35} {vA:>10.3f} {vB:>10.3f} {'+' if d>=0 else ''}{d:>9.2f}%")

    print(f"\n  DESGLOS ENERGIA:")
    print(f"    Via 1 shear:       {stB['via1_acc']/1e6:+.1f} MJ")
    print(f"    Via 2 gravitatori: {stB['via2_acc']/1e6:+.1f} MJ")
    print(f"    Via 3 ressonància: {stB['via3_acc']/1e6:+.1f} MJ")
    print(f"    Mode escut:        {stB['shield_acc']/1e6:+.1f} MJ")
    mask = rec['B']['E_cicle']!=0
    if mask.sum()>0:
        wu2 = len(t_rec)//5
        m2  = mask & (np.arange(len(t_rec))>wu2)
        Em  = np.mean(rec['B']['E_cicle'][m2]) if m2.sum()>0 else 0
        print(f"\n  E_transfer/cicle: {Em:+.4f} MJ  {'DONANT ✓' if Em>0 else 'absorbedor'}")

    dP = np.mean(rec['B']['P_tot'][wu:]-rec['A']['P'][wu:])
    print(f"\n  ΔP net: {dP:+.3f} MW  |  Per molí: {dP/N_ZZ*1000:+.1f} kW")
    x  = stB['_policy_x']
    pw,me,kb,nm,ki = decode_policy(x)
    print(f"  Política apresa: p={pw:.2f}  m={me/1e3:.1f}t  K_bus={kb:.4f}  "
          f"N={nm}  K_ind={ki:.4f}")
    print(f"  Gradient estimat: {np.round(stB['_policy_grad'],4)}")
    print(f"  TSR PID K_bus final: {stB['_K_bus_tsr']:.5f}")
    print(f"  Var(τ_Q3): {stB['tau_var']/1e12:.2f}×10¹² N²m²  "
          f"(límit: {TAU_VAR_MAX/1e12:.2f})")
    print(f"  Mode escut: {100*np.mean(rec['B']['shield'][wu:]):.1f}% temps")
    lam_B = np.mean(rec['B']['lambda'][wu:])
    print(f"  λ_camp final: {lam_B:.3f}  target: [{LAMBDA_MIN},{LAMBDA_MAX}]  "
          f"{'✓ dins' if LAMBDA_MIN<=lam_B<=LAMBDA_MAX else '✗ fora'}")


# ═══════════════════════════════════════════════════════════════════════════
# GRÀFICS
# ═══════════════════════════════════════════════════════════════════════════

def plot_results(t_rec, v_rec, rec, stB_ref=None):
    fig = plt.figure(figsize=(22, 20))
    fig.patch.set_facecolor('#0d1117')
    gs  = gridspec.GridSpec(5, 4, figure=fig, hspace=0.52, wspace=0.40)
    BG='#161b22'; TXT='#c9d1d9'
    CA='#378ADD'; CB='#D85A30'; CV='#5DCAA5'; CQ='#FAC775'
    CG='#97C459'; CR='#ED93B1'; CM='#AFA9EC'; CO='#EF9F27'

    def mkax(r,c): return fig.add_subplot(gs[r,c])
    def style(a, title, xl='Temps [s]', yl=''):
        a.set_facecolor(BG); a.set_title(title,color=TXT,fontsize=9,pad=4)
        a.set_xlabel(xl,color=TXT,fontsize=8); a.set_ylabel(yl,color=TXT,fontsize=8)
        a.tick_params(colors=TXT,labelsize=7)
        for sp in a.spines.values(): sp.set_color('#30363d')
        a.grid(True,color='#21262d',linewidth=0.4)

    wu = len(t_rec)//5

    # 1. Potència
    a1 = mkax(0,slice(0,3))
    a1.plot(t_rec,rec['A']['P'],    color=CA,lw=0.7,alpha=0.6,label='A — ZZ sols')
    a1.plot(t_rec,rec['B']['P'],    color=CB,lw=0.7,alpha=0.6,label='B — ZZ')
    a1.plot(t_rec,rec['B']['P_tot'],color=CG,lw=1.1,          label='B_tot — ZZ+Q3')
    a1.fill_between(t_rec,rec['A']['P'],rec['B']['P_tot'],
                    where=rec['B']['P_tot']>=rec['A']['P'],alpha=0.15,color=CG)
    a1.fill_between(t_rec,rec['A']['P'],rec['B']['P_tot'],
                    where=rec['B']['P_tot']<rec['A']['P'], alpha=0.10,color=CA)
    sh = rec['B']['shield']>0.5
    if sh.sum()>0:
        a1.fill_between(t_rec,0,rec['B']['P_tot'].max()*1.05,
                        where=sh,alpha=0.06,color=CM,label='Escut')
    a1.set_xlim(t_rec[0],t_rec[-1])
    a1.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT,ncol=4)
    style(a1,f'P total v4 — PG-finit + FFT + TSR-PID + K_ind-adapt + càrregues',yl='MW')

    # 2. Vent
    a2 = mkax(0,3)
    a2.plot(t_rec,v_rec,color=CV,lw=0.8)
    a2.axhline(V_RAT,   color='white',lw=0.4,ls='--',alpha=0.4)
    a2.axhline(V_SHIELD,color=CM,    lw=0.6,ls='--',alpha=0.6,label='v_escudo')
    a2.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a2,'Vent hub',yl='m/s')

    # 3. λ camp + target
    a3 = mkax(1,0)
    a3.plot(t_rec,rec['A']['lambda'],color=CA,lw=0.8,label='A')
    a3.plot(t_rec,rec['B']['lambda'],color=CB,lw=0.9,label='B (TSR-PID)')
    a3.axhline(LAMBDA_OPT,color=CQ,lw=1,ls='--',alpha=0.8,label=f'λ_opt')
    a3.axhspan(LAMBDA_MIN,LAMBDA_MAX,alpha=0.08,color=CG,label='target')
    a3.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a3,'λ camp — TSR keeper PID (5)',yl='λ')

    # 4. K_bus TSR PID
    a4 = mkax(1,1)
    a4.plot(t_rec,rec['B']['K_bus_tsr'],color=CV,lw=0.9,label='K_bus TSR')
    a4.plot(t_rec,rec['B']['K_q3'],     color=CQ,lw=0.9,label='K_Q3_BUS')
    a4.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a4,'K_bus TSR-PID + K_Q3_BUS',yl='K')

    # 5. K_ind adaptatiu
    a5 = mkax(1,2)
    a5.plot(t_rec,rec['B']['K_ind'],color=CR,lw=0.9)
    a5.axhline(K_Q3_IND_0,color='white',lw=0.4,ls='--',alpha=0.4)
    style(a5,'K_Q3_IND adaptatiu (3)',yl='K_ind')

    # 6. Mode resonant FFT
    a6 = mkax(1,3)
    a6.plot(t_rec,rec['B']['res_mode'],color=CG,lw=0.8,drawstyle='steps-post')
    a6.set_yticks(RES_MODES);  a6.set_ylim(min(RES_MODES)-0.5,max(RES_MODES)+0.5)
    style(a6,'Mode resonant FFT (4)',yl='N')

    # 7. ΔP net
    a7 = mkax(2,slice(0,3))
    diff = rec['B']['P_tot']-rec['A']['P']
    a7.fill_between(t_rec,diff,0,where=diff>=0,color=CG,alpha=0.30,label='B_tot millor')
    a7.fill_between(t_rec,diff,0,where=diff<0, color=CA,alpha=0.18,label='A millor')
    a7.plot(t_rec,diff,color='white',lw=0.4,alpha=0.4)
    win=max(2,int(30/(DT*4)))
    if win<len(diff):
        mov=np.convolve(diff,np.ones(win)/win,mode='same')
        a7.plot(t_rec,mov,color=CQ,lw=1.8,label='Mitj. 30s')
    a7.axhline(0,color='#30363d',lw=0.8)
    a7.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a7,'ΔP = P_Btot − P_A',yl='ΔP [MW]')

    # 8. Cost J política
    a8 = mkax(2,3)
    a8.plot(t_rec,rec['B']['J_cost'],color=CO,lw=0.9)
    a8.axhline(0,color='white',lw=0.4,ls='--',alpha=0.4)
    style(a8,'Cost J política (1)',yl='J')

    # 9. p_win + m_exc
    a9 = mkax(3,0)
    ax9b = a9.twinx()
    a9.plot(t_rec,rec['B']['p_win'],color=CQ,lw=0.9,label='p_win')
    ax9b.plot(t_rec,rec['B']['m_exc']/1e3,color=CR,lw=0.9,label='m_exc')
    ax9b.set_ylabel('m [t]',color=CR,fontsize=7);  ax9b.tick_params(colors=CR,labelsize=7)
    a9.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT,loc='upper left')
    ax9b.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT,loc='lower right')
    style(a9,'p_win i m_exc (política apresa)',yl='p')

    # 10. Var(τ) càrregues
    a10 = mkax(3,1)
    tau_v_norm = rec['B']['tau_var'] / max(rec['B']['tau_var'].max(),1)
    a10.plot(t_rec,tau_v_norm,color=CR,lw=0.8)
    a10.axhline(TAU_VAR_MAX/max(rec['B']['tau_var'].max(),1),
                color='white',lw=0.6,ls='--',alpha=0.5,label='límit struct.')
    a10.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a10,'Var(τ_Q3) càrregues (2)',yl='norm.')

    # 11. Ordre r + sync
    a11 = mkax(3,2)
    a11.plot(t_rec,rec['A']['r'],color=CA,lw=0.8,label='r A')
    a11.plot(t_rec,rec['B']['r'],color=CB,lw=0.9,label='r B')
    a11.set_ylim(0,1.05)
    ax11b = a11.twinx()
    ax11b.plot(t_rec,rec['B']['sync'],color=CG,lw=0.7,alpha=0.6,label='sync B')
    ax11b.set_ylim(0,N_ZZ+2);  ax11b.tick_params(colors=CG,labelsize=7)
    a11.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT,loc='upper left')
    style(a11,'Ordre r + Sync',yl='r')

    # 12. Desv ω
    a12 = mkax(3,3)
    a12.plot(t_rec,rec['A']['std_om'],color=CA,lw=0.8,label='A')
    a12.plot(t_rec,rec['B']['std_om'],color=CB,lw=0.8,label='B')
    a12.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a12,'Desv. ω',yl='rad/s')

    # 13. Cp(λ) punts treball
    a13 = mkax(4,slice(0,2))
    lam_arr=np.linspace(1,18,400)
    a13.plot(lam_arr,[cp_lambda(l) for l in lam_arr],color=CV,lw=1.5)
    a13.axvline(LAMBDA_OPT,color=CQ,lw=0.8,ls='--',alpha=0.7,label='λ_opt')
    a13.axvspan(LAMBDA_MIN,LAMBDA_MAX,alpha=0.10,color=CG,label='target TSR')
    wu2=len(t_rec)//5
    lA=np.mean(rec['A']['lambda'][wu2:]); lB=np.mean(rec['B']['lambda'][wu2:])
    a13.axvline(lA,color=CA,lw=1.2,alpha=0.8,label=f'λ_A={lA:.2f}')
    a13.axvline(lB,color=CB,lw=1.2,alpha=0.8,label=f'λ_B={lB:.2f}')
    a13.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a13,'Cp(λ) — target TSR PID',xl='λ',yl='Cp')

    # 14. Taula resum millores
    a14 = mkax(4,slice(2,4))
    a14.set_facecolor(BG);  a14.axis('off')
    wu3=len(t_rec)//5
    dP=np.mean(rec['B']['P_tot'][wu3:]-rec['A']['P'][wu3:])
    lam_f=np.mean(rec['B']['lambda'][wu3:])
    via1_tmp = stB_ref["via1_acc"]/1e6 if stB_ref else 0
    via2_tmp = stB_ref["via2_acc"]/1e6 if stB_ref else 0
    items = [
        ("MILLORES v4 ACTIVES", CQ, True),
        (f"1. Política PG-finit: η={ETA_POLICY}  ε={EPS_GRAD}", TXT, False),
        (f"2. Càrregues: τ_var_max={TAU_VAR_MAX/1e12:.1f}×10¹²", TXT, False),
        (f"3. K_ind adapt: γ=2.5 cluster_balance", TXT, False),
        (f"4. FFT ressonància cada {T_FFT:.0f}s", TXT, False),
        (f"5. TSR PID Kp={TSR_KP} Ki={TSR_KI} Kd={TSR_KD}", TXT, False),
        ("", TXT, False),
        ("RESULTAT FINAL", CQ, True),
        (f"ΔP net: {dP:+.3f} MW", CG if dP>0 else CR, False),
        (f"Per molí: {dP/N_ZZ*1000:+.0f} kW", CG if dP>0 else CR, False),
        (f"λ_camp: {lam_f:.3f}  target [{LAMBDA_MIN},{LAMBDA_MAX}]",
         CG if LAMBDA_MIN<=lam_f<=LAMBDA_MAX else CO, False),
        (f"Shear: +{via1_tmp:.0f} MJ | Grav: {via2_tmp:.0f} MJ",
         TXT, False),
    ]
    # obtenim via1/via2 del stB (globals de la funció)
    via1 = stB_ref["via1_acc"]/1e6 if stB_ref else 0
    via2 = stB_ref["via2_acc"]/1e6 if stB_ref else 0
    items[-1] = (f"Shear: +{via1_tmp:.0f} MJ | Grav: {via2_tmp:.0f} MJ", TXT, False)
    for li,(ln,col,bold) in enumerate(items):
        a14.text(0.04, 0.97-li*0.077, ln, color=col, fontsize=8.5,
                 fontfamily='monospace',
                 fontweight='bold' if bold else 'normal',
                 transform=a14.transAxes)

    fig.suptitle(
        f'Gemell Digital v4 — {N_ZZ}×ZZ(3A)+Q3  |  '
        'PG-finit · FFT-mode · TSR-PID · K_ind-adapt · load-constraint  |  '
        f'T={T_SIM:.0f}s',
        color=TXT, fontsize=9, y=0.999)

    out='/mnt/user-data/outputs/model_44_Quijote_zypyzape_v4.png'
    plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=fig.get_facecolor())
    print(f"  Gràfic: {out}")
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    t_rec, v_rec, rec, stA, stB = run_simulation(verbose=True)
    print_summary(t_rec, rec, stB)
    print("\n  Generant gràfics v4...")
    plot_results(t_rec, v_rec, rec, stB)
    import shutil
    shutil.copy(__file__, '/mnt/user-data/outputs/model_44_Quijote_zypyzape.py')
    print("  ✓ v4 completada.")
