"""
model_44_Quijote_zypyzape.py  v2.1
═══════════════════════════════════════════════════════════════════════════════
Gemell digital: 44 molins ZypyZape (3 aspas, NREL 5MW) + Quijote (3 aspas)
Bus comú d'inèrcia sintètica — Lògica gemell_3vs7_rigoros adaptada a 3A

MILLORES v2.1 (basades en revisió externa + diagnòstic de signe):
  1. K_Q3_BUS=0.05 ↓ / K_Q3_IND=0.09 ↑  → hurto local molí-a-molí
  2. α=0.22 (moderat)                     → menys cost energètic inercia
  3. Hurto gravitatori DIRECCIONAL corregit:
       τ_transfer = K · |sin(3θ_Q3)| · sign(ω_bus − ω_Q3) · I_Q3(θ)
     La finestra |sin(3θ)| = "quan les aspas estan en posició favorable"
     sign(ω_bus−ω_Q3) = "sempre empeny cap al bus, no contra"
     Això és el hurto gravitatori real: transferència direccional + modulació
  4. Controlador adaptatiu sobre K_Q3_BUS (ΔP finestra 30s)
  5. Escaneig paramètric: K_bus × K_ind × α  (mapa de fase ΔP)

Comparació:
  [A] Camp sol : 44 ZypyZape sense Quijote
  [B] Camp+Q3  : 44 ZypyZape + Quijote adaptatiu

Autor: Víctor Manzanares Alberola — EPSA / UPV (Alcoi)
Repositori: github.com/espiradesombra/claude
═══════════════════════════════════════════════════════════════════════════════
FÍSICA DEL HURTO GRAVITATORI (versió corregida)
────────────────────────────────────────────────
Quijote actua com a "volandum orbital": absorbeix energia de la corrent quan
la seva fase és desfavorable (aspas en posició de drag) i la retorna al bus
quan és favorable (aspas en posició de lift màxim).

  I_Q3(θ) = I₀ · (1 + α·sin(3θ))       ← inercia sintètica Fe+oli
  τ_transfer = K · |sin(3θ)| · sign(ω_bus − ω_Q3) · I_Q3(θ)

La diferència clau respecte a Kuramoto simple:
  - Kuramoto: sin(θ_j − θ_i)           → simètric, pot cancel·lar-se
  - Hurto:    |sin(3θ)| · sign(Δω)     → sempre direccional, mai negatiu net

Injecció als molins individuals (3A acoblades):
  τ_Q3→i = K_ind · |sin(3θ_Q3 − 3θ_i)| · sign(ω_bus − ω_i) · I_Q3
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

# ─── PARÀMETRES NREL 5MW ─────────────────────────────────────────────────────
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

N_ZZ       = 44
N_ASPAS    = 3
K_BUS      = 0.08

# ─── QUIJOTE v2.1 ─────────────────────────────────────────────────────────────
I0_Q3      = J_NREL * 1.25
ALPHA_Q3   = 0.22
K_Q3_BUS   = 0.05
K_Q3_IND   = 0.09

# ─── CONTROLADOR ─────────────────────────────────────────────────────────────
ADAPT_WIN  = 30.0
ADAPT_KMIN = 0.02
ADAPT_KMAX = 0.15

DT         = 0.05
T_SIM      = 600.0
T_STEPS    = int(T_SIM / DT)

WIND_BASE      = 9.0
WIND_TURB      = 0.12
WIND_GUST_T    = 200.0
WIND_GUST_V    = 14.5
WIND_GUST_DUR  = 60.0


# ═══════════════════════════════════════════════════════════════════════════════
def cp_lambda(lam):
    return CP_MAX * np.exp(-((lam - LAMBDA_OPT) / SIGMA_CP) ** 2)

def tau_aero(omega, v_eff, R=R_NREL, rho=RHO):
    omega = max(omega, 0.01);  v_eff = max(v_eff, 0.1)
    lam = omega * R / v_eff
    P   = 0.5 * rho * np.pi * R**2 * v_eff**3 * max(cp_lambda(lam), 0.0)
    return P / omega

def tau_loss(omega, I, k_loss=0.03):
    return k_loss * I * omega

def wind_profile(t_arr):
    rng  = np.random.default_rng(42)
    turb = rng.normal(0, WIND_TURB * WIND_BASE, len(t_arr))
    kern = np.ones(max(1, int(5/DT))) / max(1, int(5/DT))
    turb = np.convolve(turb, kern, mode='same')
    gust = WIND_GUST_V * np.exp(-0.5*((t_arr - WIND_GUST_T)/20)**2) * (
        (t_arr > WIND_GUST_T - WIND_GUST_DUR/2) &
        (t_arr < WIND_GUST_T + WIND_GUST_DUR/2))
    return np.clip(WIND_BASE + turb + gust, V_CUT_IN, V_CUT_OUT)

def init_state(seed=7):
    rng = np.random.default_rng(seed)
    omega_nat = OMEGA_RAT * (1.0 + rng.uniform(-0.08, 0.08, N_ZZ))
    return {
        'theta':     rng.uniform(0, 2*np.pi, N_ZZ),
        'omega':     omega_nat * rng.uniform(0.85, 1.05, N_ZZ),
        'omega_nat': omega_nat,
        'I':         np.full(N_ZZ, J_NREL),
        'theta_Q3':  0.0,
        'omega_Q3':  OMEGA_RAT * 0.9,
        'E_cicle':   0.0,
        '_E_acc':    0.0,
        '_last_rev': 0.0,
        'K_Q3_BUS':  K_Q3_BUS,
        '_dP_buf':   [],
    }


# ═══════════════════════════════════════════════════════════════════════════════
def step_ZZ_only(state, v_wind):
    th = state['theta'];  om = state['omega'];  I = state['I']
    omega_bus = np.sum(I * om) / np.sum(I)
    d_om = np.zeros(N_ZZ)
    for i in range(N_ZZ):
        v_eff = v_wind * (0.92 + 0.08*np.sin(i*2*np.pi/N_ZZ + om[i]*0.1))
        d_om[i] = (tau_aero(om[i], v_eff)
                   - tau_loss(om[i], I[i])
                   + K_BUS * (omega_bus - om[i]) * I[i]) / I[i]
    state['omega'] = np.maximum(0, om + d_om * DT)
    state['theta'] = th + state['omega'] * DT
    return state


def step_ZZ_Q3(state, v_wind, P_A_now,
               kbus_ovr=None, kind_ovr=None, alpha_ovr=None):
    th    = state['theta'];  om = state['omega'];  I = state['I']
    th_Q3 = state['theta_Q3'];  om_Q3 = state['omega_Q3']
    Kbus  = kbus_ovr  if kbus_ovr  is not None else state['K_Q3_BUS']
    Kind  = kind_ovr  if kind_ovr  is not None else K_Q3_IND
    alpha = alpha_ovr if alpha_ovr is not None else ALPHA_Q3

    # Inercia sintètica Q3
    I_Q3 = I0_Q3 * (1.0 + alpha * np.sin(N_ASPAS * th_Q3))

    # Velocitat bus (inclou Q3)
    omega_bus = (np.sum(I * om) + I_Q3 * om_Q3) / (np.sum(I) + I_Q3)

    # ── Molins ZZ: bus + injecció directa Q3 ─────────────────────────────────
    d_om = np.zeros(N_ZZ)
    for i in range(N_ZZ):
        v_eff = v_wind * (0.92 + 0.08*np.sin(i*2*np.pi/N_ZZ + om[i]*0.1))
        ta = tau_aero(om[i], v_eff)
        tl = tau_loss(om[i], I[i])
        tb = K_BUS * (omega_bus - om[i]) * I[i]
        # Hurto local Q3→molí: direccional + finestra 3A
        win_i    = abs(np.sin(N_ASPAS * th_Q3 - N_ASPAS * th[i]))
        tau_q3_i = Kind * win_i * np.sign(omega_bus - om[i]) * I_Q3
        d_om[i]  = (ta - tl + tb + tau_q3_i) / I[i]

    # ── Quijote: hurto gravitatori direccional ────────────────────────────────
    # |sin(3θ)| = finestra de fase favorable (3 aspas)
    # sign(ω_bus − ω_Q3) = direcció correcta sempre
    win_Q3       = abs(np.sin(N_ASPAS * th_Q3))
    tau_transfer = Kbus * win_Q3 * np.sign(omega_bus - om_Q3) * I_Q3

    ta_Q3   = tau_aero(om_Q3, v_wind * 1.05, R=R_NREL*1.1)
    tl_Q3   = tau_loss(om_Q3, I_Q3, k_loss=0.025)
    d_om_Q3 = (ta_Q3 - tl_Q3 + tau_transfer) / I_Q3

    # Energia per cicle
    state['_E_acc'] += tau_transfer * om_Q3 * DT
    if th_Q3 - state['_last_rev'] >= 2*np.pi:
        state['E_cicle']   = state['_E_acc']
        state['_E_acc']    = 0.0
        state['_last_rev'] = th_Q3

    # Integració
    state['omega']    = np.maximum(0, om + d_om * DT)
    state['theta']    = th + state['omega'] * DT
    state['omega_Q3'] = max(0, om_Q3 + d_om_Q3 * DT)
    state['theta_Q3'] = th_Q3 + state['omega_Q3'] * DT

    # ── Controlador adaptatiu ─────────────────────────────────────────────────
    if kbus_ovr is None:
        P_B_now = power_MW(state, v_wind)
        state['_dP_buf'].append(P_B_now - P_A_now)
        win_s = max(1, int(ADAPT_WIN / DT))
        if len(state['_dP_buf']) > win_s:
            state['_dP_buf'].pop(0)
        avg_dP = np.mean(state['_dP_buf'])
        r_now  = order_parameter(state)
        K_cur  = state['K_Q3_BUS']
        if avg_dP < 0:
            K_cur *= 0.98
        elif avg_dP > 0 and r_now > 0.15:
            K_cur *= 1.01
        state['K_Q3_BUS'] = float(np.clip(K_cur, ADAPT_KMIN, ADAPT_KMAX))

    return state


# ═══════════════════════════════════════════════════════════════════════════════
def power_MW(state, v_wind):
    total = 0.0
    for i in range(N_ZZ):
        v_eff = v_wind * (0.92 + 0.08*np.sin(i*2*np.pi/N_ZZ))
        om    = state['omega'][i]
        lam   = om * R_NREL / max(v_eff, 0.1)
        total += 0.5*RHO*np.pi*R_NREL**2*v_eff**3*max(cp_lambda(lam), 0)
    return total / 1e6

def order_parameter(state):
    return abs(np.mean(np.exp(1j * (state['theta'] % (2*np.pi)))))

def omega_std(state):
    return float(np.std(state['omega']))

def sync_count(state, tol=0.15):
    om_mean = np.mean(state['omega'])
    return int(np.sum(np.abs(state['omega'] - om_mean) < tol))


# ═══════════════════════════════════════════════════════════════════════════════
def run_simulation(verbose=True):
    if verbose:
        print("═"*70)
        print("  model_44_Quijote_zypyzape.py  v2.1")
        print("  Hurto gravitatori DIRECCIONAL: |sin(3θ)|·sign(Δω)")
        print("  Víctor Manzanares Alberola — EPSA/UPV Alcoi")
        print("═"*70)
        print(f"\n  α={ALPHA_Q3}  K_bus={K_Q3_BUS}  K_ind={K_Q3_IND}  "
              f"T={T_SIM:.0f}s  dt={DT}s  N={N_ZZ}\n")

    t_arr = np.arange(T_STEPS) * DT
    v_arr = wind_profile(t_arr)
    rec   = {s: {k: [] for k in ['P','r','std_om','sync','E_cicle','K_q3']}
             for s in ['A','B']}

    stA = init_state(7)
    stB = init_state(7)
    stB['theta'] = stA['theta'].copy()
    stB['omega'] = stA['omega'].copy()
    stB['omega_nat'] = stA['omega_nat'].copy()

    log_iv = max(1, T_STEPS // 20)
    for k in range(T_STEPS):
        v       = v_arr[k]
        P_A_now = power_MW(stA, v)
        stA     = step_ZZ_only(stA, v)
        stB     = step_ZZ_Q3(stB, v, P_A_now)

        if k % 4 == 0:
            for lbl, st in [('A', stA), ('B', stB)]:
                rec[lbl]['P'].append(power_MW(st, v))
                rec[lbl]['r'].append(order_parameter(st))
                rec[lbl]['std_om'].append(omega_std(st))
                rec[lbl]['sync'].append(sync_count(st))
                rec[lbl]['E_cicle'].append(
                    st['E_cicle']/1e6 if lbl=='B' else 0.0)
                rec[lbl]['K_q3'].append(st.get('K_Q3_BUS', K_Q3_BUS))

        if verbose and k % log_iv == 0:
            print(f"  {100*k/T_STEPS:5.1f}%  t={k*DT:6.1f}s  "
                  f"PA={power_MW(stA,v):6.1f}  PB={power_MW(stB,v):6.1f}  "
                  f"rA={order_parameter(stA):.3f}  rB={order_parameter(stB):.3f}  "
                  f"K_q3={stB['K_Q3_BUS']:.4f}")

    t_rec = t_arr[::4][:len(rec['A']['P'])]
    for lbl in ['A','B']:
        for k2 in rec[lbl]:
            rec[lbl][k2] = np.array(rec[lbl][k2])
    return t_rec, v_arr[::4][:len(t_rec)], rec, stA, stB


# ═══════════════════════════════════════════════════════════════════════════════
def param_scan():
    K_bus_vals = [0.02, 0.04, 0.06, 0.10]
    K_ind_vals = [0.04, 0.09, 0.14]
    alpha_vals = [0.15, 0.22, 0.35]
    T_S        = 300.0
    STEPS_S    = int(T_S / DT)
    t_s        = np.arange(STEPS_S) * DT
    v_s        = wind_profile(t_s)
    results    = []
    total      = len(K_bus_vals)*len(K_ind_vals)*len(alpha_vals)
    done       = 0
    print(f"\n  Escaneig paramètric: {total} combinacions × {T_S:.0f}s...")

    for kb, ki, al in product(K_bus_vals, K_ind_vals, alpha_vals):
        stA = init_state(7)
        stB = init_state(7)
        stB['theta'] = stA['theta'].copy()
        stB['omega'] = stA['omega'].copy()
        stB['omega_nat'] = stA['omega_nat'].copy()
        PA_l, PB_l, rB_l = [], [], []

        for k in range(STEPS_S):
            v       = v_s[k]
            PA_now  = power_MW(stA, v)
            stA     = step_ZZ_only(stA, v)
            stB     = step_ZZ_Q3(stB, v, PA_now,
                                 kbus_ovr=kb, kind_ovr=ki, alpha_ovr=al)
            if k % 4 == 0:
                PA_l.append(power_MW(stA, v))
                PB_l.append(power_MW(stB, v))
                rB_l.append(order_parameter(stB))

        wu = len(PA_l) // 5
        dP = np.mean(np.array(PB_l[wu:]) - np.array(PA_l[wu:]))
        rB = np.mean(rB_l[wu:])
        results.append((kb, ki, al, dP, rB))
        done += 1
        if done % 6 == 0:
            print(f"    {done}/{total}  kb={kb:.2f} ki={ki:.2f} α={al:.2f}  "
                  f"ΔP={dP:+.2f} MW  r_B={rB:.3f}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
def print_summary(t_rec, rec):
    print("\n" + "═"*70)
    print("  RESUM COMPARATIU  v2.1")
    print("═"*70)
    wu = len(t_rec) // 5
    metrics = [
        ('Potència mitjana [MW]',       'P',      np.mean),
        ('Potència màxima [MW]',        'P',      np.max),
        ('Desv. estàndard P [MW]',      'P',      np.std),
        ('Ordre sincronització r',      'r',      np.mean),
        ('Molins sincronitzats (mitj)', 'sync',   np.mean),
        ('Desv. ω [rad/s]',             'std_om', np.mean),
    ]
    print(f"\n  {'Mètrica':<38} {'Esc. A':>10} {'Esc. B':>10} {'Millora':>10}")
    print(f"  {'-'*38} {'-'*10} {'-'*10} {'-'*10}")
    for label, key, fn in metrics:
        vA = fn(rec['A'][key][wu:])
        vB = fn(rec['B'][key][wu:])
        d  = 100*(vB-vA)/(abs(vA)+1e-12)
        print(f"  {label:<38} {vA:>10.3f} {vB:>10.3f} {'+' if d>=0 else ''}{d:>9.2f}%")

    mask = rec['B']['E_cicle'][wu:] != 0
    if mask.sum() > 0:
        E_m = np.mean(rec['B']['E_cicle'][wu:][mask])
        print(f"\n  E_cicle Q3: {E_m:+.4f} MJ/rev  "
              f"{'→ DONANT ✓' if E_m>0 else '→ absorbedor'}")

    dP_m = np.mean(rec['B']['P'][wu:] - rec['A']['P'][wu:])
    print(f"  Guany net: {dP_m:+.3f} MW  |  Per molí: {dP_m/N_ZZ*1000:+.1f} kW")
    kb = rec['B']['K_q3']
    print(f"  K_Q3_BUS: inici={kb[0]:.4f} → final={kb[-1]:.4f}  "
          f"rang=[{kb.min():.4f}, {kb.max():.4f}]")


# ═══════════════════════════════════════════════════════════════════════════════
def plot_results(t_rec, v_rec, rec, scan_results):
    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor('#0e1117')
    gs  = gridspec.GridSpec(4, 4, figure=fig, hspace=0.52, wspace=0.42)
    BG  = '#161b22';  TXT = '#c9d1d9'
    CA  = '#378ADD';  CB  = '#D85A30'
    CV  = '#5DCAA5';  CQ  = '#FAC775'

    def mkax(r, c): return fig.add_subplot(gs[r, c])
    def style(a, title, xl='Temps [s]', yl=''):
        a.set_facecolor(BG);  a.set_title(title, color=TXT, fontsize=9, pad=5)
        a.set_xlabel(xl, color=TXT, fontsize=8)
        a.set_ylabel(yl, color=TXT, fontsize=8)
        a.tick_params(colors=TXT, labelsize=7)
        for sp in a.spines.values(): sp.set_color('#30363d')
        a.grid(True, color='#21262d', linewidth=0.4)

    wu = len(t_rec) // 5

    # 1. Potència
    a1 = mkax(0, slice(0,3))
    a1.plot(t_rec, rec['A']['P'], color=CA, lw=0.9, alpha=0.7, label='A — sense Quijote')
    a1.plot(t_rec, rec['B']['P'], color=CB, lw=0.9, alpha=0.9, label='B — Quijote adaptatiu v2.1')
    a1.fill_between(t_rec, rec['A']['P'], rec['B']['P'],
                    where=rec['B']['P']>=rec['A']['P'], alpha=0.18, color=CB, label='Guany B>A')
    a1.fill_between(t_rec, rec['A']['P'], rec['B']['P'],
                    where=rec['B']['P']<rec['A']['P'], alpha=0.10, color=CA, label='Cost B<A')
    a1.axvspan(t_rec[0], t_rec[wu], color='white', alpha=0.02)
    a1.set_xlim(t_rec[0], t_rec[-1])
    a1.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT, ncol=4)
    style(a1, f'Potència total — {N_ZZ}×ZZ(3A) + Quijote direccional', yl='MW')

    # 2. Vent
    a2 = mkax(0,3)
    a2.plot(t_rec, v_rec, color=CV, lw=0.8)
    a2.axhline(V_RAT, color='white', lw=0.4, ls='--', alpha=0.4)
    style(a2, 'Vent', yl='m/s')

    # 3. K adaptatiu
    a3 = mkax(1,0)
    a3.plot(t_rec, rec['B']['K_q3'], color=CQ, lw=0.9)
    a3.axhline(K_Q3_BUS, color='white', lw=0.5, ls='--', alpha=0.4, label=f'K0={K_Q3_BUS}')
    a3.set_ylim(ADAPT_KMIN*0.9, ADAPT_KMAX*1.05)
    a3.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a3, 'K_Q3_BUS adaptatiu', yl='K')

    # 4. Ordre r
    a4 = mkax(1,1)
    a4.plot(t_rec, rec['A']['r'], color=CA, lw=0.9, label='A')
    a4.plot(t_rec, rec['B']['r'], color=CB, lw=0.9, label='B')
    a4.set_ylim(0, 1.05)
    a4.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a4, 'Ordre Kuramoto r', yl='r')

    # 5. Sync
    a5 = mkax(1,2)
    a5.plot(t_rec, rec['A']['sync'], color=CA, lw=0.9, drawstyle='steps-post', label='A')
    a5.plot(t_rec, rec['B']['sync'], color=CB, lw=0.9, drawstyle='steps-post', label='B')
    a5.axhline(N_ZZ, color='white', lw=0.4, ls=':', alpha=0.3)
    a5.set_ylim(0, N_ZZ+2)
    a5.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a5, f'Molins sincronitzats /{N_ZZ}', yl='N')

    # 6. Hurto
    a6 = mkax(1,3)
    mask = rec['B']['E_cicle'] != 0
    if mask.sum() > 0:
        colors_e = np.where(rec['B']['E_cicle'][mask]>=0, CB, CA)
        a6.scatter(t_rec[mask], rec['B']['E_cicle'][mask], c=colors_e, s=5, alpha=0.7)
        a6.axhline(0, color='white', lw=0.5, ls='--', alpha=0.4)
        em = np.mean(rec['B']['E_cicle'][mask])
        a6.axhline(em, color=CQ, lw=1.2, ls='--',
                   label=f'mitj={em:+.3f} MJ')
        a6.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a6, 'Hurto gravitatori (E/cicle)', yl='MJ/rev')

    # 7. ΔP
    a7 = mkax(2, slice(0,3))
    diff = rec['B']['P'] - rec['A']['P']
    a7.fill_between(t_rec, diff, 0, where=diff>=0, color=CB, alpha=0.30, label='B millor')
    a7.fill_between(t_rec, diff, 0, where=diff<0,  color=CA, alpha=0.20, label='A millor')
    a7.plot(t_rec, diff, color='white', lw=0.5, alpha=0.4)
    win = max(2, int(30/(DT*4)))
    if win < len(diff):
        mov = np.convolve(diff, np.ones(win)/win, mode='same')
        a7.plot(t_rec, mov, color=CQ, lw=1.8, alpha=0.85, label='Mitj. 30s')
    a7.axhline(0, color='#30363d', lw=0.8)
    a7.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a7, 'ΔP = P_B − P_A  (guany Quijote)', yl='ΔP [MW]')

    # 8. Desv ω
    a8 = mkax(2,3)
    a8.plot(t_rec, rec['A']['std_om'], color=CA, lw=0.9, label='A')
    a8.plot(t_rec, rec['B']['std_om'], color=CB, lw=0.9, label='B')
    a8.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a8, 'Dispersió ω', yl='rad/s')

    # 9–11. Mapes de fase escaneig
    alpha_vals_u = sorted(set(r[2] for r in scan_results))
    K_bus_u = sorted(set(r[0] for r in scan_results))
    K_ind_u = sorted(set(r[1] for r in scan_results))

    for ai, al in enumerate(alpha_vals_u[:3]):
        a_sc = mkax(3, ai)
        mat  = np.zeros((len(K_ind_u), len(K_bus_u)))
        for rs in scan_results:
            if abs(rs[2]-al) < 1e-6:
                bi = K_bus_u.index(rs[0])
                ii = K_ind_u.index(rs[1])
                mat[ii, bi] = rs[3]
        vmax = max(abs(mat).max(), 0.5)
        im   = a_sc.imshow(mat, aspect='auto', cmap='RdYlGn',
                           vmin=-vmax, vmax=vmax, origin='lower')
        a_sc.set_xticks(range(len(K_bus_u)))
        a_sc.set_xticklabels([f'{v:.2f}' for v in K_bus_u], color=TXT, fontsize=7)
        a_sc.set_yticks(range(len(K_ind_u)))
        a_sc.set_yticklabels([f'{v:.2f}' for v in K_ind_u], color=TXT, fontsize=7)
        for ii2 in range(len(K_ind_u)):
            for bi2 in range(len(K_bus_u)):
                v = mat[ii2, bi2]
                a_sc.text(bi2, ii2, f'{v:+.1f}', ha='center', va='center',
                          color='white', fontsize=7, fontweight='bold')
        a_sc.set_facecolor(BG)
        a_sc.set_title(f'Mapa ΔP [MW]  α={al:.2f}', color=TXT, fontsize=9, pad=5)
        a_sc.set_xlabel('K_Q3_BUS', color=TXT, fontsize=8)
        a_sc.set_ylabel('K_Q3_IND', color=TXT, fontsize=8)
        a_sc.tick_params(colors=TXT, labelsize=7)
        for sp in a_sc.spines.values(): sp.set_color('#30363d')
        plt.colorbar(im, ax=a_sc, fraction=0.046, pad=0.04).ax.tick_params(
            colors=TXT, labelsize=7)

    # 12. Millors paràmetres (text)
    a12 = mkax(3, 3)
    a12.set_facecolor(BG);  a12.axis('off')
    best  = max(scan_results, key=lambda x: x[3])
    worst = min(scan_results, key=lambda x: x[3])
    lines = [
        "MILLORS PARÀMETRES",
        f"K_bus = {best[0]:.2f}",
        f"K_ind = {best[1]:.2f}",
        f"α     = {best[2]:.2f}",
        f"ΔP    = {best[3]:+.2f} MW",
        f"r_B   = {best[4]:.3f}",
        "",
        "PITJORS",
        f"K_bus = {worst[0]:.2f}",
        f"K_ind = {worst[1]:.2f}",
        f"α     = {worst[2]:.2f}",
        f"ΔP    = {worst[3]:+.2f} MW",
    ]
    for li, ln in enumerate(lines):
        color = CQ if li == 0 else (CA if li >= 8 else TXT)
        a12.text(0.05, 0.95 - li*0.073, ln, transform=a12.transAxes,
                 color=color, fontsize=9,
                 fontfamily='monospace')

    fig.suptitle(
        f'Gemell Digital v2.1 — {N_ZZ}×ZZ(3A) + Quijote  |  '
        f'Hurto: |sin(3θ)|·sign(Δω)  |  T={T_SIM:.0f}s',
        color=TXT, fontsize=10, y=0.998)

    out = '/mnt/user-data/outputs/model_44_Quijote_zypyzape_v2.png'
    plt.savefig(out, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    print(f"  Gràfic guardat: {out}")
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    t_rec, v_rec, rec, stA, stB = run_simulation(verbose=True)
    print_summary(t_rec, rec)
    scan = param_scan()
    best = max(scan, key=lambda x: x[3])
    print(f"\n  MILLORS: K_bus={best[0]:.2f}  K_ind={best[1]:.2f}  "
          f"α={best[2]:.2f}  →  ΔP={best[3]:+.3f} MW  r={best[4]:.3f}")
    print("\n  Generant gràfics...")
    plot_results(t_rec, v_rec, rec, scan)
    print("\n  ✓ v2.1 completada.")
