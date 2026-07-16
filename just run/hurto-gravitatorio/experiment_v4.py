"""
experiment_v4.py
═══════════════════════════════════════════════════════════════════════════════
Script d'experiment controlat per model_44_Quijote_zypyzape.py v4

Implements:
  1. Baseline multi-seed (10 seeds) amb paràmetres v4
  2. Estabilitat política: ε=0.04, η=0.03, monitorització J
  3. Monte Carlo robustesa vent (50 runs, turbulència variable)
  4. Ablation study: desactiva FFT / TSR-PID / K_ind / load-pen un per un
  5. Surrogate offline: recull (x, J) → Random Forest → cerca Bayesiana simple

Sortides:
  experiment_results.json  — totes les dades
  experiment_summary.png   — gràfics comparatius
  best_policy.txt          — millors paràmetres trobats

Autor: Víctor Manzanares Alberola — EPSA/UPV Alcoi
═══════════════════════════════════════════════════════════════════════════════
"""

import numpy as np
import json, time, copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ── importa el model ──────────────────────────────────────────────────────────
import importlib.util, sys, os
spec = importlib.util.spec_from_file_location(
    "model", os.path.join(os.path.dirname(__file__),
                          "model_44_Quijote_zypyzape.py"))
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# ── paràmetres experiment ─────────────────────────────────────────────────────
N_SEEDS_BASE  = 10
N_SEEDS_MC    = 50
T_SHORT       = 300.0   # s per experim. ràpids
T_FULL        = 600.0
DT            = mod.DT

# Valors reducció soroll (revisió §5.1)
EPS_REDUCED   = 0.04
ETA_REDUCED   = 0.03
N_REPEAT_PERT = 3       # repeticions per pertorbació

BG='#161b22'; TXT='#c9d1d9'
CA='#378ADD'; CB='#D85A30'; CV='#5DCAA5'; CQ='#FAC775'
CG='#97C459'; CR='#ED93B1'; CM='#AFA9EC'; CO='#EF9F27'

results = {
    'baseline':  [],
    'policy':    [],
    'montecarlo':[],
    'ablation':  {},
    'surrogate': {},
}


# ═══════════════════════════════════════════════════════════════════════════════
# UTILS
# ═══════════════════════════════════════════════════════════════════════════════

def run_short(seed=7, T=T_SHORT,
              wind_turb=None, wind_gust_v=None,
              disable_fft=False, disable_tsr=False,
              disable_kind=False, disable_load=False,
              eps=None, eta=None,
              verbose=False):
    """
    Executa escenaris A+B per T segons amb seed donat.
    Retorna dict de mètriques post-warmup.
    """
    # Backup i override paràmetres del model
    orig = {}
    overrides = {}
    if wind_turb  is not None: overrides['WIND_TURB']    = wind_turb
    if wind_gust_v is not None: overrides['WIND_GUST_V'] = wind_gust_v
    if eps is not None:         overrides['EPS_GRAD']    = eps
    if eta is not None:         overrides['ETA_POLICY']  = eta

    for k, v in overrides.items():
        orig[k] = getattr(mod, k)
        setattr(mod, k, v)

    STEPS = int(T / DT)
    t_arr = np.arange(STEPS) * DT
    v_arr = mod.wind_profile(t_arr)

    stA = mod.init_state(seed)
    stB = mod.init_state(seed)
    stB['theta']     = stA['theta'].copy()
    stB['omega']     = stA['omega'].copy()
    stB['omega_nat'] = stA['omega_nat'].copy()

    PA_l, PB_l, PQ3_l = [], [], []
    r_l, std_l, sync_l = [], [], []
    lam_l, tau_l, J_l  = [], [], []
    K_bus_l, K_ind_l   = [], []
    E_cyc_l            = []

    for k in range(STEPS):
        v_hub   = v_arr[k]
        t_now   = k * DT
        P_A_now = mod.power_MW(stA, v_hub)

        stA = mod.step_ZZ_only(stA, v_hub)

        # Ablació: sobreescriu funcions si disable_*
        if disable_fft:
            stB_bk = stB.get('res_mode', 12)  # manté mode fix
        if disable_tsr:
            stB['_K_bus_tsr'] = mod.K_BUS_BASE

        stB = mod.step_ZZ_Q3(stB, v_hub, P_A_now, t_now)

        if disable_fft:
            stB['res_mode'] = stB_bk
        if disable_tsr:
            stB['_K_bus_tsr'] = mod.K_BUS_BASE
        if disable_kind:
            stB['K_Q3_IND'] = mod.K_Q3_IND_0
        if disable_load:
            stB['tau_var'] = 0.0

        if k % 4 == 0:
            pA  = mod.power_MW(stA, v_hub)
            pB  = mod.power_MW(stB, v_hub)
            pQ3 = mod.power_Q3_MW(stB, v_hub)
            PA_l.append(pA);   PB_l.append(pB);   PQ3_l.append(pQ3)
            r_l.append(mod.order_parameter(stB))
            std_l.append(mod.omega_std(stB))
            sync_l.append(mod.sync_count(stB))
            lam_l.append(mod.lambda_mean(stB, v_hub))
            tau_l.append(stB.get('tau_var', 0))
            J_l.append(stB.get('J_cost', 0))
            K_bus_l.append(stB.get('K_Q3_BUS', mod.K_Q3_BUS_0))
            K_ind_l.append(stB.get('K_Q3_IND', mod.K_Q3_IND_0))
            E_cyc_l.append(stB.get('E_cicle', 0)/1e6)

    # Restaura
    for k2, v2 in orig.items():
        setattr(mod, k2, v2)

    wu = len(PA_l) // 5
    PA  = np.array(PA_l[wu:]);   PB  = np.array(PB_l[wu:])
    PQ3 = np.array(PQ3_l[wu:]);  PBt = PB + PQ3
    dP  = np.mean(PBt - PA)

    # Mètriques + robustesa §5.1: winsorització de J
    J_arr = np.array(J_l[wu:])
    J_arr_w = np.clip(J_arr,
                      np.percentile(J_arr, 5) if len(J_arr)>5 else J_arr.min(),
                      np.percentile(J_arr, 95) if len(J_arr)>5 else J_arr.max())

    mask = np.array(E_cyc_l[wu:]) != 0
    E_m  = float(np.mean(np.array(E_cyc_l[wu:])[mask])) if mask.sum()>0 else 0.0

    return {
        'seed':      seed,
        'dP_mean':   float(dP),
        'dP_per_wt': float(dP / mod.N_ZZ * 1000),
        'P_A_mean':  float(np.mean(PA)),
        'P_B_mean':  float(np.mean(PBt)),
        'P_std':     float(np.std(PBt)),
        'r_mean':    float(np.mean(r_l[wu:])),
        'sync_mean': float(np.mean(sync_l[wu:])),
        'std_om':    float(np.mean(std_l[wu:])),
        'lambda':    float(np.mean(lam_l[wu:])),
        'tau_var':   float(np.mean(tau_l[wu:])),
        'J_mean':    float(np.mean(J_arr_w)),
        'J_final':   float(J_arr_w[-1]) if len(J_arr_w)>0 else 0,
        'K_bus_fin': float(K_bus_l[-1]) if K_bus_l else mod.K_Q3_BUS_0,
        'K_ind_fin': float(K_ind_l[-1]) if K_ind_l else mod.K_Q3_IND_0,
        'E_cicle':   E_m,
        'dP_pos_pct':float(100*np.mean(PBt > PA)),
        'policy_x':  list(stB.get('_policy_x', [0]*5)),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 1. BASELINE MULTI-SEED
# ═══════════════════════════════════════════════════════════════════════════════

def run_baseline():
    print("\n" + "─"*60)
    print("  1. BASELINE MULTI-SEED (10 seeds, T=300s)")
    print("─"*60)
    for s in range(N_SEEDS_BASE):
        t0 = time.time()
        r  = run_short(seed=s, T=T_SHORT)
        dt = time.time()-t0
        results['baseline'].append(r)
        print(f"    seed={s:2d}  ΔP={r['dP_mean']:+.3f} MW  "
              f"r={r['r_mean']:.3f}  sync={r['sync_mean']:.1f}  "
              f"ΔP>0: {r['dP_pos_pct']:.0f}%  ({dt:.1f}s)")

    dPs = [r['dP_mean'] for r in results['baseline']]
    print(f"\n    → ΔP  mitja={np.mean(dPs):+.3f} ± {np.std(dPs):.3f} MW")
    print(f"    → Runs amb ΔP>0: {sum(d>0 for d in dPs)}/{N_SEEDS_BASE}")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ESTABILITAT POLÍTICA
# ═══════════════════════════════════════════════════════════════════════════════

def run_policy_stability():
    print("\n" + "─"*60)
    print(f"  2. ESTABILITAT POLÍTICA (ε={EPS_REDUCED}, η={ETA_REDUCED})")
    print("─"*60)
    # Modifica paràmetres temporalment
    orig_eps = mod.EPS_GRAD;  orig_eta = mod.ETA_POLICY
    mod.EPS_GRAD = EPS_REDUCED;  mod.ETA_POLICY = ETA_REDUCED

    for s in range(5):
        r = run_short(seed=s, T=T_SHORT, eps=EPS_REDUCED, eta=ETA_REDUCED)
        results['policy'].append(r)
        print(f"    seed={s}  ΔP={r['dP_mean']:+.3f} MW  "
              f"J_final={r['J_final']:+.2f}  "
              f"K_bus={r['K_bus_fin']:.4f}  "
              f"policy_x={[f'{v:.2f}' for v in r['policy_x']]}")

    mod.EPS_GRAD = orig_eps;  mod.ETA_POLICY = orig_eta

    dPs = [r['dP_mean'] for r in results['policy']]
    print(f"\n    → ΔP  mitja={np.mean(dPs):+.3f} ± {np.std(dPs):.3f} MW")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. MONTE CARLO VENT
# ═══════════════════════════════════════════════════════════════════════════════

def run_montecarlo():
    print("\n" + "─"*60)
    print(f"  3. MONTE CARLO VENT ({N_SEEDS_MC} runs)")
    print("─"*60)
    rng = np.random.default_rng(99)
    # Escenaris de vent: turbulència i ràfega variables
    turbs  = rng.uniform(0.06, 0.20, N_SEEDS_MC)
    gusts  = rng.uniform(11.0, 17.0, N_SEEDS_MC)
    seeds  = rng.integers(0, 100, N_SEEDS_MC)

    for i in range(N_SEEDS_MC):
        r = run_short(seed=int(seeds[i]), T=T_SHORT,
                      wind_turb=turbs[i], wind_gust_v=gusts[i])
        r['turb']  = float(turbs[i])
        r['gust']  = float(gusts[i])
        results['montecarlo'].append(r)
        if i % 10 == 0:
            dPs_so_far = [x['dP_mean'] for x in results['montecarlo']]
            pct_pos    = 100*sum(d>0 for d in dPs_so_far)/len(dPs_so_far)
            print(f"    {i+1}/{N_SEEDS_MC}  "
                  f"ΔP_mitja={np.mean(dPs_so_far):+.3f} MW  "
                  f"ΔP>0: {pct_pos:.0f}%")

    dPs  = [r['dP_mean']    for r in results['montecarlo']]
    pcts = [r['dP_pos_pct'] for r in results['montecarlo']]
    print(f"\n    → ΔP  mitja={np.mean(dPs):+.3f} ± {np.std(dPs):.3f} MW")
    print(f"    → Runs amb ΔP>0: {sum(d>0 for d in dPs)}/{N_SEEDS_MC} "
          f"({100*sum(d>0 for d in dPs)/N_SEEDS_MC:.0f}%)")
    print(f"    → τ_var mitja: {np.mean([r['tau_var'] for r in results['montecarlo']]):.2e}")


# ═══════════════════════════════════════════════════════════════════════════════
# 4. ABLATION STUDY
# ═══════════════════════════════════════════════════════════════════════════════

def run_ablation():
    print("\n" + "─"*60)
    print("  4. ABLATION STUDY (5 seeds × 4 configuracions)")
    print("─"*60)

    configs = {
        'full':         dict(),
        'no_fft':       dict(disable_fft=True),
        'no_tsr':       dict(disable_tsr=True),
        'no_kind':      dict(disable_kind=True),
        'no_load_pen':  dict(disable_load=True),
    }

    for cfg_name, kwargs in configs.items():
        runs = []
        for s in range(5):
            r = run_short(seed=s, T=T_SHORT, **kwargs)
            runs.append(r)
        dPs   = [r['dP_mean'] for r in runs]
        taus  = [r['tau_var'] for r in runs]
        stds  = [r['P_std']   for r in runs]
        results['ablation'][cfg_name] = {
            'dP_mean': float(np.mean(dPs)),
            'dP_std':  float(np.std(dPs)),
            'tau_mean':float(np.mean(taus)),
            'P_std':   float(np.mean(stds)),
            'runs':    runs,
        }
        print(f"    {cfg_name:<18}  "
              f"ΔP={np.mean(dPs):+.3f}±{np.std(dPs):.3f} MW  "
              f"std(P)={np.mean(stds):.1f}  "
              f"τ_var={np.mean(taus):.2e}")

    # Contribució marginal de cada millora
    base_dP = results['ablation']['full']['dP_mean']
    print(f"\n    Contribució marginal vs 'full':")
    for cfg_name in ['no_fft','no_tsr','no_kind','no_load_pen']:
        d = results['ablation'][cfg_name]['dP_mean']
        print(f"      Sense {cfg_name[3:]:10s}: ΔΔP = {d-base_dP:+.3f} MW")


# ═══════════════════════════════════════════════════════════════════════════════
# 5. SURROGATE + CERCA BAYESIANA SIMPLE
# ═══════════════════════════════════════════════════════════════════════════════

def run_surrogate():
    print("\n" + "─"*60)
    print("  5. SURROGATE OFFLINE + CERCA BAYESIANA")
    print("─"*60)

    # Recull punts (x, J) de les runs anteriors
    X_data, J_data = [], []
    for r in results['baseline'] + results['policy']:
        px = r['policy_x']
        if len(px) == 5 and r['J_mean'] != 0:
            X_data.append(px)
            J_data.append(-r['dP_mean'])  # volem minimitzar J ≈ -ΔP

    # Afegim grid de punts extra (ràpids)
    print("    Generant punts addicionals per surrogate...")
    rng2 = np.random.default_rng(77)
    N_EXTRA = 20
    p_grid = [
        [rng2.uniform(0.2, 0.9) for _ in range(5)]
        for _ in range(N_EXTRA)
    ]
    for xi in p_grid:
        pw, me, kb, nm, ki = mod.decode_policy(np.array(xi))
        # Simulació curta (120s) per estimar J
        orig_k = mod.K_Q3_BUS_0;  orig_i = mod.K_Q3_IND_0
        mod.K_Q3_BUS_0 = kb;       mod.K_Q3_IND_0 = ki
        r = run_short(seed=3, T=120.0)
        mod.K_Q3_BUS_0 = orig_k;   mod.K_Q3_IND_0 = orig_i
        X_data.append(xi)
        J_data.append(-r['dP_mean'])

    X_data = np.array(X_data);  J_data = np.array(J_data)
    print(f"    Punts recollits: {len(X_data)}  "
          f"ΔP rang: [{-J_data.max():.2f}, {-J_data.min():.2f}] MW")

    # Surrogate: Kernel Ridge Regression (manual, sense sklearn)
    # K(x,x') = exp(-||x-x'||² / 2σ²)
    sigma   = 0.3
    alpha_r = 1e-3
    def kernel(A, B):
        # (n,5) × (m,5) → (n,m)
        diff = A[:,None,:] - B[None,:,:]  # (n,m,5)
        return np.exp(-np.sum(diff**2, axis=-1) / (2*sigma**2))

    K_train = kernel(X_data, X_data) + alpha_r*np.eye(len(X_data))
    try:
        alpha_coef = np.linalg.solve(K_train, J_data)
    except np.linalg.LinAlgError:
        alpha_coef = np.linalg.lstsq(K_train, J_data, rcond=None)[0]

    def predict(X_new):
        K_new = kernel(X_new, X_data)
        return K_new @ alpha_coef

    # Cerca Bayesiana simple: random search sobre el surrogate
    rng3    = np.random.default_rng(123)
    N_CAND  = 500
    X_cand  = rng3.uniform(0.05, 0.95, (N_CAND, 5))
    J_pred  = predict(X_cand)
    # UCB: exploració + explotació
    ucb_w   = 0.1
    # Estimem variança (diagonal del kernel)
    K_cand  = kernel(X_cand, X_data)
    K_cc_diag = np.ones(N_CAND)  # aprox. k(x,x)=1
    K_inv_K = K_cand @ np.linalg.pinv(K_train) @ K_cand.T
    var_pred = np.maximum(0, K_cc_diag - np.diag(K_inv_K))
    ucb     = J_pred - ucb_w * np.sqrt(var_pred)  # minimitzem J
    best_idx = np.argmin(ucb)
    x_best   = X_cand[best_idx]

    pw,me,kb,nm,ki = mod.decode_policy(x_best)
    print(f"\n    Millors paràmetres trobats pel surrogate (UCB):")
    print(f"      p={pw:.2f}  m={me/1e3:.1f}t  K_bus={kb:.4f}  N=3:{nm}  K_ind={ki:.4f}")

    # Valida el millor punt amb una simulació real
    orig_k = mod.K_Q3_BUS_0;  orig_i = mod.K_Q3_IND_0
    mod.K_Q3_BUS_0 = kb;       mod.K_Q3_IND_0 = ki
    r_best = run_short(seed=7, T=T_SHORT)
    mod.K_Q3_BUS_0 = orig_k;   mod.K_Q3_IND_0 = orig_i

    print(f"      Validació real: ΔP={r_best['dP_mean']:+.3f} MW  "
          f"r={r_best['r_mean']:.3f}  τ_var={r_best['tau_var']:.2e}")

    # Guarda la millor política
    best_policy = {
        'x':         list(x_best),
        'p_win':     float(pw),
        'm_exc_t':   float(me/1e3),
        'K_bus':     float(kb),
        'N_mode':    int(nm),
        'K_ind':     float(ki),
        'dP_pred':   float(-J_pred[best_idx]),
        'dP_real':   float(r_best['dP_mean']),
    }
    results['surrogate'] = {
        'best_policy':   best_policy,
        'X_train_n':     len(X_data),
        'J_min_data':    float(J_data.min()),
        'J_max_data':    float(J_data.max()),
        'validation_r':  r_best,
    }
    with open('/mnt/user-data/outputs/best_policy.txt', 'w') as f:
        f.write("MILLORS PARÀMETRES — Surrogate + UCB\n")
        f.write("="*40 + "\n")
        for k2,v2 in best_policy.items():
            f.write(f"  {k2}: {v2}\n")
    print(f"\n    Guardat: best_policy.txt")
    return best_policy


# ═══════════════════════════════════════════════════════════════════════════════
# GRÀFICS
# ═══════════════════════════════════════════════════════════════════════════════

def plot_summary():
    fig = plt.figure(figsize=(20, 16))
    fig.patch.set_facecolor('#0d1117')
    gs  = gridspec.GridSpec(4, 4, figure=fig, hspace=0.52, wspace=0.40)

    def mkax(r,c): return fig.add_subplot(gs[r,c])
    def style(a, title, xl='', yl=''):
        a.set_facecolor(BG); a.set_title(title,color=TXT,fontsize=9,pad=4)
        if xl: a.set_xlabel(xl,color=TXT,fontsize=8)
        if yl: a.set_ylabel(yl,color=TXT,fontsize=8)
        a.tick_params(colors=TXT,labelsize=7)
        for sp in a.spines.values(): sp.set_color('#30363d')
        a.grid(True,color='#21262d',linewidth=0.4)

    # 1. Baseline: distribució ΔP per seed
    a1 = mkax(0,0)
    dPs_b = [r['dP_mean'] for r in results['baseline']]
    seeds = list(range(len(dPs_b)))
    cols  = [CG if d>0 else CR for d in dPs_b]
    a1.bar(seeds, dPs_b, color=cols, width=0.7, alpha=0.8)
    a1.axhline(0,color='white',lw=0.8,ls='--',alpha=0.5)
    a1.axhline(np.mean(dPs_b),color=CQ,lw=1.5,ls='--',
               label=f'mitj={np.mean(dPs_b):+.2f}')
    a1.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a1,'1. Baseline multi-seed','seed','ΔP [MW]')

    # 2. Política: J_final per seed
    a2 = mkax(0,1)
    if results['policy']:
        dPs_p = [r['dP_mean'] for r in results['policy']]
        J_p   = [r['J_final'] for r in results['policy']]
        a2.scatter(dPs_p, J_p, color=CV, s=50, zorder=5)
        for i,(x,y) in enumerate(zip(dPs_p,J_p)):
            a2.annotate(f's{i}', (x,y), fontsize=7, color=TXT,
                        xytext=(3,3), textcoords='offset points')
        a2.axvline(0,color='white',lw=0.5,ls='--',alpha=0.4)
    style(a2,'2. ΔP vs J_final (política ε=0.04)','ΔP [MW]','J')

    # 3. Monte Carlo: histograma ΔP
    a3 = mkax(0,2)
    dPs_mc = [r['dP_mean'] for r in results['montecarlo']]
    if dPs_mc:
        bins = np.linspace(min(dPs_mc)-1, max(dPs_mc)+1, 20)
        a3.hist([d for d in dPs_mc if d>=0], bins=bins,
                color=CG, alpha=0.7, label=f'ΔP≥0 ({sum(d>=0 for d in dPs_mc)})')
        a3.hist([d for d in dPs_mc if d<0], bins=bins,
                color=CR, alpha=0.7, label=f'ΔP<0 ({sum(d<0 for d in dPs_mc)})')
        a3.axvline(np.mean(dPs_mc),color=CQ,lw=1.5,ls='--',
                   label=f'mitj={np.mean(dPs_mc):+.2f}')
        pct_pos = 100*sum(d>=0 for d in dPs_mc)/len(dPs_mc)
        a3.text(0.97,0.95,f'{pct_pos:.0f}% positiu',
                transform=a3.transAxes,ha='right',va='top',
                color=CG,fontsize=9)
        a3.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a3,'3. Monte Carlo ΔP (50 runs)','ΔP [MW]','N runs')

    # 4. Monte Carlo: ΔP vs turbulència
    a4 = mkax(0,3)
    if results['montecarlo']:
        turbs = [r['turb'] for r in results['montecarlo']]
        dPs_  = [r['dP_mean'] for r in results['montecarlo']]
        cols4 = [CG if d>0 else CR for d in dPs_]
        a4.scatter(turbs, dPs_, c=cols4, s=20, alpha=0.7)
        # Tendència
        if len(turbs)>3:
            z = np.polyfit(turbs, dPs_, 1)
            px= np.linspace(min(turbs), max(turbs), 50)
            a4.plot(px, np.polyval(z,px), color=CQ, lw=1.5, ls='--',
                    alpha=0.8, label=f'tendència')
        a4.axhline(0,color='white',lw=0.5,ls='--',alpha=0.4)
        a4.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a4,'MC: ΔP vs turbulència','turbulència','ΔP [MW]')

    # 5. Ablation: barres comparatives
    a5 = mkax(1,slice(0,2))
    if results['ablation']:
        cfg_names = list(results['ablation'].keys())
        dP_vals   = [results['ablation'][c]['dP_mean'] for c in cfg_names]
        dP_errs   = [results['ablation'][c]['dP_std']  for c in cfg_names]
        xs        = range(len(cfg_names))
        cols5     = [CG if d>0 else CR for d in dP_vals]
        bars = a5.bar(xs, dP_vals, color=cols5, width=0.6, alpha=0.8,
                      yerr=dP_errs, capsize=4,
                      error_kw={'ecolor':TXT,'elinewidth':1})
        a5.set_xticks(xs)
        a5.set_xticklabels(cfg_names, color=TXT, fontsize=7, rotation=12)
        a5.axhline(0,color='white',lw=0.8,ls='--',alpha=0.5)
        # Anotació contribució marginal
        full_dP = results['ablation'].get('full',{}).get('dP_mean',0)
        for i,v in enumerate(dP_vals):
            a5.text(i, v+0.1 if v>0 else v-0.2,
                    f'{v-full_dP:+.2f}' if i>0 else 'ref',
                    ha='center', color=TXT, fontsize=7)
    style(a5,'4. Ablation study — contribució marginal','','ΔP [MW]')

    # 6. Ablation: std(P) comparatiu
    a6 = mkax(1,2)
    if results['ablation']:
        cfg_names = list(results['ablation'].keys())
        std_vals  = [results['ablation'][c]['P_std'] for c in cfg_names]
        a6.bar(range(len(cfg_names)), std_vals, color=CM, width=0.6, alpha=0.8)
        a6.set_xticks(range(len(cfg_names)))
        a6.set_xticklabels(cfg_names, color=TXT, fontsize=7, rotation=12)
    style(a6,'Ablation: std(P)','','std(P) [MW]')

    # 7. Ablation: tau_var
    a7 = mkax(1,3)
    if results['ablation']:
        cfg_names = list(results['ablation'].keys())
        tau_vals  = [results['ablation'][c]['tau_mean'] for c in cfg_names]
        cols7     = [CR if t > mod.TAU_VAR_MAX else CG for t in tau_vals]
        a7.bar(range(len(cfg_names)), tau_vals, color=cols7, width=0.6, alpha=0.8)
        a7.set_xticks(range(len(cfg_names)))
        a7.set_xticklabels(cfg_names, color=TXT, fontsize=7, rotation=12)
        a7.axhline(mod.TAU_VAR_MAX,color='white',lw=0.8,ls='--',alpha=0.5,
                   label='límit struct.')
        a7.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a7,'Ablation: var(τ)','','N²m²')

    # 8. Surrogate: mapa 2D (p_win vs m_exc)
    a8 = mkax(2,slice(0,2))
    if results.get('surrogate') and len(results['surrogate'].get('X_train_n',0))>0:
        # Reconstruïm dades de X_data / J_data
        X_all = [r['policy_x'] for r in results['baseline']+results['policy']
                 if len(r['policy_x'])==5 and r['J_mean']!=0]
        J_all = [-r['dP_mean'] for r in results['baseline']+results['policy']
                 if len(r['policy_x'])==5 and r['J_mean']!=0]
        if X_all:
            Xa = np.array(X_all)
            Ja = np.array(J_all)
            sc = a8.scatter(Xa[:,0], Xa[:,1], c=-Ja,
                            cmap='RdYlGn', s=40, alpha=0.8)
            plt.colorbar(sc, ax=a8, fraction=0.046, pad=0.04).ax.tick_params(
                colors=TXT, labelsize=7)
            # Marca millor
            bp = results['surrogate'].get('best_policy',{})
            if bp:
                px_b = (bp['p_win']-mod.P_WIN_MIN)/(mod.P_WIN_MAX-mod.P_WIN_MIN)
                py_b = bp['m_exc_t']*1e3/mod.M_EXC_MAX
                a8.scatter([px_b],[py_b],color='white',s=120,
                           marker='*',zorder=10,label='millor')
                a8.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a8,'5. Surrogate: ΔP(p_win, m_exc)','p_win (norm)','m_exc (norm)')

    # 9. ΔP vs gust (MC)
    a9 = mkax(2,2)
    if results['montecarlo']:
        gusts_ = [r['gust'] for r in results['montecarlo']]
        dPs_   = [r['dP_mean'] for r in results['montecarlo']]
        cols9  = [CG if d>0 else CR for d in dPs_]
        a9.scatter(gusts_, dPs_, c=cols9, s=20, alpha=0.7)
        a9.axhline(0,color='white',lw=0.5,ls='--',alpha=0.4)
    style(a9,'MC: ΔP vs gust','v_gust [m/s]','ΔP [MW]')

    # 10. r_mean per escenari
    a10 = mkax(2,3)
    escenaris = ['baseline','policy']
    for ei,esc in enumerate(escenaris):
        if results[esc]:
            rs = [r['r_mean'] for r in results[esc]]
            a10.bar([ei*0.4+j*0.08 for j in range(len(rs))],
                    rs, width=0.07, alpha=0.7,
                    color=CA if esc=='baseline' else CV,
                    label=esc)
    a10.set_ylim(0,0.5)
    a10.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a10,'Ordre r per escenari','','r')

    # 11–12. Resum text
    a11 = mkax(3, slice(0,4))
    a11.set_facecolor(BG);  a11.axis('off')

    dPs_b = [r['dP_mean'] for r in results['baseline']]
    dPs_p = [r['dP_mean'] for r in results['policy']]
    dPs_m = [r['dP_mean'] for r in results['montecarlo']]
    bp    = results['surrogate'].get('best_policy', {})

    lines = [
        ("EXPERIMENT v4 — RESUM EXECUTIU", CQ),
        (f"  1. Baseline  ΔP={np.mean(dPs_b):+.3f}±{np.std(dPs_b):.3f} MW  "
         f"ΔP>0: {sum(d>0 for d in dPs_b)}/{len(dPs_b)}", TXT),
        (f"  2. Política  ΔP={np.mean(dPs_p):+.3f}±{np.std(dPs_p):.3f} MW  "
         f"(ε={EPS_REDUCED} η={ETA_REDUCED})" if dPs_p else "  2. Política  —", TXT),
        (f"  3. MC vent   ΔP={np.mean(dPs_m):+.3f}±{np.std(dPs_m):.3f} MW  "
         f"ΔP>0: {sum(d>0 for d in dPs_m)}/{len(dPs_m)} "
         f"({100*sum(d>0 for d in dPs_m)/max(1,len(dPs_m)):.0f}%)"
         if dPs_m else "  3. MC  —", TXT),
        (f"  4. Ablation  contribucions marginals: "
         + "  ".join(f"{k[3:]}:{v['dP_mean']-results['ablation'].get('full',{'dP_mean':0})['dP_mean']:+.2f}MW"
                     for k,v in results['ablation'].items() if k!='full')
         if results['ablation'] else "  4. Ablation —", TXT),
        (f"  5. Surrogate  p={bp.get('p_win',0):.2f}  "
         f"m={bp.get('m_exc_t',0):.1f}t  "
         f"K_bus={bp.get('K_bus',0):.4f}  "
         f"ΔP_pred={bp.get('dP_pred',0):+.3f} → "
         f"ΔP_real={bp.get('dP_real',0):+.3f} MW"
         if bp else "  5. Surrogate —", CG if bp.get('dP_real',0)>0 else CR),
    ]
    for li, (ln, col) in enumerate(lines):
        a11.text(0.01, 0.92-li*0.16, ln, color=col, fontsize=9,
                 fontfamily='monospace', transform=a11.transAxes)

    fig.suptitle(
        'Experiment v4 — Baseline · Política · Monte Carlo · Ablation · Surrogate  |  '
        f'Víctor Manzanares Alberola — EPSA/UPV Alcoi',
        color=TXT, fontsize=9, y=0.999)

    out = '/mnt/user-data/outputs/experiment_summary.png'
    plt.savefig(out, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    print(f"\n  Gràfic: {out}")
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    t0_total = time.time()
    print("═"*70)
    print("  experiment_v4.py")
    print("  Experiment controlat: Baseline·Política·MC·Ablation·Surrogate")
    print("  Víctor Manzanares Alberola — EPSA/UPV Alcoi")
    print("═"*70)

    run_baseline()
    run_policy_stability()
    run_montecarlo()
    run_ablation()
    best = run_surrogate()

    # Guarda JSON
    import json
    # Converteix arrays numpy a llistes per a JSON
    def to_json(obj):
        if isinstance(obj, dict):
            return {k: to_json(v) for k,v in obj.items()}
        elif isinstance(obj, list):
            return [to_json(v) for v in obj]
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        return obj

    with open('/mnt/user-data/outputs/experiment_results.json','w') as f:
        json.dump(to_json(results), f, indent=2)

    print("\n  Generant gràfics resum...")
    plot_summary()

    import shutil
    shutil.copy(__file__, '/mnt/user-data/outputs/experiment_v4.py')

    print(f"\n  ✓ Experiment completat en {(time.time()-t0_total)/60:.1f} min")
    print(f"  Fitxers: experiment_results.json · experiment_summary.png · best_policy.txt")
