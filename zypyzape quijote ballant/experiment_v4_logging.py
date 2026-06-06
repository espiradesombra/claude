#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
experiment_v4_logging.py
═══════════════════════════════════════════════════════════════════════════════
Runner d'experiments per model_44_Quijote_zypyzape.py v4

Característiques:
  - Logging JSON per run (metadades + mètriques)
  - Checkpointing automàtic cada N runs
  - Paral·lelització via multiprocessing.Pool
  - 5 fases: Baseline · Política · Monte Carlo · Ablation · Surrogate
  - Recuperació automàtica si la sessió cau

Ús:
  python3 experiment_v4_logging.py
  python3 experiment_v4_logging.py --quick     # versió ràpida (T=120s)
  python3 experiment_v4_logging.py --workers 4

Autor: Víctor Manzanares Alberola — EPSA/UPV Alcoi
═══════════════════════════════════════════════════════════════════════════════
"""

import os, sys, time, json, traceback, argparse
from datetime import datetime
from multiprocessing import Pool, cpu_count
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ─── Carrega el model des del mateix directori ────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_MODEL_PATH = os.path.join(_HERE, "model_44_Quijote_zypyzape.py")
if not os.path.exists(_MODEL_PATH):
    raise FileNotFoundError(f"Model no trobat: {_MODEL_PATH}")

import importlib.util
_spec = importlib.util.spec_from_file_location("model", _MODEL_PATH)
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURACIÓ
# ═══════════════════════════════════════════════════════════════════════════════

OUT_DIR            = os.path.join(_HERE, "outputs", "experiments")
os.makedirs(OUT_DIR, exist_ok=True)

CHECKPOINT_FILE    = os.path.join(OUT_DIR, "checkpoint_latest.json")
RESULTS_FILE       = os.path.join(OUT_DIR, "experiment_results.json")
BEST_POLICY_FILE   = os.path.join(OUT_DIR, "best_policy.txt")
SUMMARY_PNG        = os.path.join(OUT_DIR, "experiment_summary.png")

CHECKPOINT_EVERY   = 10      # guarda checkpoint cada N runs

# Paràmetres per defecte (sobreescribibles per CLI)
N_SEEDS_BASE       = 10
N_SEEDS_MC         = 50
T_SHORT            = 300.0
T_FULL             = 600.0
N_WORKERS          = max(1, min(cpu_count(), 6))

EPS_REDUCED        = 0.04
ETA_REDUCED        = 0.03
N_REPEAT_PERT      = 3

# Paleta
BG='#161b22'; TXT='#c9d1d9'
CA='#378ADD'; CB='#D85A30'; CV='#5DCAA5'; CQ='#FAC775'
CG='#97C459'; CR='#ED93B1'; CM='#AFA9EC'


# ═══════════════════════════════════════════════════════════════════════════════
# UTILITATS
# ═══════════════════════════════════════════════════════════════════════════════

def ts():
    return datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")

def _json_default(o):
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, (np.integer, np.floating)): return float(o)
    return str(o)

def safe_write(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=2, default=_json_default)
    os.replace(tmp, path)

def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        try:
            with open(CHECKPOINT_FILE) as f:
                return json.load(f)
        except Exception:
            return None
    return None

def save_checkpoint(data):
    safe_write(CHECKPOINT_FILE, data)

def consolidate():
    """Llegeix tots els JSONs de run i consolida en un únic fitxer."""
    runs = []
    for fn in sorted(os.listdir(OUT_DIR)):
        if fn.startswith("run_") and fn.endswith(".json"):
            try:
                with open(os.path.join(OUT_DIR, fn)) as f:
                    runs.append(json.load(f))
            except Exception:
                pass
    safe_write(RESULTS_FILE, {"runs": runs, "generated_at": ts()})
    if len(runs) % CHECKPOINT_EVERY == 0:
        save_checkpoint({"n_runs": len(runs), "updated": ts()})
    return runs


# ═══════════════════════════════════════════════════════════════════════════════
# NUCLi D'EXECUCIÓ (compatible amb multiprocessing)
# ═══════════════════════════════════════════════════════════════════════════════

def _run_single(args):
    """
    Executa un sol run i desa el JSON de resultats.
    args: dict amb claus phase, seed, T, overrides, extra
    """
    t0    = time.time()
    seed  = int(args.get("seed", 0))
    phase = args.get("phase", "baseline")
    T     = float(args.get("T", T_SHORT))
    ov    = args.get("overrides", {}) or {}
    extra = args.get("extra", {}) or {}

    run_id   = f"{phase}_{ts()}_s{seed}"
    out_file = os.path.join(OUT_DIR, f"run_{run_id}.json")

    result = {"run_id": run_id, "phase": phase, "seed": seed,
              "overrides": ov, "extra": extra,
              "start": ts(), "status": "failed",
              "error": None, "metrics": None, "runtime_s": None}

    # Guarda i aplica overrides al mòdul
    orig = {}
    try:
        for k, v in ov.items():
            if hasattr(mod, k):
                orig[k] = getattr(mod, k)
                setattr(mod, k, v)

        metrics = _simulate(seed=seed, T=T,
                            disable_fft=extra.get("disable_fft", False),
                            disable_tsr=extra.get("disable_tsr", False),
                            disable_kind=extra.get("disable_kind", False),
                            disable_load=extra.get("disable_load", False))
        result["metrics"] = metrics
        result["status"]  = "ok"
    except Exception as e:
        result["error"] = {"msg": str(e), "tb": traceback.format_exc()}
    finally:
        for k, v in orig.items():
            setattr(mod, k, v)
        result["runtime_s"] = round(time.time() - t0, 2)
        result["end"]        = ts()
        try:
            safe_write(out_file, result)
        except Exception:
            pass
    return result


def _simulate(seed=7, T=T_SHORT,
              disable_fft=False, disable_tsr=False,
              disable_kind=False, disable_load=False):
    """Executa escenaris A+B i retorna mètriques."""
    DT    = mod.DT
    STEPS = int(T / DT)
    t_arr = np.arange(STEPS) * DT
    v_arr = mod.wind_profile(t_arr)

    stA = mod.init_state(seed)
    stB = mod.init_state(seed)
    stB['theta']     = stA['theta'].copy()
    stB['omega']     = stA['omega'].copy()
    stB['omega_nat'] = stA['omega_nat'].copy()

    PA_l=[]; PB_l=[]; PQ3_l=[]
    r_l=[]; std_l=[]; sync_l=[]; lam_l=[]
    tau_l=[]; J_l=[]; Kb_l=[]; Ki_l=[]; Ec_l=[]

    for k in range(STEPS):
        v    = v_arr[k]; t_now = k * DT
        PA_n = mod.power_MW(stA, v)
        stA  = mod.step_ZZ_only(stA, v)

        if disable_tsr:
            stB['_K_bus_tsr'] = mod.K_BUS_BASE
        if disable_kind:
            stB['K_Q3_IND']   = mod.K_Q3_IND_0
        if disable_load:
            stB['tau_var']    = 0.0

        stB = mod.step_ZZ_Q3(stB, v, PA_n, t_now)

        if disable_fft:
            stB['res_mode'] = 12

        if k % 4 == 0:
            PA_l.append(mod.power_MW(stA, v))
            pB  = mod.power_MW(stB, v)
            pQ3 = mod.power_Q3_MW(stB, v)
            PB_l.append(pB); PQ3_l.append(pQ3)
            r_l.append(mod.order_parameter(stB))
            std_l.append(mod.omega_std(stB))
            sync_l.append(mod.sync_count(stB))
            lam_l.append(mod.lambda_mean(stB, v))
            tau_l.append(stB.get('tau_var', 0))
            J_l.append(stB.get('J_cost', 0))
            Kb_l.append(stB.get('K_Q3_BUS', mod.K_Q3_BUS_0))
            Ki_l.append(stB.get('K_Q3_IND', mod.K_Q3_IND_0))
            Ec_l.append(stB.get('E_cicle', 0) / 1e6)

    wu  = len(PA_l) // 5
    PA  = np.array(PA_l[wu:])
    PBt = np.array(PB_l[wu:]) + np.array(PQ3_l[wu:])
    dP  = float(np.mean(PBt - PA)) if PA.size > 0 else 0.0

    J   = np.array(J_l[wu:])
    if J.size > 4:
        J = np.clip(J, np.percentile(J, 5), np.percentile(J, 95))

    mask = np.array(Ec_l[wu:]) != 0
    E_m  = float(np.mean(np.array(Ec_l[wu:])[mask])) if mask.sum() > 0 else 0.0

    return {
        "seed":        seed,
        "dP_mean":     dP,
        "dP_per_kw":   dP / mod.N_ZZ * 1000,
        "P_A_mean":    float(np.mean(PA)) if PA.size > 0 else 0.0,
        "P_B_mean":    float(np.mean(PBt)) if PBt.size > 0 else 0.0,
        "P_std":       float(np.std(PBt)) if PBt.size > 0 else 0.0,
        "r_mean":      float(np.mean(r_l[wu:])) if len(r_l) > wu else 0.0,
        "sync_mean":   float(np.mean(sync_l[wu:])) if len(sync_l) > wu else 0.0,
        "std_om":      float(np.mean(std_l[wu:])) if len(std_l) > wu else 0.0,
        "lambda":      float(np.mean(lam_l[wu:])) if len(lam_l) > wu else 0.0,
        "tau_var":     float(np.mean(tau_l[wu:])) if len(tau_l) > wu else 0.0,
        "J_mean":      float(np.mean(J)) if J.size > 0 else 0.0,
        "J_final":     float(J[-1]) if J.size > 0 else 0.0,
        "K_bus_fin":   float(Kb_l[-1]) if Kb_l else mod.K_Q3_BUS_0,
        "K_ind_fin":   float(Ki_l[-1]) if Ki_l else mod.K_Q3_IND_0,
        "E_cicle":     E_m,
        "dP_pos_pct":  float(100 * np.mean(PBt > PA)) if PA.size > 0 else 0.0,
        "policy_x":    list(stB.get("_policy_x", [0]*5)),
        "via1_MJ":     stB.get("via1_acc", 0) / 1e6,
        "via2_MJ":     stB.get("via2_acc", 0) / 1e6,
        "via3_MJ":     stB.get("via3_acc", 0) / 1e6,
        "shield_MJ":   stB.get("shield_acc", 0) / 1e6,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# RUNNER PARAL·LEL
# ═══════════════════════════════════════════════════════════════════════════════

def _run_parallel(arg_list, workers=1):
    results = []
    if workers <= 1:
        for a in arg_list:
            results.append(_run_single(a))
            consolidate()
    else:
        with Pool(processes=min(workers, len(arg_list))) as pool:
            for r in pool.imap_unordered(_run_single, arg_list):
                results.append(r)
                consolidate()
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# FASES D'EXPERIMENT
# ═══════════════════════════════════════════════════════════════════════════════

def phase_baseline(n=N_SEEDS_BASE, T=T_SHORT, workers=N_WORKERS):
    print(f"\n── 1. BASELINE ({n} seeds, T={T:.0f}s) ──")
    args = [{"phase":"baseline","seed":s,"T":T,"overrides":{}} for s in range(n)]
    res  = _run_parallel(args, workers)
    dPs  = [r["metrics"]["dP_mean"] for r in res if r["status"]=="ok"]
    print(f"   ΔP = {np.mean(dPs):+.3f} ± {np.std(dPs):.3f} MW  "
          f"ΔP>0: {sum(d>0 for d in dPs)}/{len(dPs)}")
    return res


def phase_policy(n=5, T=T_SHORT, workers=1):
    print(f"\n── 2. POLÍTICA ESTABILITAT (ε={EPS_REDUCED}, η={ETA_REDUCED}) ──")
    args = [{"phase":"policy","seed":s,"T":T,
             "overrides":{"EPS_GRAD":EPS_REDUCED,"ETA_POLICY":ETA_REDUCED}}
            for s in range(n)]
    res  = _run_parallel(args, workers)
    dPs  = [r["metrics"]["dP_mean"] for r in res if r["status"]=="ok"]
    Js   = [r["metrics"]["J_final"] for r in res if r["status"]=="ok"]
    print(f"   ΔP = {np.mean(dPs):+.3f} MW  J_final mitj = {np.mean(Js):+.2f}")
    return res


def phase_montecarlo(n=N_SEEDS_MC, T=T_SHORT, workers=N_WORKERS):
    print(f"\n── 3. MONTE CARLO ({n} runs) ──")
    rng   = np.random.default_rng(99)
    turbs = rng.uniform(0.06, 0.20, n)
    gusts = rng.uniform(11.0, 17.0, n)
    seeds = rng.integers(0, 10000, n)
    args  = [{"phase":"montecarlo","seed":int(seeds[i]),"T":T,
              "overrides":{"WIND_TURB":float(turbs[i]),
                           "WIND_GUST_V":float(gusts[i])},
              "extra":{"turb":float(turbs[i]),"gust":float(gusts[i])}}
             for i in range(n)]
    res   = _run_parallel(args, workers)
    dPs   = [r["metrics"]["dP_mean"] for r in res if r["status"]=="ok"]
    pct   = 100 * sum(d>0 for d in dPs) / max(len(dPs),1)
    print(f"   ΔP = {np.mean(dPs):+.3f} ± {np.std(dPs):.3f} MW  "
          f"ΔP>0: {sum(d>0 for d in dPs)}/{len(dPs)} ({pct:.0f}%)")
    return res


def phase_ablation(T=T_SHORT, workers=2):
    print(f"\n── 4. ABLATION STUDY ──")
    configs = {
        "full":        {},
        "no_fft":      {"disable_fft":True},
        "no_tsr":      {"disable_tsr":True},
        "no_kind":     {"disable_kind":True},
        "no_load_pen": {"disable_load":True},
    }
    abl = {}
    for name, extra in configs.items():
        args = [{"phase":f"ablation_{name}","seed":s,"T":T,
                 "overrides":{},"extra":extra}
                for s in range(5)]
        res  = _run_parallel(args, workers)
        dPs  = [r["metrics"]["dP_mean"] for r in res if r["status"]=="ok"]
        taus = [r["metrics"]["tau_var"]  for r in res if r["status"]=="ok"]
        abl[name] = {"dP_mean":float(np.mean(dPs)),
                     "dP_std": float(np.std(dPs)),
                     "tau_mean":float(np.mean(taus))}
        print(f"   {name:<18}  ΔP={np.mean(dPs):+.3f}±{np.std(dPs):.3f} MW  "
              f"τ_var={np.mean(taus):.2e}")

    full_dP = abl["full"]["dP_mean"]
    print("   Contribucions marginals:")
    for k,v in abl.items():
        if k != "full":
            print(f"     -{k[3:]:12s}: {v['dP_mean']-full_dP:+.3f} MW")
    return abl


def phase_surrogate(base_runs, T=120.0, workers=2):
    print(f"\n── 5. SURROGATE + UCB ──")
    # Punt de partida: policy_x de les runs baseline
    X, J = [], []
    for r in base_runs:
        if r.get("status")=="ok":
            px = r["metrics"]["policy_x"]
            if len(px)==5:
                X.append(px)
                J.append(-r["metrics"]["dP_mean"])

    # Punts extra aleatoris
    rng = np.random.default_rng(77)
    extra_args = []
    extra_xs   = []
    for _ in range(20):
        x = rng.uniform(0.1, 0.9, 5)
        extra_xs.append(x)
        pw, me, kb, nm, ki = mod.decode_policy(x)
        extra_args.append({"phase":"surrogate_extra","seed":999,"T":T,
                           "overrides":{"K_Q3_BUS_0":float(kb),
                                        "K_Q3_IND_0":float(ki)}})
    extra_res = _run_parallel(extra_args, workers)
    for i, r in enumerate(extra_res):
        if r.get("status")=="ok":
            X.append(list(extra_xs[i]))
            J.append(-r["metrics"]["dP_mean"])

    if len(X) < 3:
        print("   Dades insuficients per al surrogate.")
        return {}

    Xa = np.array(X); Ja = np.array(J)
    print(f"   Punts recollits: {len(Xa)}  "
          f"ΔP rang: [{-Ja.max():.2f}, {-Ja.min():.2f}] MW")

    # Kernel Ridge Regression simple
    sigma = 0.3; lam = 1e-3
    def K(A, B):
        d = A[:,None,:] - B[None,:,:]
        return np.exp(-np.sum(d**2, axis=-1) / (2*sigma**2))

    Kt = K(Xa, Xa) + lam*np.eye(len(Xa))
    try:
        ac = np.linalg.solve(Kt, Ja)
    except Exception:
        ac = np.linalg.lstsq(Kt, Ja, rcond=None)[0]

    # UCB sobre candidats aleatoris
    rng2   = np.random.default_rng(123)
    Xc     = rng2.uniform(0.05, 0.95, (500, 5))
    Jp     = K(Xc, Xa) @ ac
    K_cc   = np.ones(len(Xc))
    K_inv  = K(Xc, Xa) @ np.linalg.pinv(Kt) @ K(Xc, Xa).T
    var_p  = np.maximum(0, K_cc - np.diag(K_inv))
    ucb    = Jp - 0.1 * np.sqrt(var_p)
    bi     = int(np.argmin(ucb))

    pw, me, kb, nm, ki = mod.decode_policy(Xc[bi])
    print(f"   UCB millor: p={pw:.2f} m={me/1e3:.1f}t "
          f"K_bus={kb:.4f} N=3:{nm} K_ind={ki:.4f}")

    # Valida
    val = _run_single({"phase":"surrogate_valid","seed":7,"T":T_SHORT,
                       "overrides":{"K_Q3_BUS_0":float(kb),
                                    "K_Q3_IND_0":float(ki)}})
    dP_val = val["metrics"]["dP_mean"] if val.get("status")=="ok" else 0.0
    print(f"   Validació real: ΔP={dP_val:+.3f} MW")

    best = {"p_win":float(pw),"m_exc_t":float(me/1e3),
            "K_bus":float(kb),"N_mode":int(nm),"K_ind":float(ki),
            "dP_real":dP_val}

    # Desa best_policy.txt
    with open(BEST_POLICY_FILE, "w") as f:
        f.write("MILLORS PARÀMETRES — Surrogate UCB\n" + "="*38 + "\n")
        for k2, v2 in best.items():
            f.write(f"  {k2}: {v2}\n")

    return {"best": best, "n_points": len(Xa)}


# ═══════════════════════════════════════════════════════════════════════════════
# GRÀFICS DE RESUM
# ═══════════════════════════════════════════════════════════════════════════════

def plot_summary(all_runs, abl, surr):
    fig = plt.figure(figsize=(20, 14))
    fig.patch.set_facecolor('#0d1117')
    gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.52, wspace=0.40)

    def mkax(r,c): return fig.add_subplot(gs[r,c])
    def style(a, title, xl='', yl=''):
        a.set_facecolor(BG); a.set_title(title, color=TXT, fontsize=9, pad=4)
        if xl: a.set_xlabel(xl, color=TXT, fontsize=8)
        if yl: a.set_ylabel(yl, color=TXT, fontsize=8)
        a.tick_params(colors=TXT, labelsize=7)
        for sp in a.spines.values(): sp.set_color('#30363d')
        a.grid(True, color='#21262d', linewidth=0.4)

    # Separa per fase
    def get_phase(phase):
        return [r for r in all_runs
                if r.get("phase")==phase and r.get("status")=="ok"]

    base_ok = get_phase("baseline")
    pol_ok  = get_phase("policy")
    mc_ok   = get_phase("montecarlo")

    # 1. Baseline bars
    a1 = mkax(0,0)
    if base_ok:
        dPs = [r["metrics"]["dP_mean"] for r in base_ok]
        a1.bar(range(len(dPs)), dPs,
               color=[CG if d>0 else CR for d in dPs], width=0.7, alpha=0.8)
        a1.axhline(0, color='white', lw=0.8, ls='--', alpha=0.5)
        a1.axhline(np.mean(dPs), color=CQ, lw=1.5, ls='--',
                   label=f'mitj={np.mean(dPs):+.2f}')
        a1.legend(fontsize=7, facecolor=BG, edgecolor='#30363d', labelcolor=TXT)
    style(a1, '1. Baseline ΔP per seed', 'seed', 'ΔP [MW]')

    # 2. Política ΔP vs J
    a2 = mkax(0,1)
    if pol_ok:
        dPs_p = [r["metrics"]["dP_mean"] for r in pol_ok]
        Js_p  = [r["metrics"]["J_final"]  for r in pol_ok]
        a2.scatter(dPs_p, Js_p, color=CV, s=60, zorder=5)
        for i,(x,y) in enumerate(zip(dPs_p,Js_p)):
            a2.annotate(f's{i}',(x,y),fontsize=7,color=TXT,
                        xytext=(3,3),textcoords='offset points')
        a2.axvline(0, color='white', lw=0.5, ls='--', alpha=0.4)
    style(a2, '2. Política: ΔP vs J_final', 'ΔP [MW]', 'J')

    # 3. MC histograma
    a3 = mkax(0,2)
    if mc_ok:
        dPs_m = [r["metrics"]["dP_mean"] for r in mc_ok]
        bins  = np.linspace(min(dPs_m)-0.5, max(dPs_m)+0.5, 18)
        a3.hist([d for d in dPs_m if d>=0], bins=bins,
                color=CG, alpha=0.75, label=f'≥0 ({sum(d>=0 for d in dPs_m)})')
        a3.hist([d for d in dPs_m if d<0],  bins=bins,
                color=CR, alpha=0.75, label=f'<0 ({sum(d<0 for d in dPs_m)})')
        a3.axvline(np.mean(dPs_m), color=CQ, lw=1.5, ls='--',
                   label=f'mitj={np.mean(dPs_m):+.2f}')
        pct = 100*sum(d>=0 for d in dPs_m)/len(dPs_m)
        a3.text(0.97,0.95,f'{pct:.0f}% >0',transform=a3.transAxes,
                ha='right',va='top',color=CG if pct>50 else CR,fontsize=9)
        a3.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a3, '3. Monte Carlo ΔP (hist.)', 'ΔP [MW]', 'N runs')

    # 4. MC ΔP ordenat
    a4 = mkax(0,3)
    if mc_ok:
        dPs_sorted = sorted([r["metrics"]["dP_mean"] for r in mc_ok])
        a4.scatter(range(len(dPs_sorted)), dPs_sorted,
                   c=[CG if d>0 else CR for d in dPs_sorted], s=20, alpha=0.8)
        a4.axhline(0, color='white', lw=0.8, ls='--', alpha=0.5)
        a4.axhline(np.mean(dPs_sorted), color=CQ, lw=1.2, ls='--')
    style(a4, 'MC ΔP ordenat', 'rang', 'ΔP [MW]')

    # 5. Ablation barres
    a5 = mkax(1, slice(0,2))
    if abl:
        names  = list(abl.keys())
        dP_v   = [abl[n]["dP_mean"] for n in names]
        dP_e   = [abl[n]["dP_std"]  for n in names]
        a5.bar(range(len(names)), dP_v,
               color=[CG if d>0 else CR for d in dP_v],
               width=0.6, alpha=0.8,
               yerr=dP_e, capsize=4,
               error_kw={"ecolor":TXT,"elinewidth":0.8})
        a5.set_xticks(range(len(names)))
        a5.set_xticklabels([n.replace('_',' ') for n in names],
                           color=TXT, fontsize=7, rotation=12)
        a5.axhline(0, color='white', lw=0.8, ls='--', alpha=0.5)
        full_dP = abl.get("full",{"dP_mean":0})["dP_mean"]
        for i,v in enumerate(dP_v):
            a5.text(i, v+0.1 if v>0 else v-0.2,
                    f'{v-full_dP:+.2f}' if i>0 else 'ref',
                    ha='center', color=TXT, fontsize=7)
    style(a5, '4. Ablation — ΔP per config.', '', 'ΔP [MW]')

    # 6. τ_var per config
    a6 = mkax(1,2)
    if abl:
        names  = list(abl.keys())
        tau_v  = [abl[n]["tau_mean"] for n in names]
        tau_lim= getattr(mod,"TAU_VAR_MAX",1.5e12)
        a6.bar(range(len(names)), tau_v,
               color=[CR if t>tau_lim else CM for t in tau_v],
               width=0.6, alpha=0.8)
        a6.set_xticks(range(len(names)))
        a6.set_xticklabels([n.replace('_',' ') for n in names],
                           color=TXT, fontsize=7, rotation=12)
        a6.axhline(tau_lim, color='white', lw=0.8, ls='--',
                   alpha=0.5, label='límit struct.')
        a6.legend(fontsize=7,facecolor=BG,edgecolor='#30363d',labelcolor=TXT)
    style(a6, 'Ablation: var(τ)', '', 'N²m²')

    # 7. Surrogate + best policy
    a7 = mkax(1,3)
    a7.set_facecolor(BG); a7.axis('off')
    bp = surr.get("best", {}) if surr else {}
    lines = [
        ("SURROGATE UCB", CQ),
        (f"p_win  = {bp.get('p_win',0):.3f}", TXT),
        (f"m_exc  = {bp.get('m_exc_t',0):.1f} t", TXT),
        (f"K_bus  = {bp.get('K_bus',0):.5f}", TXT),
        (f"N_mode = 3:{bp.get('N_mode',12)}", TXT),
        (f"K_ind  = {bp.get('K_ind',0):.5f}", TXT),
        ("", TXT),
        (f"ΔP_val = {bp.get('dP_real',0):+.3f} MW",
         CG if bp.get('dP_real',0)>0 else CR),
        (f"Punts training: {surr.get('n_points',0)}", TXT),
    ]
    for li,(ln,col) in enumerate(lines):
        a7.text(0.06, 0.93-li*0.10, ln, color=col, fontsize=9,
                fontfamily='monospace', transform=a7.transAxes)

    # 8. Vies d'energia (mitja de runs base)
    a8 = mkax(2,slice(0,2))
    if base_ok:
        v1s = [r["metrics"]["via1_MJ"] for r in base_ok]
        v2s = [r["metrics"]["via2_MJ"] for r in base_ok]
        v3s = [r["metrics"]["via3_MJ"] for r in base_ok]
        shs = [r["metrics"]["shield_MJ"] for r in base_ok]
        labels = ['Via 1\nshear','Via 2\ngrav.','Via 3\nres.','Escut']
        means  = [np.mean(v1s),np.mean(v2s),np.mean(v3s),np.mean(shs)]
        stds   = [np.std(v1s), np.std(v2s), np.std(v3s), np.std(shs)]
        cols8  = [CG if m>0 else CR for m in means]
        a8.bar(range(4), means, color=cols8, width=0.6, alpha=0.8,
               yerr=stds, capsize=4,
               error_kw={"ecolor":TXT,"elinewidth":0.8})
        a8.set_xticks(range(4))
        a8.set_xticklabels(labels, color=TXT, fontsize=8)
        a8.axhline(0,color='white',lw=0.8,ls='--',alpha=0.5)
    style(a8, 'Desglos energètic vies (mitja baseline)', '', 'MJ')

    # 9. Resum executiu text
    a9 = mkax(2,slice(2,4))
    a9.set_facecolor(BG); a9.axis('off')
    dPs_b = [r["metrics"]["dP_mean"] for r in base_ok] if base_ok else [0]
    dPs_m = [r["metrics"]["dP_mean"] for r in mc_ok]   if mc_ok  else [0]
    pos_b = sum(d>0 for d in dPs_b); pos_m = sum(d>0 for d in dPs_m)
    sumlines = [
        ("RESUM EXPERIMENT v4", CQ),
        (f"1. Baseline  ΔP={np.mean(dPs_b):+.3f}±{np.std(dPs_b):.3f} MW  "
         f"ΔP>0:{pos_b}/{len(dPs_b)}", TXT),
        (f"3. MC (runs) ΔP={np.mean(dPs_m):+.3f}±{np.std(dPs_m):.3f} MW  "
         f"ΔP>0:{pos_m}/{len(dPs_m)} ({100*pos_m/max(1,len(dPs_m)):.0f}%)", TXT),
    ]
    if abl:
        full_dP = abl.get("full",{"dP_mean":0})["dP_mean"]
        sumlines.append(("4. Ablation (Δ vs full):", CQ))
        for k,v in abl.items():
            if k!="full":
                sumlines.append((f"   -{k[3:]:12s}: {v['dP_mean']-full_dP:+.3f} MW", TXT))
    if bp:
        sumlines.append(("5. Surrogate:", CQ))
        sumlines.append((f"   ΔP_real = {bp.get('dP_real',0):+.3f} MW",
                         CG if bp.get('dP_real',0)>0 else CR))

    for li,(ln,col) in enumerate(sumlines):
        a9.text(0.02, 0.97-li*0.082, ln, color=col, fontsize=8.5,
                fontfamily='monospace', transform=a9.transAxes)

    fig.suptitle(
        'Experiment v4 — Baseline · Política · MC · Ablation · Surrogate  |  '
        'Víctor Manzanares Alberola — EPSA/UPV Alcoi',
        color=TXT, fontsize=9, y=0.999)

    plt.savefig(SUMMARY_PNG, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"\n   Gràfic guardat: {SUMMARY_PNG}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Experiment runner v4")
    parser.add_argument("--quick",   action="store_true",
                        help="Versió ràpida: T=120s, 4 seeds, 10 MC")
    parser.add_argument("--workers", type=int, default=N_WORKERS)
    parser.add_argument("--t-short", type=float, default=T_SHORT)
    parser.add_argument("--no-mc",   action="store_true")
    parser.add_argument("--no-surr", action="store_true")
    args = parser.parse_args()

    if args.quick:
        T     = 120.0;  n_base = 4;  n_mc = 10
    else:
        T     = args.t_short;  n_base = N_SEEDS_BASE;  n_mc = N_SEEDS_MC

    w = args.workers

    print("═"*70)
    print("  experiment_v4_logging.py")
    print(f"  T={T:.0f}s  workers={w}  seeds={n_base}  MC={n_mc}")
    print(f"  OUT_DIR: {OUT_DIR}")
    print("═"*70)

    cp = load_checkpoint()
    if cp:
        print(f"  Checkpoint trobat: {cp.get('n_runs',0)} runs previes.")

    t0 = time.time()

    base = phase_baseline(n=n_base, T=T, workers=w)
    pol  = phase_policy(n=min(5,n_base), T=T, workers=1)

    mc   = []
    if not args.no_mc:
        mc = phase_montecarlo(n=n_mc, T=T, workers=w)

    abl  = phase_ablation(T=T, workers=max(1,w//2))

    surr = {}
    if not args.no_surr:
        surr = phase_surrogate(base, T=min(T,120.0), workers=max(1,w//2))

    # Consolida tots els runs
    all_runs = consolidate()

    print(f"\n  Generant gràfics...")
    plot_summary(all_runs, abl, surr)

    elapsed = (time.time()-t0)/60
    print(f"\n═══════════════════════════════════")
    print(f"  ✓ Completat en {elapsed:.1f} min")
    print(f"  {RESULTS_FILE}")
    print(f"  {SUMMARY_PNG}")
    if surr.get("best"):
        print(f"  {BEST_POLICY_FILE}")
    print("═══════════════════════════════════")

if __name__ == "__main__":
    main()
