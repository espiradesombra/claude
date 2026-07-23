"""
gemelo_xfi_avion_4turbinas.py
=============================
XFI — Experimental Flight Infinite (gemelo conceptual VMA)

Idea original (Víctor Manzanares):
  - Avió: 3 motors ZypyZape (gen / thr / buf)
  - Variant: 4 turbines (cotxe/camió o craft)

Ús:
  python gemelo_xfi_avion_4turbinas.py          # N=4
  python gemelo_xfi_avion_4turbinas.py --n 3    # N=3 (avió clàssic del xat)
  python gemelo_xfi_avion_4turbinas.py --compare

Sortides:
  xfi_avion_dinamica_N3.png / N4.png
  mètriques a consola
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


G = 9.81
RHO = 1.225

M_CRAFT = 25.0
S_WING = 1.2
CD0 = 0.08
CL_MAX = 1.2
CL_CRUISE = 0.6

R_TURB = 0.35
A_TURB = math.pi * R_TURB ** 2
CP = 0.40
J_ROTOR = 0.8
OMEGA_MAX = 400.0
ETA_ELEC = 0.88
ETA_PROP = 0.75

F_ZZ = 0.35
P_ZZ_FRAC = 0.12

DT = 0.02
T_TOTAL = 120.0


def mission_phase(t: float) -> str:
    period = 30.0
    x = (t % period) / period
    if x < 0.30:
        return "climb"
    if x < 0.45:
        return "cruise"
    if x < 0.80:
        return "dive"
    return "pullout"


def wind_rel(v_air: float, v_craft: float, phase: str) -> float:
    if phase == "dive":
        return abs(v_air) + 0.35 * abs(v_craft)
    if phase == "climb":
        return max(0.0, abs(v_air) - 0.15 * abs(v_craft))
    return abs(v_air) + 0.05 * abs(v_craft)


@dataclass
class XFIState:
    n_motors: int = 3
    h: float = 200.0
    v: float = 22.0
    gamma: float = 0.0
    omega: np.ndarray = field(default_factory=lambda: np.ones(3) * 180.0)
    e_buffer: float = 0.0
    e_chem: float = 2.0e5
    t: float = 0.0
    hist_t: List[float] = field(default_factory=list)
    hist_h: List[float] = field(default_factory=list)
    hist_v: List[float] = field(default_factory=list)
    hist_ebuf: List[float] = field(default_factory=list)
    hist_pgen: List[float] = field(default_factory=list)
    hist_pthr: List[float] = field(default_factory=list)
    hist_phase: List[str] = field(default_factory=list)
    hist_roles: List[List[str]] = field(default_factory=list)
    hist_n_gen: List[int] = field(default_factory=list)
    hist_n_thr: List[int] = field(default_factory=list)


def rotor_ke(omega: np.ndarray) -> float:
    return 0.5 * J_ROTOR * float(np.sum(omega ** 2))


def motor_roles(t: float, n: int) -> List[str]:
    """
    N=3 (avió XFI del xat): un gen, un thr, un buf a 120° — sempre.
    N=4: roda gen/thr/buf amb fractions de 2π/3 (pot haver-hi 2 thr).
    """
    if n == 3:
        # Desfasatge 120° fix: index per temps
        # slot 0→gen, 1→thr, 2→buf, rotació cada 1/(3*f_zz)
        base = int(math.floor(3 * F_ZZ * t)) % 3
        order = ["gen", "thr", "buf"]
        return [order[(base + i) % 3] for i in range(3)]

    roles = []
    phase = 2 * math.pi * F_ZZ * t
    for i in range(n):
        ang = (phase + 2 * math.pi * i / n) % (2 * math.pi)
        if ang < 2 * math.pi / 3:
            roles.append("gen")
        elif ang < 4 * math.pi / 3:
            roles.append("thr")
        else:
            roles.append("buf")
    if "thr" not in roles and n >= 1:
        roles[0] = "thr"
    if "gen" not in roles and n >= 2:
        roles[1] = "gen"
    return roles


def step(state: XFIState, dt: float = DT) -> None:
    n = state.n_motors
    phase = mission_phase(state.t)
    roles = motor_roles(state.t, n)

    if phase == "climb":
        gamma_cmd = math.radians(12.0)
    elif phase == "dive":
        gamma_cmd = math.radians(-18.0)
    elif phase == "pullout":
        gamma_cmd = math.radians(5.0)
    else:
        gamma_cmd = math.radians(1.0)

    state.gamma += 0.4 * (gamma_cmd - state.gamma) * dt

    v = max(state.v, 3.0)
    q = 0.5 * RHO * v * v
    cl = CL_CRUISE * (1.0 + 0.5 * math.sin(state.gamma))
    cl = max(-CL_MAX, min(CL_MAX, cl))
    cd = CD0 + 0.05 * cl * cl
    D = q * S_WING * cd

    v_rel = wind_rel(v, v, phase)
    p_disc = 0.5 * RHO * A_TURB * (v_rel ** 3) * CP

    p_gen_total = 0.0
    p_thr_total = 0.0
    omega_new = state.omega.copy()

    for i, role in enumerate(roles):
        if role == "gen":
            p_in = p_disc * ETA_ELEC
            t_load = p_in / max(omega_new[i], 1.0)
            alpha = (t_load / J_ROTOR) - 0.02 * omega_new[i]
            omega_new[i] = min(OMEGA_MAX, max(5.0, omega_new[i] + alpha * dt))
            p_gen_total += p_in
            state.e_chem += p_in * dt * 0.15
        elif role == "thr":
            p_want = 400.0 + 200.0 * max(0.0, math.sin(state.gamma + 0.2))
            p_avail = min(p_want, state.e_chem / max(dt, 1e-3) * 0.5 + 50.0)
            e_rot = 0.5 * J_ROTOR * omega_new[i] ** 2
            take = min(e_rot * 0.05, p_avail * dt)
            if e_rot > 1.0:
                omega_new[i] = math.sqrt(max(25.0, 2.0 * (e_rot - take) / J_ROTOR))
            state.e_chem = max(0.0, state.e_chem - p_avail * dt)
            p_thr_total += p_avail * ETA_PROP
        else:
            j = (i + 1) % n
            mean = 0.5 * (omega_new[i] + omega_new[j])
            k = P_ZZ_FRAC
            omega_new[i] = (1 - k) * omega_new[i] + k * mean
            omega_new[j] = (1 - k) * omega_new[j] + k * mean
            omega_new[i] *= 0.999
            omega_new[j] *= 0.999

    state.omega = omega_new
    state.e_buffer = rotor_ke(state.omega)

    thrust = p_thr_total / v
    weight_along = -M_CRAFT * G * math.sin(state.gamma)
    f_net = thrust - D + weight_along
    a = f_net / M_CRAFT
    state.v = max(5.0, state.v + a * dt)

    vh = state.v * math.sin(state.gamma)
    state.h = max(5.0, state.h + vh * dt)
    state.e_chem = max(0.0, state.e_chem - 15.0 * dt)

    state.t += dt
    state.hist_t.append(state.t)
    state.hist_h.append(state.h)
    state.hist_v.append(state.v)
    state.hist_ebuf.append(state.e_buffer)
    state.hist_pgen.append(p_gen_total)
    state.hist_pthr.append(p_thr_total)
    state.hist_phase.append(phase)
    state.hist_roles.append(roles)
    state.hist_n_gen.append(roles.count("gen"))
    state.hist_n_thr.append(roles.count("thr"))


def simulate(n_motors: int = 3, t_total: float = T_TOTAL) -> XFIState:
    st = XFIState(
        n_motors=n_motors,
        omega=np.ones(n_motors) * 180.0,
    )
    steps = int(t_total / DT)
    for _ in range(steps):
        step(st)
    return st


def metrics(st: XFIState) -> dict:
    h = np.array(st.hist_h)
    v = np.array(st.hist_v)
    ebuf = np.array(st.hist_ebuf)
    pgen = np.array(st.hist_pgen)
    pthr = np.array(st.hist_pthr)
    # alçada neta: mitjana final 20s vs inicial 5s
    h0 = float(np.mean(h[: max(1, int(5 / DT))]))
    hf = float(np.mean(h[-max(1, int(20 / DT)) :]))
    return {
        "n": st.n_motors,
        "h_min": float(h.min()),
        "h_max": float(h.max()),
        "h0": h0,
        "hf": hf,
        "dh": hf - h0,
        "v_min": float(v.min()),
        "v_max": float(v.max()),
        "v_mean": float(v.mean()),
        "ebuf_min": float(ebuf.min()),
        "ebuf_max": float(ebuf.max()),
        "ebuf_mean": float(ebuf.mean()),
        "pgen_mean": float(pgen.mean()),
        "pthr_mean": float(pthr.mean()),
        "e_chem_final": float(st.e_chem),
        "n_gen_mean": float(np.mean(st.hist_n_gen)),
        "n_thr_mean": float(np.mean(st.hist_n_thr)),
        "ratio_gen_thr": float(pgen.mean() / max(pthr.mean(), 1e-9)),
    }


def plot_report(st: XFIState, out_dir: str, tag: Optional[str] = None) -> str:
    n = st.n_motors
    tag = tag or f"N{n}"
    t = np.array(st.hist_t)
    h = np.array(st.hist_h)
    v = np.array(st.hist_v)
    ebuf = np.array(st.hist_ebuf)
    pgen = np.array(st.hist_pgen)
    pthr = np.array(st.hist_pthr)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), facecolor="#0f1419")
    for ax in axes.ravel():
        ax.set_facecolor("#1a222c")
        ax.tick_params(colors="#c8d0d8")
        for spine in ax.spines.values():
            spine.set_color("#3a4550")
        ax.xaxis.label.set_color("#a8b0b8")
        ax.yaxis.label.set_color("#a8b0b8")
        ax.title.set_color("#e8eef4")
        ax.grid(True, alpha=0.2, color="#5a6570")

    axes[0, 0].plot(t, h, color="#5ec8ff", lw=1.6)
    axes[0, 0].set_title("Alçada h(t)")
    axes[0, 0].set_xlabel("t [s]")
    axes[0, 0].set_ylabel("h [m]")

    axes[0, 1].plot(t, v, color="#7dffa0", lw=1.6)
    axes[0, 1].set_title("Velocitat aire v(t)")
    axes[0, 1].set_xlabel("t [s]")
    axes[0, 1].set_ylabel("v [m/s]")

    axes[1, 0].plot(t, ebuf / 1000.0, color="#ffb86c", lw=1.6)
    axes[1, 0].set_title("Buffer cinètic rotores")
    axes[1, 0].set_xlabel("t [s]")
    axes[1, 0].set_ylabel("E [kJ]")

    axes[1, 1].plot(t, pgen / 1000.0, color="#c792ea", lw=1.2, label="P_gen")
    axes[1, 1].plot(t, pthr / 1000.0, color="#ff6b6b", lw=1.2, label="P_thrust")
    axes[1, 1].set_title(f"Potències — N={n}, f_zz={F_ZZ} Hz")
    axes[1, 1].set_xlabel("t [s]")
    axes[1, 1].set_ylabel("P [kW]")
    axes[1, 1].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8")

    fig.suptitle(
        f"XFI — Avió/drone ZypyZape  ·  N={n} motors\n"
        f"Rols: gen / thr / buf  ·  NO perpetuum (η<1)",
        color="#f0f4f8",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(out_dir, f"xfi_avion_dinamica_{tag}.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_compare(st3: XFIState, st4: XFIState, out_dir: str) -> str:
    t3 = np.array(st3.hist_t)
    t4 = np.array(st4.hist_t)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), facecolor="#0f1419")
    for ax in axes.ravel():
        ax.set_facecolor("#1a222c")
        ax.tick_params(colors="#c8d0d8")
        for spine in ax.spines.values():
            spine.set_color("#3a4550")
        ax.xaxis.label.set_color("#a8b0b8")
        ax.yaxis.label.set_color("#a8b0b8")
        ax.title.set_color("#e8eef4")
        ax.grid(True, alpha=0.2, color="#5a6570")

    axes[0, 0].plot(t3, st3.hist_h, color="#5ec8ff", lw=1.5, label="N=3")
    axes[0, 0].plot(t4, st4.hist_h, color="#ff9f43", lw=1.5, label="N=4", alpha=0.85)
    axes[0, 0].set_title("Alçada")
    axes[0, 0].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8")

    axes[0, 1].plot(t3, st3.hist_v, color="#5ec8ff", lw=1.5, label="N=3")
    axes[0, 1].plot(t4, st4.hist_v, color="#ff9f43", lw=1.5, label="N=4", alpha=0.85)
    axes[0, 1].set_title("Velocitat")
    axes[0, 1].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8")

    axes[1, 0].plot(t3, np.array(st3.hist_ebuf) / 1000, color="#5ec8ff", lw=1.5, label="N=3")
    axes[1, 0].plot(t4, np.array(st4.hist_ebuf) / 1000, color="#ff9f43", lw=1.5, label="N=4", alpha=0.85)
    axes[1, 0].set_title("Buffer rotors [kJ]")
    axes[1, 0].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8")

    axes[1, 1].plot(t3, np.array(st3.hist_pgen) / 1000, color="#5ec8ff", lw=1.2, label="Pgen N=3")
    axes[1, 1].plot(t3, np.array(st3.hist_pthr) / 1000, color="#5ec8ff", lw=1.2, ls="--", label="Pthr N=3")
    axes[1, 1].plot(t4, np.array(st4.hist_pgen) / 1000, color="#ff9f43", lw=1.2, label="Pgen N=4")
    axes[1, 1].plot(t4, np.array(st4.hist_pthr) / 1000, color="#ff9f43", lw=1.2, ls="--", label="Pthr N=4")
    axes[1, 1].set_title("Potències [kW]")
    axes[1, 1].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=8)

    fig.suptitle("XFI comparativa N=3 (avió) vs N=4", color="#f0f4f8", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(out_dir, "xfi_avion_comparativa_N3_vs_N4.png")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def print_metrics(m: dict) -> None:
    print(f"  N motors     = {m['n']}")
    print(f"  h  min/max   = {m['h_min']:.1f} / {m['h_max']:.1f} m")
    print(f"  h0 → hf (Δh) = {m['h0']:.1f} → {m['hf']:.1f}  (Δ={m['dh']:+.1f} m)")
    print(f"  v  min/max/μ = {m['v_min']:.1f} / {m['v_max']:.1f} / {m['v_mean']:.1f} m/s")
    print(f"  Ebuf min/max = {m['ebuf_min']:.0f} / {m['ebuf_max']:.0f} J  (μ={m['ebuf_mean']:.0f})")
    print(f"  Pgen / Pthr  = {m['pgen_mean']:.1f} / {m['pthr_mean']:.1f} W  (ratio={m['ratio_gen_thr']:.2f})")
    print(f"  #gen / #thr  = {m['n_gen_mean']:.2f} / {m['n_thr_mean']:.2f}  (mitjana temporal)")
    print(f"  e_chem final = {m['e_chem_final']:.0f} J")


def main() -> None:
    ap = argparse.ArgumentParser(description="XFI gemelo avió ZypyZape")
    ap.add_argument("--n", type=int, default=3, help="nombre de motors (3=avió xat, 4=variant)")
    ap.add_argument("--compare", action="store_true", help="simula N=3 i N=4 i compara")
    ap.add_argument("--t", type=float, default=T_TOTAL, help="temps de sim [s]")
    args = ap.parse_args()

    out_dir = os.path.dirname(os.path.abspath(__file__))

    if args.compare:
        print("=" * 60)
        print("XFI comparativa N=3 vs N=4")
        print("=" * 60)
        st3 = simulate(3, args.t)
        st4 = simulate(4, args.t)
        m3, m4 = metrics(st3), metrics(st4)
        p3 = plot_report(st3, out_dir, "N3")
        p4 = plot_report(st4, out_dir, "N4")
        pc = plot_compare(st3, st4, out_dir)
        print("\n--- N=3 (avió XFI) ---")
        print_metrics(m3)
        print("\n--- N=4 ---")
        print_metrics(m4)
        print("\n--- Delta (N3 - N4) ---")
        print(f"  Δh_range   = {(m3['h_max']-m3['h_min']) - (m4['h_max']-m4['h_min']):+.1f} m")
        print(f"  Δv_mean    = {m3['v_mean']-m4['v_mean']:+.2f} m/s")
        print(f"  Δebuf_mean = {m3['ebuf_mean']-m4['ebuf_mean']:+.0f} J")
        print(f"  Δe_chem    = {m3['e_chem_final']-m4['e_chem_final']:+.0f} J")
        print(f"  Δratio P   = {m3['ratio_gen_thr']-m4['ratio_gen_thr']:+.2f}")
        print(f"\nfigures:\n  {p3}\n  {p4}\n  {pc}")
        return

    n = max(2, min(8, args.n))
    print("=" * 60)
    print(f"XFI gemelo — N={n} motors ZypyZape")
    print("=" * 60)
    st = simulate(n, args.t)
    m = metrics(st)
    path = plot_report(st, out_dir, f"N{n}")
    print_metrics(m)
    print(f"  figura      = {path}")
    if n == 3:
        print("\n  [N=3] Arquitectura del xat: 1 gen + 1 thr + 1 buf a 120°.")
        print("  Baixada capta, picat accelera, pujada gasta buffer/química.")


if __name__ == "__main__":
    main()
