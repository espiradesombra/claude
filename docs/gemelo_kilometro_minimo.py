"""
gemelo_kilometro_minimo.py
==========================
Kilómetro — gemelo mínimo con balance energético CERRADO (1ª ley).

Modelo:
  Masa m en círculo de radio R (ángulo θ).
  I = m R²
  τ_g    = -m g R cos(θ)     →  potencial V = m g R sin(θ)
  τ_fric = -c(θ,ω) · ω       →  disipa; fracción η_reg → E_útil
  τ_ext  = actuador          →  inyecta o frena

Identidad mecánica (exacta en el discreto si usamos el mismo τ_net):
  ΔE_mec = W_g + W_fric + W_ext
  E_fric_diss = -W_fric
  E_util + E_heat_fric = E_fric_diss

En régimen que DRENA el volante (baja E_cin):
  E_util puede ser > E_ext  sin magia: se gasta la batería cinética inicial.
  η_ciclo = E_util / (E_ext + max(0, -ΔE_mec_sin_pot?)) se reporta con cuidado.

Uso:
  python gemelo_kilometro_minimo.py
  python gemelo_kilometro_minimo.py --T 40 --eta 0.7 --mode maintain
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = 9.81


@dataclass
class KMParams:
    m: float = 80.0
    R: float = 2.5
    c_base: float = 2.5
    c_harvest: float = 16.0
    eta_reg: float = 0.65
    tau_ext_max: float = 100.0
    omega_target: float = 1.8   # rad/s  (modo maintain)
    omega0: float = 3.5         # arranque con "carga" cinética
    mode: str = "drain"         # drain | maintain
    dt: float = 0.001
    T: float = 30.0


@dataclass
class KMState:
    t: float = 0.0
    theta: float = 0.0
    omega: float = 3.5
    # trabajos acumulados [J]  W = ∫ τ ω dt
    W_g: float = 0.0
    W_fric: float = 0.0
    W_ext: float = 0.0
    E_util: float = 0.0
    E_heat: float = 0.0
    hist_t: List[float] = field(default_factory=list)
    hist_omega: List[float] = field(default_factory=list)
    hist_Ekin: List[float] = field(default_factory=list)
    hist_Epot: List[float] = field(default_factory=list)
    hist_Emec: List[float] = field(default_factory=list)
    hist_P_util: List[float] = field(default_factory=list)
    hist_P_ext: List[float] = field(default_factory=list)


def inertia(p: KMParams) -> float:
    return p.m * p.R * p.R


def pot(theta: float, p: KMParams) -> float:
    return p.m * G * p.R * math.sin(theta)


def torques(theta: float, omega: float, p: KMParams) -> Tuple[float, float, float]:
    tau_g = -p.m * G * p.R * math.cos(theta)

    # Cosecha: más fricción cuando τ_g acelera a favor de ω
    assisting = (tau_g * omega) > 0.0 and abs(omega) > 0.4
    c = p.c_harvest if assisting else p.c_base
    tau_fric = -c * omega

    tau_ext = 0.0
    if p.mode == "maintain":
        # control P sobre |ω| hacia omega_target
        err = p.omega_target - abs(omega)
        direction = 1.0 if omega >= 0 else -1.0
        if abs(omega) < 1e-4:
            direction = 1.0
        if err > 0:
            tau_ext = direction * min(p.tau_ext_max, 25.0 + 40.0 * err)
        elif err < -0.3:
            # freno suave del actuador
            tau_ext = -direction * min(p.tau_ext_max * 0.4, 15.0 * (-err))
    else:
        # drain: solo rescata si se para casi del todo
        if abs(omega) < 0.35:
            direction = 1.0 if omega >= 0 else -1.0
            if abs(omega) < 1e-4:
                direction = 1.0
            tau_ext = direction * 35.0

    return tau_g, tau_fric, tau_ext


def step(st: KMState, p: KMParams) -> None:
    I = inertia(p)
    tau_g, tau_fric, tau_ext = torques(st.theta, st.omega, p)
    tau_net = tau_g + tau_fric + tau_ext

    # semi-implícito simple: α con ω actual
    alpha = tau_net / I
    omega_mid = st.omega + 0.5 * alpha * p.dt  # para potencias más estables
    # recompute torques lightly with omega_mid for work (opcional: keep simple)
    P_g = tau_g * st.omega
    P_fric = tau_fric * st.omega
    P_ext = tau_ext * st.omega

    st.W_g += P_g * p.dt
    st.W_fric += P_fric * p.dt
    st.W_ext += P_ext * p.dt

    # regeneración solo si fricción extrae (P_fric < 0)
    if P_fric < 0:
        harvested = -P_fric * p.dt
        st.E_util += p.eta_reg * harvested
        st.E_heat += (1.0 - p.eta_reg) * harvested
    # si actuador frena (P_ext < 0), eso también es calor (no regeneramos actuador aquí)
    if P_ext < 0:
        st.E_heat += -P_ext * p.dt

    st.omega += alpha * p.dt
    st.theta += st.omega * p.dt
    # no wrap agresivo para pot continua; unwrap libre
    st.t += p.dt

    Ekin = 0.5 * I * st.omega ** 2
    Epot = pot(st.theta, p)
    st.hist_t.append(st.t)
    st.hist_omega.append(st.omega)
    st.hist_Ekin.append(Ekin)
    st.hist_Epot.append(Epot)
    st.hist_Emec.append(Ekin + Epot)
    st.hist_P_util.append(max(0.0, -P_fric) * p.eta_reg)
    st.hist_P_ext.append(max(0.0, P_ext))


def simulate(p: KMParams) -> KMState:
    st = KMState(omega=p.omega0)
    n = int(p.T / p.dt)
    for _ in range(n):
        step(st, p)
        if abs(st.omega) > 80:  # safety
            break
    return st


def analyze(st: KMState, p: KMParams) -> dict:
    Emec0 = st.hist_Emec[0]
    Emec1 = st.hist_Emec[-1]
    dEmec = Emec1 - Emec0
    Ekin0, Ekin1 = st.hist_Ekin[0], st.hist_Ekin[-1]
    Epot0, Epot1 = st.hist_Epot[0], st.hist_Epot[-1]
    dEkin = Ekin1 - Ekin0
    dEpot = Epot1 - Epot0

    # Identidades (mecánica clásica, τ_g = -∂V/∂θ):
    #   ΔE_kin = W_g + W_fric + W_ext
    #   ΔE_pot = -W_g
    #   ΔE_mec = W_fric + W_ext          ← residual debe ~0
    residual_mec = dEmec - (st.W_fric + st.W_ext)
    residual_kin = dEkin - (st.W_g + st.W_fric + st.W_ext)
    residual_pot = dEpot - (-st.W_g)

    E_fric_diss = max(0.0, -st.W_fric)
    if len(st.hist_t) > 2:
        trap = getattr(np, "trapezoid", None) or np.trapz
        E_ext_pos = float(trap(st.hist_P_ext, st.hist_t))
    else:
        E_ext_pos = 0.0

    # E_util sale de disipar fricción regenerada; el "pago" es:
    # actuador positivo + vaciado de energía mecánica (drain)
    drain = max(0.0, -dEmec)
    paid = E_ext_pos + drain
    eta_paid = st.E_util / paid if paid > 1e-9 else float("nan")

    return {
        "dEmec": dEmec,
        "dEkin": dEkin,
        "dEpot": dEpot,
        "W_g": st.W_g,
        "W_fric": st.W_fric,
        "W_ext": st.W_ext,
        "residual": residual_mec,
        "residual_kin": residual_kin,
        "residual_pot": residual_pot,
        "E_util": st.E_util,
        "E_heat": st.E_heat,
        "E_fric_diss": E_fric_diss,
        "E_ext_pos": E_ext_pos,
        "drain": drain,
        "paid": paid,
        "eta_paid": eta_paid,
        "eta_reg": p.eta_reg,
        "P_util_mean": st.E_util / max(st.t, 1e-9),
        "P_ext_mean": E_ext_pos / max(st.t, 1e-9),
        "omega_final": st.omega,
        "Emec0": Emec0,
        "Emec1": Emec1,
    }


def plot_report(st: KMState, a: dict, p: KMParams, path: str) -> None:
    t = np.array(st.hist_t)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), facecolor="#0f1419")
    for ax in axes.ravel():
        ax.set_facecolor("#1a222c")
        ax.tick_params(colors="#c8d0d8")
        for sp in ax.spines.values():
            sp.set_color("#3a4550")
        ax.xaxis.label.set_color("#a8b0b8")
        ax.yaxis.label.set_color("#a8b0b8")
        ax.title.set_color("#e8eef4")
        ax.grid(True, alpha=0.2, color="#5a6570")

    axes[0, 0].plot(t, st.hist_omega, color="#5ec8ff", lw=1.2)
    axes[0, 0].set_title(f"ω(t)  mode={p.mode}")
    axes[0, 0].set_xlabel("t [s]")
    axes[0, 0].set_ylabel("rad/s")

    axes[0, 1].plot(t, st.hist_Ekin, color="#ffb86c", lw=1.1, label="E_cin")
    axes[0, 1].plot(t, st.hist_Epot, color="#7dffa0", lw=1.0, alpha=0.8, label="E_pot")
    axes[0, 1].plot(t, st.hist_Emec, color="#c792ea", lw=1.2, label="E_mec")
    axes[0, 1].set_title("Energía mecánica (batería cinética)")
    axes[0, 1].set_xlabel("t [s]")
    axes[0, 1].set_ylabel("J")
    axes[0, 1].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=8)

    axes[1, 0].plot(t, st.hist_P_util, color="#c792ea", lw=1.0, label="P_útil")
    axes[1, 0].plot(t, st.hist_P_ext, color="#ff6b6b", lw=1.0, label="P_ext>0")
    axes[1, 0].set_title("Potencias")
    axes[1, 0].set_xlabel("t [s]")
    axes[1, 0].set_ylabel("W")
    axes[1, 0].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=8)

    labels = ["E_ext\n(actuador+)", "drain\n(-ΔE_mec)", "E_util", "E_heat"]
    vals = [a["E_ext_pos"], a["drain"], a["E_util"], a["E_heat"]]
    colors = ["#ff6b6b", "#5ec8ff", "#7dffa0", "#ffb86c"]
    bars = axes[1, 1].bar(labels, vals, color=colors, edgecolor="#3a4550")
    axes[1, 1].set_title(
        f"Balance  ·  res(ΔEmec-W_f-W_e)={a['residual']:.2f} J  ·  η_paid={a['eta_paid']:.2f}"
    )
    axes[1, 1].set_ylabel("J")
    for b, v in zip(bars, vals):
        axes[1, 1].text(b.get_x() + b.get_width() / 2, b.get_height(),
                        f"{v:.0f}", ha="center", va="bottom", fontsize=8, color="#c8d0d8")

    fig.suptitle(
        "Kilómetro mínimo — 1ª ley  ΔE_mec = W_fric + W_ext  (W_g cancela en pot)\n"
        f"m={p.m} kg  R={p.R} m  η_reg={p.eta_reg}  ·  "
        f"W_g={a['W_g']:.1f} J  residual={a['residual']:.2f} J",
        color="#f0f4f8", fontsize=11, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=30.0)
    ap.add_argument("--eta", type=float, default=0.65)
    ap.add_argument("--m", type=float, default=80.0)
    ap.add_argument("--R", type=float, default=2.5)
    ap.add_argument("--mode", choices=["drain", "maintain"], default="drain")
    ap.add_argument("--omega0", type=float, default=3.5)
    args = ap.parse_args()

    p = KMParams(T=args.T, eta_reg=args.eta, m=args.m, R=args.R,
                 mode=args.mode, omega0=args.omega0)
    st = simulate(p)
    a = analyze(st, p)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "kilometro_minimo_balance.png")
    plot_report(st, a, p, out)

    print("=" * 62)
    print("KILÓMETRO MÍNIMO — balance cerrado (1ª ley)")
    print("=" * 62)
    print(f"  mode           = {p.mode}")
    print(f"  T              = {st.t:.2f} s")
    print(f"  ω0 → ωf        = {p.omega0:.2f} → {a['omega_final']:.2f} rad/s")
    print(f"  Emec0 → Emec1  = {a['Emec0']:.1f} → {a['Emec1']:.1f} J")
    print(f"  ΔE_mec         = {a['dEmec']:.1f} J")
    print(f"  W_g            = {a['W_g']:.1f} J")
    print(f"  W_fric         = {a['W_fric']:.1f} J")
    print(f"  W_ext (neto)   = {a['W_ext']:.1f} J")
    print(f"  residual ΔEmec = {a['residual']:.3f} J   (ΔEmec - W_fric - W_ext → 0)")
    print(f"  residual ΔEkin = {a['residual_kin']:.3f} J")
    print(f"  residual ΔEpot = {a['residual_pot']:.3f} J   (ΔEpot + W_g → 0)")
    print(f"  E_fric_diss    = {a['E_fric_diss']:.1f} J")
    print(f"  E_util         = {a['E_util']:.1f} J   (Pμ={a['P_util_mean']:.1f} W)")
    print(f"  E_heat         = {a['E_heat']:.1f} J")
    print(f"  E_ext+         = {a['E_ext_pos']:.1f} J   (Pμ={a['P_ext_mean']:.1f} W)")
    print(f"  drain (-ΔEmec) = {a['drain']:.1f} J")
    print(f"  paid=E_ext++drain = {a['paid']:.1f} J")
    print(f"  η_paid=E_util/paid = {a['eta_paid']:.3f}   (≤ η_reg={p.eta_reg})")
    print(f"  figura         = {out}")
    print()
    if abs(a["residual"]) > 5.0:
        print("  AVISO: residual mecánico grande → revisar integración.")
    if a["eta_paid"] > p.eta_reg + 0.05:
        print("  AVISO: η_paid > η_reg → revisar contabilidad.")
    print("  Mensaje: E_útil sale de vaciar inercia + regen de fricción;")
    print("  W_g mueve E entre cinética y potencial (no crea stock neto).")
    print("  Kilómetro = convertidor/batería, no motor perpetuo.")


if __name__ == "__main__":
    main()
