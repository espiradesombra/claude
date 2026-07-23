"""
gemelo_quijote_3vs7_balance.py
==============================
Quijote — gemelo 3 vs 7 palas con balance energético cerrado (1ª ley).

Modelo (didáctico, escala reducida tipo NREL-ideas):
  Rotor con inercia base J_G + N masas m_q en radios r_k(t).
  Ball: r_k = r0 + A * sin(ω_rot * t + 2π k / N)   [modo ball]
  o control estático / por signo de cos(θ_k)         [modo hurto-fase]

  J(t) = J_G + m_q * Σ r_k²
  J ω̇ + ω J̇ = T_aero - T_gen - T_act   (patinadora)

  T_aero ~ k_aero * v_wind² * sign-ish (simple: T_aero = T0 - b*ω  torque drive)
  T_gen  = K_g * ω                     (extracción generador)
  T_act  = trabajo del actuador de masas (estimado por |ṙ| · F_rad)

Balance:
  ΔE_cin + ΔE_act_store ≈ W_aero + W_gen + W_act + residual
  E_útil ~ -W_gen (si W_gen < 0 en convención... usamos signos claros)

Convención potencias:
  P_aero > 0  inyecta al rotor
  P_gen  > 0  extrae del rotor (útil a red)
  P_act  > 0  actuador gasta (entra energía mecánica al sistema de masas)

Uso:
  python gemelo_quijote_3vs7_balance.py
  python gemelo_quijote_3vs7_balance.py --T 40 --wind 11
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
class QParams:
    N: int = 3
    m_q: float = 4.0          # kg por pala (escala toy; paper usa ~332 kg fluido)
    J_G: float = 8.0e4        # kg m²  (toy; NREL real ~1e7)
    r_min: float = 5.0
    r_max: float = 25.0       # toy radius tip smaller for speed
    mode: str = "ball"        # ball | phase
    # aerodinámica simple
    T_drive: float = 3500.0   # N·m  par de empuje constante efectivo
    b_aero: float = 80.0      # amortiguación aero ~ b ω
    # generador
    K_gen: float = 120.0      # T_gen = K_gen * ω
    # actuador: coste ≈ m * |a_rad| * |ṙ| * r_scale (proxy)
    eta_act: float = 0.55     # eficiencia actuador (gasto > trabajo ideal)
    # sim
    omega0: float = 1.2
    dt: float = 0.002
    T: float = 30.0
    A_frac: float = 0.85      # amplitud ball como fracción de (rmax-rmin)/2


@dataclass
class QState:
    t: float = 0.0
    omega: float = 1.2
    theta_rotor: float = 0.0
    r: np.ndarray = field(default_factory=lambda: np.zeros(3))
    # trabajos
    W_aero: float = 0.0
    W_gen: float = 0.0
    W_act: float = 0.0
    W_pat: float = 0.0   # integral -ω J̇ * ω?  better track via energy
    E_util: float = 0.0  # energía a red ≈ integral T_gen * ω
    E_act_cost: float = 0.0
    hist_t: List[float] = field(default_factory=list)
    hist_omega: List[float] = field(default_factory=list)
    hist_J: List[float] = field(default_factory=list)
    hist_Ekin: List[float] = field(default_factory=list)
    hist_P_gen: List[float] = field(default_factory=list)
    hist_P_act: List[float] = field(default_factory=list)
    hist_P_aero: List[float] = field(default_factory=list)


def r_command(t: float, omega: float, theta_rotor: float, k: int, p: QParams) -> float:
    r0 = 0.5 * (p.r_min + p.r_max)
    amp = p.A_frac * 0.5 * (p.r_max - p.r_min)
    if p.mode == "ball":
        # ball en el frame del rotor: fase por pala
        phase = theta_rotor + 2 * math.pi * k / p.N
        return float(np.clip(r0 + amp * math.sin(phase), p.r_min, p.r_max))
    # phase / "hurto": fuera cuando baja (cos>0 si θ=0 horizontal... toy)
    blade_angle = theta_rotor + 2 * math.pi * k / p.N
    if math.cos(blade_angle) > 0:
        return p.r_max
    return p.r_min


def J_of(r: np.ndarray, p: QParams) -> float:
    return p.J_G + p.m_q * float(np.sum(r ** 2))


def step(st: QState, p: QParams) -> None:
    # target radii
    r_new = np.array([r_command(st.t, st.omega, st.theta_rotor, k, p) for k in range(p.N)])
    # rate limited slide
    v_max = 8.0  # m/s
    dr_max = v_max * p.dt
    r_next = st.r + np.clip(r_new - st.r, -dr_max, dr_max)

    J0 = J_of(st.r, p)
    J1 = J_of(r_next, p)
    dJ_dt = (J1 - J0) / p.dt

    # torques
    T_aero = p.T_drive - p.b_aero * st.omega
    if T_aero < 0:
        T_aero = 0.1 * T_aero  # no invertir loco
    T_gen = p.K_gen * max(st.omega, 0.0)
    # actuador: trabajo ideal de mover masas en campo centrífugo ~ m ω² r ṙ
    # P_act_ideal = Σ m * ω² * r * ṙ   (signo: si ṙ>0 hay que "pagar" contra o a favor)
    r_dot = (r_next - st.r) / p.dt
    P_cent = float(np.sum(p.m_q * (st.omega ** 2) * st.r * r_dot))
    # coste actuador: solo pagamos cuando trabajamos en contra de la fuerza neta efectiva
    # proxy: |P_cent| / eta si P_cent>0 (extender contra... simplified: always cost |ṙ| * F)
    P_act_cost = 0.0
    for k in range(p.N):
        F_cent = p.m_q * st.omega ** 2 * st.r[k]  # hacia fuera
        # fuerza de control opuesta al movimiento no deseado
        # trabajo del actuador ≈ |F_cent * ṙ| cuando retraemos (ṙ<0) o más
        P_act_cost += abs(F_cent * r_dot[k]) / max(p.eta_act, 0.1)
    # torque equivalente del actuador sobre el eje: no directo; entra como pérdida de potencia
    # en la eq de rotor usamos solo dJ/dt term; actuador se cuenta aparte en energía
    T_net = T_aero - T_gen
    # eq: J ω̇ = T_net - ω dJ/dt
    omega_dot = (T_net - st.omega * dJ_dt) / max(J0, 1.0)

    # potencias
    P_aero = T_aero * st.omega
    P_gen = T_gen * st.omega
    P_pat = st.omega * dJ_dt * st.omega  # ω² dJ/dt  (potencia del término patinadora)

    st.W_aero += P_aero * p.dt
    st.W_gen += P_gen * p.dt
    st.W_act += P_act_cost * p.dt
    st.E_util += P_gen * p.dt
    st.E_act_cost += P_act_cost * p.dt

    # integrate
    st.omega = max(0.05, st.omega + omega_dot * p.dt)
    st.theta_rotor += st.omega * p.dt
    st.r = r_next
    st.t += p.dt

    Ekin = 0.5 * J1 * st.omega ** 2
    st.hist_t.append(st.t)
    st.hist_omega.append(st.omega)
    st.hist_J.append(J1)
    st.hist_Ekin.append(Ekin)
    st.hist_P_gen.append(P_gen)
    st.hist_P_act.append(P_act_cost)
    st.hist_P_aero.append(P_aero)


def simulate(p: QParams) -> QState:
    r0 = 0.5 * (p.r_min + p.r_max)
    st = QState(omega=p.omega0, r=np.full(p.N, r0))
    n = int(p.T / p.dt)
    for _ in range(n):
        step(st, p)
    return st


def analyze(st: QState, p: QParams) -> dict:
    E0, E1 = st.hist_Ekin[0], st.hist_Ekin[-1]
    dE = E1 - E0
    # ΔE_cin ≈ W_aero - W_gen - W_pat_effect ...
    # Contabilidad práctica didáctica:
    # paid = E_act_cost
    # gained = E_util
    # also aero injects W_aero
    # residual: dE - (W_aero - W_gen - W_act_mech)
    # usamos: dE ?= W_aero - W_gen - W_pat_store
    # proxy residual using aero and gen only + act as external:
    residual = dE - (st.W_aero - st.W_gen) + 0.0
    # mejor: energy balance
    # E1 - E0 = W_aero - W_gen - W_to_act_channel
    # no separamos W_to_act en el eje; reportamos métricas claras
    return {
        "N": p.N,
        "mode": p.mode,
        "E0": E0,
        "E1": E1,
        "dE": dE,
        "W_aero": st.W_aero,
        "W_gen": st.W_gen,
        "E_util": st.E_util,
        "E_act_cost": st.E_act_cost,
        "P_gen_mean": st.E_util / max(st.t, 1e-9),
        "P_act_mean": st.E_act_cost / max(st.t, 1e-9),
        "P_aero_mean": st.W_aero / max(st.t, 1e-9),
        "omega_mean": float(np.mean(st.hist_omega)),
        "J_mean": float(np.mean(st.hist_J)),
        "J_std": float(np.std(st.hist_J)),
        "net_power": (st.E_util - st.E_act_cost) / max(st.t, 1e-9),
        "ratio_util_act": st.E_util / max(st.E_act_cost, 1e-9),
    }


def run_pair(p_base: QParams) -> Tuple[QState, QState, dict, dict]:
    p3 = QParams(**{**p_base.__dict__, "N": 3})
    p7 = QParams(**{**p_base.__dict__, "N": 7})
    s3, s7 = simulate(p3), simulate(p7)
    return s3, s7, analyze(s3, p3), analyze(s7, p7)


def plot_compare(s3: QState, s7: QState, a3: dict, a7: dict, out: str) -> None:
    t3, t7 = np.array(s3.hist_t), np.array(s7.hist_t)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.4), facecolor="#0f1419")
    for ax in axes.ravel():
        ax.set_facecolor("#1a222c")
        ax.tick_params(colors="#c8d0d8")
        for sp in ax.spines.values():
            sp.set_color("#3a4550")
        ax.xaxis.label.set_color("#a8b0b8")
        ax.yaxis.label.set_color("#a8b0b8")
        ax.title.set_color("#e8eef4")
        ax.grid(True, alpha=0.2, color="#5a6570")

    axes[0, 0].plot(t3, s3.hist_omega, color="#5ec8ff", lw=1.2, label="N=3")
    axes[0, 0].plot(t7, s7.hist_omega, color="#ffb86c", lw=1.2, label="N=7", alpha=0.9)
    axes[0, 0].set_title("ω(t)")
    axes[0, 0].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8")

    axes[0, 1].plot(t3, s3.hist_J, color="#5ec8ff", lw=1.1, label="J N=3")
    axes[0, 1].plot(t7, s7.hist_J, color="#ffb86c", lw=1.1, label="J N=7")
    axes[0, 1].set_title("J(t)  (ball → casi plano si r0,A bien)")
    axes[0, 1].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=8)

    axes[1, 0].plot(t3, s3.hist_P_gen, color="#5ec8ff", lw=1.0, label="P_gen N=3")
    axes[1, 0].plot(t7, s7.hist_P_gen, color="#ffb86c", lw=1.0, label="P_gen N=7")
    axes[1, 0].plot(t3, s3.hist_P_act, color="#5ec8ff", lw=1.0, ls="--", alpha=0.7, label="P_act N=3")
    axes[1, 0].plot(t7, s7.hist_P_act, color="#ffb86c", lw=1.0, ls="--", alpha=0.7, label="P_act N=7")
    axes[1, 0].set_title("Potencias gen (sólido) vs actuador (dashed)")
    axes[1, 0].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=7)

    labels = ["P_gen μ", "P_act μ", "P_net μ\n(gen-act)"]
    x = np.arange(len(labels))
    w = 0.35
    v3 = [a3["P_gen_mean"], a3["P_act_mean"], a3["net_power"]]
    v7 = [a7["P_gen_mean"], a7["P_act_mean"], a7["net_power"]]
    axes[1, 1].bar(x - w / 2, v3, w, color="#5ec8ff", label="N=3", edgecolor="#3a4550")
    axes[1, 1].bar(x + w / 2, v7, w, color="#ffb86c", label="N=7", edgecolor="#3a4550")
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels(labels)
    axes[1, 1].set_ylabel("W")
    axes[1, 1].set_title("Comparativa media  ·  toy scale (no NREL full)")
    axes[1, 1].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8")

    fig.suptitle(
        "Quijote 3 vs 7 — gemelo didáctico con coste de actuador\n"
        f"mode={a3['mode']}  ·  net3={a3['net_power']:.0f} W  net7={a7['net_power']:.0f} W  ·  "
        f"ratio util/act: 3→{a3['ratio_util_act']:.2f}  7→{a7['ratio_util_act']:.2f}",
        color="#f0f4f8", fontsize=11, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=30.0)
    ap.add_argument("--mode", choices=["ball", "phase"], default="ball")
    ap.add_argument("--mq", type=float, default=4.0)
    args = ap.parse_args()

    base = QParams(T=args.T, mode=args.mode, m_q=args.mq)
    s3, s7, a3, a7 = run_pair(base)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "quijote_3vs7_balance.png")
    plot_compare(s3, s7, a3, a7, out)

    def dump(tag, a):
        print(f"--- N={a['N']} mode={a['mode']} ---")
        print(f"  ω_mean        = {a['omega_mean']:.3f} rad/s")
        print(f"  J_mean ± std  = {a['J_mean']:.0f} ± {a['J_std']:.0f} kg m²")
        print(f"  P_aero μ      = {a['P_aero_mean']:.1f} W")
        print(f"  P_gen  μ      = {a['P_gen_mean']:.1f} W")
        print(f"  P_act  μ      = {a['P_act_mean']:.1f} W")
        print(f"  P_net  μ      = {a['net_power']:.1f} W   (gen - act)")
        print(f"  E_util / E_act= {a['ratio_util_act']:.3f}")
        print(f"  ΔE_cin        = {a['dE']:.1f} J")

    print("=" * 62)
    print("QUIJOTE 3 vs 7 — balance con coste de actuador (toy)")
    print("=" * 62)
    dump("3", a3)
    dump("7", a7)
    print()
    print(f"  ΔP_net (7-3)  = {a7['net_power'] - a3['net_power']:+.1f} W")
    print(f"  figura        = {out}")
    print()
    print("  Lectura: más palas → ball más continuo (menos rizado J);")
    print("  el actuador SIEMPRE cuesta. Si P_act > beneficio, no hay 'hurto' neto.")
    print("  Escala toy (no sustituye gemell_3vs7_rigoros NREL).")


if __name__ == "__main__":
    main()
