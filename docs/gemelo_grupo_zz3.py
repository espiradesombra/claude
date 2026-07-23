"""
gemelo_grupo_zz3.py
===================
Grupo de 3 unidades estilo ZypyZape (roles gen / thr / buf cíclicos)
+ opcional "Kilómetro-like" buffer en cada nodo (E_buf).

Modelo de red minimalista (no CFD):
  Cada unidad i tiene ω_i, J fijo, acoplamiento Kuramoto-like en par:
    T_coup_i = K * Σ sin(θ_j - θ_i)
  Cicle ZZ: rol r_i(t) ∈ {gen, thr, buf} a 120°
    gen: extrae P_gen = c_g * ω²
    thr: inyecta par de "empuje" (simula recibir energía de red o de buf)
    buf: frena/acelera suave hacia ω_ref del grupo

  Swing de frecuencia de bus (1 bus):
    2 H S d(Δf)/dt = P_in - P_load - D Δf

Uso:
  python gemelo_grupo_zz3.py
  python gemelo_grupo_zz3.py --K 50 --T 40
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, field
from typing import List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class GParams:
    n: int = 3
    J: float = 5000.0
    K_coup: float = 80.0       # N·m
    c_gen: float = 40.0        # P = c ω²
    T_thr: float = 60.0        # N·m empuje en rol thr
    b_buf: float = 30.0        # freno buf ~ b (ω - ω_ref)
    f_zz: float = 0.35         # Hz rotación de roles
    # red
    H: float = 4.0
    S_base: float = 5e5        # VA toy
    D: float = 0.02
    P_load: float = 1.2e4      # W
    # sim
    dt: float = 0.002
    T: float = 35.0
    omega0: float = 1.5
    drop_at: float = 12.0      # s  escalón de carga
    drop_dP: float = 4.0e3     # W  extra load


@dataclass
class GState:
    t: float = 0.0
    theta: np.ndarray = field(default_factory=lambda: np.zeros(3))
    omega: np.ndarray = field(default_factory=lambda: np.ones(3) * 1.5)
    df: float = 0.0            # Hz desviación
    E_gen: float = 0.0
    E_thr: float = 0.0
    hist_t: List[float] = field(default_factory=list)
    hist_omega: List[List[float]] = field(default_factory=list)
    hist_f: List[float] = field(default_factory=list)
    hist_Pgen: List[float] = field(default_factory=list)
    hist_roles: List[List[str]] = field(default_factory=list)


def roles(t: float, n: int, f_zz: float) -> List[str]:
    order = ["gen", "thr", "buf"]
    base = int(math.floor(n * f_zz * t)) % n
    return [order[(base + i) % n] for i in range(n)]


def step(st: GState, p: GParams) -> None:
    n = p.n
    r = roles(st.t, n, p.f_zz)
    omega_ref = float(np.mean(st.omega))
    P_gen_sum = 0.0
    P_thr_sum = 0.0

    T_net = np.zeros(n)
    for i in range(n):
        # acoplamiento
        T_c = 0.0
        for j in range(n):
            if i == j:
                continue
            T_c += p.K_coup * math.sin(st.theta[j] - st.theta[i])
        T_role = 0.0
        if r[i] == "gen":
            # freno generador
            T_role = -p.c_gen * st.omega[i]
            P_gen_sum += max(0.0, -T_role * st.omega[i])
        elif r[i] == "thr":
            T_role = p.T_thr
            P_thr_sum += T_role * st.omega[i]
        else:
            T_role = -p.b_buf * (st.omega[i] - omega_ref)
        T_net[i] = T_c + T_role

    # integrar rotores
    for i in range(n):
        alpha = T_net[i] / p.J
        st.omega[i] = max(0.2, st.omega[i] + alpha * p.dt)
        st.theta[i] += st.omega[i] * p.dt

    # bus frecuencia (toy): P_gen_sum alimenta, thr se ve como consumo del bus
    P_load = p.P_load
    if st.t >= p.drop_at:
        P_load += p.drop_dP
    # inercia sintética del grupo
    H_zz = 0.5 * p.J * float(np.sum(st.omega ** 2)) / p.S_base
    H_eff = p.H + H_zz
    P_balance = P_gen_sum - P_thr_sum * 0.15 - P_load  # thr parcialmente de bus
    # swing in Δf (Hz): 2 H S / f0 * df/dt ≈ P  →  usamos f0=50
    f0 = 50.0
    ddf = (P_balance * f0) / (2.0 * H_eff * p.S_base) - p.D * st.df
    st.df += ddf * p.dt

    st.E_gen += P_gen_sum * p.dt
    st.E_thr += P_thr_sum * p.dt
    st.t += p.dt

    st.hist_t.append(st.t)
    st.hist_omega.append(st.omega.copy().tolist())
    st.hist_f.append(f0 + st.df)
    st.hist_Pgen.append(P_gen_sum)
    st.hist_roles.append(r)


def simulate(p: GParams) -> GState:
    st = GState(
        theta=np.linspace(0, 2 * math.pi, p.n, endpoint=False),
        omega=np.ones(p.n) * p.omega0,
    )
    n = int(p.T / p.dt)
    for _ in range(n):
        step(st, p)
    return st


def plot_report(st: GState, p: GParams, out: str) -> None:
    t = np.array(st.hist_t)
    om = np.array(st.hist_omega)
    f = np.array(st.hist_f)
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

    cols = ["#5ec8ff", "#ffb86c", "#7dffa0"]
    for i in range(p.n):
        axes[0, 0].plot(t, om[:, i], color=cols[i % 3], lw=1.1, label=f"ω{i}")
    axes[0, 0].set_title("ω_i(t) — 3 unidades ZZ")
    axes[0, 0].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=8)

    axes[0, 1].plot(t, f, color="#c792ea", lw=1.3)
    axes[0, 1].axvline(p.drop_at, color="#ff6b6b", ls="--", lw=1, label="escalón carga")
    axes[0, 1].set_title("frecuencia bus f(t) [Hz]")
    axes[0, 1].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=8)

    axes[1, 0].plot(t, st.hist_Pgen, color="#5ec8ff", lw=1.1)
    axes[1, 0].set_title("P_gen total del grupo [W]")
    axes[1, 0].set_xlabel("t [s]")

    # role timeline (numeric encode)
    role_map = {"gen": 2, "thr": 1, "buf": 0}
    R = np.array([[role_map[r] for r in row] for row in st.hist_roles])
    im = axes[1, 1].imshow(R.T, aspect="auto", origin="lower",
                           extent=[t[0], t[-1], -0.5, p.n - 0.5],
                           cmap="coolwarm", interpolation="nearest")
    axes[1, 1].set_title("Roles (0=buf, 1=thr, 2=gen)")
    axes[1, 1].set_xlabel("t [s]")
    axes[1, 1].set_ylabel("unidad")
    fig.colorbar(im, ax=axes[1, 1], fraction=0.046)

    fmin = float(np.min(f))
    fig.suptitle(
        f"Grupo ZZ ×3  ·  f_zz={p.f_zz} Hz  K={p.K_coup}  ·  f_min={fmin:.3f} Hz tras escalón",
        color="#f0f4f8", fontsize=11, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=35.0)
    ap.add_argument("--K", type=float, default=80.0)
    ap.add_argument("--fzz", type=float, default=0.35)
    args = ap.parse_args()
    p = GParams(T=args.T, K_coup=args.K, f_zz=args.fzz)
    st = simulate(p)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "grupo_zz3_balance.png")
    plot_report(st, p, out)
    f = np.array(st.hist_f)
    print("=" * 60)
    print("GRUPO ZYPYZAPE ×3 — roles gen/thr/buf + bus")
    print("=" * 60)
    print(f"  f_zz           = {p.f_zz} Hz")
    print(f"  K_coup         = {p.K_coup}")
    print(f"  f min / max    = {f.min():.4f} / {f.max():.4f} Hz")
    print(f"  nadir post-drop≈ {f[int(p.drop_at/p.dt):].min():.4f} Hz")
    print(f"  E_gen total    = {st.E_gen:.1f} J")
    print(f"  figura         = {out}")
    print()
    print("  Lectura: al escalón de carga, H_zz de los 3 rotores amortigua df/dt.")
    print("  Roles a 120° reparte quién genera y quién empuja (estilo XFI/ZZ).")


if __name__ == "__main__":
    main()
