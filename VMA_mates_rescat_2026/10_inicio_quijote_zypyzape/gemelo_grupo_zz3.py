"""
gemelo_grupo_zz3.py
===================
Grupo de 3 unidades ZypyZape (roles gen/thr/buf) + bus con H más realista.

Mejoras de calibración:
  - H_sys y S_base coherentes con potencia del grupo
  - escalón de carga relativo (fracción de P_load)
  - droop / amortiguación para nadir creíble (~décimas de Hz, no 4 Hz)

Uso:
  python gemelo_grupo_zz3.py
  python gemelo_grupo_zz3.py --compare-H
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


@dataclass
class GParams:
    n: int = 3
    J: float = 2.0e4           # kg m² por máquina (más inercia local)
    K_coup: float = 120.0
    c_gen: float = 55.0
    T_thr: float = 80.0
    b_buf: float = 40.0
    f_zz: float = 0.35
    # red — calibrada
    H: float = 5.0             # s  inercia sistema
    f0: float = 50.0
    S_base: float = 2.0e5      # VA  base del bus toy
    D: float = 1.5             # amortiguación (1/s efectiva en Δf)
    P_load: float = 8.0e3      # W
    drop_at: float = 12.0
    drop_frac: float = 0.15    # +15% carga (antes era salto enorme)
    # sim
    dt: float = 0.002
    T: float = 35.0
    omega0: float = 1.6
    use_zz: bool = True        # False = sin ciclo de roles (solo acoplamiento)


@dataclass
class GState:
    t: float = 0.0
    theta: np.ndarray = field(default_factory=lambda: np.zeros(3))
    omega: np.ndarray = field(default_factory=lambda: np.ones(3) * 1.6)
    df: float = 0.0
    E_gen: float = 0.0
    hist_t: List[float] = field(default_factory=list)
    hist_omega: List[List[float]] = field(default_factory=list)
    hist_f: List[float] = field(default_factory=list)
    hist_Pgen: List[float] = field(default_factory=list)
    hist_Hzz: List[float] = field(default_factory=list)
    hist_roles: List[List[str]] = field(default_factory=list)


def roles(t: float, n: int, f_zz: float, use_zz: bool) -> List[str]:
    if not use_zz:
        return ["gen"] * n  # todos generan un poco (baseline)
    order = ["gen", "thr", "buf"]
    base = int(math.floor(n * f_zz * t)) % n
    return [order[(base + i) % n] for i in range(n)]


def step(st: GState, p: GParams) -> None:
    n = p.n
    r = roles(st.t, n, p.f_zz, p.use_zz)
    omega_ref = float(np.mean(st.omega))
    P_gen_sum = 0.0
    P_mech_to_bus = 0.0

    T_net = np.zeros(n)
    for i in range(n):
        T_c = 0.0
        for j in range(n):
            if i == j:
                continue
            T_c += p.K_coup * math.sin(st.theta[j] - st.theta[i])
        T_role = 0.0
        if r[i] == "gen":
            T_role = -p.c_gen * st.omega[i]
            P = max(0.0, -T_role * st.omega[i])
            P_gen_sum += P
            P_mech_to_bus += P
        elif r[i] == "thr":
            # thr toma energía del bus (consume) y acelera rotor
            T_role = p.T_thr * (1.0 - 0.3 * max(0.0, -st.df))  # más thr si f baja
            P_mech_to_bus -= T_role * st.omega[i] * 0.5  # parte del thr viene del bus
        else:
            T_role = -p.b_buf * (st.omega[i] - omega_ref)
        # pequeño drive base para no parar
        T_drive = 25.0
        T_net[i] = T_c + T_role + T_drive

    for i in range(n):
        alpha = T_net[i] / p.J
        st.omega[i] = max(0.3, min(4.0, st.omega[i] + alpha * p.dt))
        st.theta[i] += st.omega[i] * p.dt

    P_load = p.P_load
    if st.t >= p.drop_at:
        P_load *= (1.0 + p.drop_frac)

    H_zz = 0.5 * p.J * float(np.sum(st.omega ** 2)) / p.S_base
    H_eff = p.H + (H_zz if p.use_zz else 0.5 * H_zz)

    P_balance = P_mech_to_bus - P_load
    # 2 H S / f0 * d(Δf)/dt = P_balance / (algo) — forma:
    # df/dt = (f0/(2 H S)) * P_balance - D*df
    ddf = (p.f0 / (2.0 * H_eff * p.S_base)) * P_balance - p.D * st.df
    st.df += ddf * p.dt
    # clamp suave
    st.df = float(np.clip(st.df, -2.5, 1.0))

    st.E_gen += P_gen_sum * p.dt
    st.t += p.dt
    st.hist_t.append(st.t)
    st.hist_omega.append(st.omega.copy().tolist())
    st.hist_f.append(p.f0 + st.df)
    st.hist_Pgen.append(P_gen_sum)
    st.hist_Hzz.append(H_zz)
    st.hist_roles.append(r)


def simulate(p: GParams) -> GState:
    st = GState(
        theta=np.linspace(0, 2 * math.pi, p.n, endpoint=False),
        omega=np.ones(p.n) * p.omega0 + 0.05 * np.random.randn(p.n),
    )
    n = int(p.T / p.dt)
    for _ in range(n):
        step(st, p)
    return st


def metrics(st: GState, p: GParams) -> dict:
    f = np.array(st.hist_f)
    i0 = int(p.drop_at / p.dt)
    post = f[i0:] if i0 < len(f) else f
    return {
        "f_min": float(np.min(f)),
        "f_max": float(np.max(f)),
        "nadir": float(np.min(post)),
        "rocof_max": float(np.max(np.abs(np.diff(f) / p.dt))) if len(f) > 2 else 0.0,
        "Hzz_mean": float(np.mean(st.hist_Hzz)),
        "E_gen": st.E_gen,
    }


def plot_report(st_zz: GState, st_base: Optional[GState], p: GParams, out: str) -> None:
    t = np.array(st_zz.hist_t)
    om = np.array(st_zz.hist_omega)
    f = np.array(st_zz.hist_f)
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.3), facecolor="#0f1419")
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
        axes[0, 0].plot(t, om[:, i], color=cols[i], lw=1.1, label=f"ω{i}")
    axes[0, 0].set_title("ω_i(t) — ZZ ON")
    axes[0, 0].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=8)

    axes[0, 1].plot(t, f, color="#c792ea", lw=1.4, label="ZZ ON")
    if st_base is not None:
        axes[0, 1].plot(st_base.hist_t, st_base.hist_f, color="#ff6b6b", lw=1.2,
                        ls="--", label="ZZ OFF (baseline)")
    axes[0, 1].axvline(p.drop_at, color="#888", ls=":", lw=1)
    axes[0, 1].axhline(50.0, color="#555", ls=":", lw=0.8)
    axes[0, 1].set_title("frecuencia bus [Hz]")
    axes[0, 1].legend(facecolor="#1a222c", edgecolor="#3a4550", labelcolor="#c8d0d8", fontsize=8)

    axes[1, 0].plot(t, st_zz.hist_Hzz, color="#7dffa0", lw=1.2)
    axes[1, 0].set_title("H_ZZ(t) [s]  inercia sintética del grupo")
    axes[1, 0].set_xlabel("t [s]")

    role_map = {"gen": 2, "thr": 1, "buf": 0}
    R = np.array([[role_map[x] for x in row] for row in st_zz.hist_roles])
    im = axes[1, 1].imshow(R.T, aspect="auto", origin="lower",
                           extent=[t[0], t[-1], -0.5, p.n - 0.5],
                           cmap="coolwarm", interpolation="nearest")
    axes[1, 1].set_title("Roles 0=buf 1=thr 2=gen")
    axes[1, 1].set_xlabel("t [s]")
    fig.colorbar(im, ax=axes[1, 1], fraction=0.046)

    m = metrics(st_zz, p)
    extra = ""
    if st_base is not None:
        mb = metrics(st_base, p)
        extra = f"  ·  nadir OFF={mb['nadir']:.3f} / ON={m['nadir']:.3f} Hz"
    fig.suptitle(
        f"Grupo ZZ×3 calibrado  ·  H={p.H}s  drop=+{100*p.drop_frac:.0f}%  "
        f"nadir={m['nadir']:.3f} Hz  RoCoFμmax≈{m['rocof_max']:.3f} Hz/s"
        + extra,
        color="#f0f4f8", fontsize=10.5, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=35.0)
    ap.add_argument("--K", type=float, default=120.0)
    ap.add_argument("--H", type=float, default=5.0)
    ap.add_argument("--drop", type=float, default=0.15, help="fracción escalón carga")
    ap.add_argument("--compare-H", action="store_true", help="ZZ ON vs OFF")
    args = ap.parse_args()

    p = GParams(T=args.T, K_coup=args.K, H=args.H, drop_frac=args.drop, use_zz=True)
    st_zz = simulate(p)
    st_base = None
    if args.compare_H:
        p_off = GParams(T=args.T, K_coup=args.K, H=args.H, drop_frac=args.drop, use_zz=False)
        st_base = simulate(p_off)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "grupo_zz3_balance.png")
    plot_report(st_zz, st_base, p, out)
    m = metrics(st_zz, p)
    print("=" * 60)
    print("GRUPO ZYPYZAPE ×3 — calibrado")
    print("=" * 60)
    print(f"  H_sys={p.H} s  S_base={p.S_base:.0f} VA  drop=+{100*p.drop_frac:.0f}%")
    print(f"  f min/max     = {m['f_min']:.4f} / {m['f_max']:.4f} Hz")
    print(f"  nadir post    = {m['nadir']:.4f} Hz")
    print(f"  RoCoF max     ≈ {m['rocof_max']:.4f} Hz/s")
    print(f"  H_ZZ mean     = {m['Hzz_mean']:.4f} s")
    if st_base is not None:
        mb = metrics(st_base, p)
        print(f"  nadir ZZ OFF  = {mb['nadir']:.4f} Hz")
        print(f"  Δnadir (ON-OFF)= {m['nadir'] - mb['nadir']:+.4f} Hz")
    print(f"  figura        = {out}")


if __name__ == "__main__":
    main()
