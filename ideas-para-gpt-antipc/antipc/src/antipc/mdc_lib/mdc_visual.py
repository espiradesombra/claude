"""
MDC visual — animación dos trenes X/Y + restos con hipotenusa y colisiones.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from .mdc_analyze import (
    _is_proper_pair,
    analyze,
    diofantic_x,
    diofantic_y,
    factors_from_xy,
)


def _curve_xy(n: int, x_min: int, x_max: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs = np.linspace(x_min, x_max, max(200, (x_max - x_min) * 40))
    ys = []
    valid_x = []
    for xf in xs:
        y = diofantic_y(n, Fraction(xf).limit_denominator(10_000))
        if y is None:
            continue
        den = float(4 * xf + 6)
        if abs(den) < 1e-9:
            continue
        valid_x.append(xf)
        ys.append(float(y))
    return np.array(valid_x), np.array(ys), xs


def _collision_points(n: int, proper_only: bool) -> list[dict]:
    r = analyze(n)
    pts: list[dict] = []
    seen: set[tuple[int, int]] = set()
    for c in r.train_x.collisions:
        if proper_only and not _is_proper_pair(c.s, c.t, n):
            continue
        key = (int(c.x) if isinstance(c.x, int) else -999, int(c.y) if isinstance(c.y, int) else -999)
        if key in seen:
            continue
        seen.add(key)
        pts.append({"x": c.x, "y": c.y, "s": c.s, "t": c.t, "k": c.k})
    return pts


def _view_window(n: int, collisions: list[dict]) -> tuple[int, int, int, int]:
    sn = math.isqrt(n)
    if collisions:
        xs = [p["x"] for p in collisions if isinstance(p["x"], int)]
        ys = [p["y"] for p in collisions if isinstance(p["y"], int)]
        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        span = max(8, max(abs(x - cx) for x in xs) + 3, max(abs(y - cy) for y in ys) + 3)
        return (
            int(cx - span),
            int(cx + span),
            int(cy - span),
            int(cy + span),
        )
    return (-2, sn + 5, -2, sn + 5)


def animate_trains(
    n: int,
    *,
    proper_only: bool = True,
    save_path: str | Path | None = None,
    show: bool = True,
    fps: int = 8,
    dpi: int = 100,
) -> Path | None:
    collisions = _collision_points(n, proper_only)
    x_lo, x_hi, y_lo, y_hi = _view_window(n, collisions)

    cx, cy, _ = _curve_xy(n, x_lo, x_hi)
    x_ints = list(range(x_lo, x_hi + 1))
    y_from_x = []
    valid_xi = []
    for xi in x_ints:
        y = diofantic_y(n, Fraction(xi))
        if y is None or y.denominator != 1:
            continue
        valid_xi.append(xi)
        y_from_x.append(int(y))

    y_ints = list(range(y_lo, y_hi + 1))
    x_from_y = []
    valid_yi = []
    for yi in y_ints:
        x = diofantic_x(n, Fraction(yi))
        if x is None or x.denominator != 1:
            continue
        valid_yi.append(yi)
        x_from_y.append(int(x))

    coll_xy = {
        (p["x"], p["y"])
        for p in collisions
        if isinstance(p["x"], int) and isinstance(p["y"], int)
    }

    fig, (ax_mdc, ax_rest) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"MDC dos trenes — n={n}  |  (2x+3)(2y+3)=n", fontsize=11)

    # --- Panel MDC ---
    ax_mdc.plot(cx, cy, color="#4a90d9", lw=1.5, alpha=0.7, label="Recta diofántica")
    ax_mdc.axhline(0, color="#ccc", lw=0.5)
    ax_mdc.axvline(0, color="#ccc", lw=0.5)
    ax_mdc.set_xlim(x_lo, x_hi)
    ax_mdc.set_ylim(y_lo, y_hi)
    ax_mdc.set_xlabel("x  (tren X →)")
    ax_mdc.set_ylabel("y  (tren Y →)")
    ax_mdc.grid(True, alpha=0.25)
    ax_mdc.set_aspect("equal", adjustable="box")

    for p in collisions:
        if not isinstance(p["x"], int) or not isinstance(p["y"], int):
            continue
        ax_mdc.scatter(
            [p["x"]],
            [p["y"]],
            s=220,
            marker="*",
            c="#e74c3c",
            zorder=5,
            edgecolors="black",
            linewidths=0.5,
        )
        ax_mdc.annotate(
            f"{p['s']}×{p['t']}",
            (p["x"], p["y"]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
            color="#c0392b",
        )

    train_x_dot, = ax_mdc.plot([], [], "o", color="#2ecc71", ms=12, label="Tren X", zorder=6)
    train_y_dot, = ax_mdc.plot([], [], "s", color="#9b59b6", ms=10, label="Tren Y", zorder=6)
    trail_x, = ax_mdc.plot([], [], "-", color="#2ecc71", alpha=0.35, lw=1)
    trail_y, = ax_mdc.plot([], [], "-", color="#9b59b6", alpha=0.35, lw=1)
    status = ax_mdc.text(0.02, 0.98, "", transform=ax_mdc.transAxes, va="top", fontsize=9)
    ax_mdc.legend(loc="upper right", fontsize=8)

    # --- Panel restos + hipotenusa ---
    i_cap = min(n - 1, max(50, math.isqrt(n) * 4))
    i_vals = np.arange(1, i_cap + 1)
    y_rest = n % i_vals
    ax_rest.plot(i_vals, y_rest, color="#34495e", lw=1.2, label="n mod i")
    ax_rest.plot([0, i_cap], [0, i_cap], "k--", lw=1, label="Hipotenusa")
    for slope in (1, 2, 3):
        ax_rest.plot(i_vals, slope * i_vals, "--", alpha=0.4, lw=0.8, label=f"p={slope}")
    ax_rest.set_xlim(0, i_cap)
    ax_rest.set_ylim(0, i_cap)
    ax_rest.set_title("Restos + hipotenusa")
    ax_rest.set_xlabel("i")
    ax_rest.set_ylabel("n mod i")
    ax_rest.legend(fontsize=7, loc="upper left")
    ax_rest.grid(True, alpha=0.2)

    frames = max(len(x_ints), len(y_ints), 1)
    hist_x: list[float] = []
    hist_y: list[float] = []

    def update(frame: int):
        if x_ints:
            xi = x_ints[min(frame, len(x_ints) - 1)]
            yf = diofantic_y(n, Fraction(xi))
            if yf is not None:
                yfl = float(yf)
                train_x_dot.set_data([xi], [yfl])
                hist_x.append(xi)
                hist_y.append(yfl)
                trail_x.set_data(hist_x, hist_y)

        if y_ints:
            yi = y_ints[min(frame, len(y_ints) - 1)]
            xf = diofantic_x(n, Fraction(yi))
            if xf is not None:
                xfl = float(xf)
                train_y_dot.set_data([xfl], [yi])

        hit = ""
        if x_ints and y_ints:
            fi = min(frame, len(x_ints) - 1)
            fj = min(frame, len(y_ints) - 1)
            xi = x_ints[fi]
            yi = y_ints[fj]
            if (xi, yi) in coll_xy:
                pair = factors_from_xy(Fraction(xi), Fraction(yi), n)
                hit = f"  COLISIÓN {pair[0]}×{pair[1]}"
            elif frame < len(x_ints):
                y = diofantic_y(n, Fraction(xi))
                if y and y.denominator == 1 and (xi, int(y)) in coll_xy:
                    hit = "  COLISIÓN (tren X)"
        status.set_text(
            f"frame {frame+1}/{frames}{hit}"
        )
        return train_x_dot, train_y_dot, trail_x, trail_y, status

    anim = FuncAnimation(fig, update, frames=frames, interval=1000 // fps, blit=False)

    out: Path | None = None
    if save_path:
        out = Path(save_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.suffix.lower() in (".gif", ".webp"):
            anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
        else:
            out = out.with_suffix(".gif")
            anim.save(str(out), writer=PillowWriter(fps=fps), dpi=dpi)
        plt.close(fig)
    elif show:
        plt.tight_layout()
        plt.show()
    else:
        plt.close(fig)

    return out


def launch_gui(n: int, proper_only: bool = True, fps: int = 8) -> None:
    """Tkinter + animación dos trenes en vivo."""
    import tkinter as tk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    collisions = _collision_points(n, proper_only)
    x_lo, x_hi, y_lo, y_hi = _view_window(n, collisions)
    cx, cy, _ = _curve_xy(n, x_lo, x_hi)
    coll_xy = {
        (p["x"], p["y"])
        for p in collisions
        if isinstance(p["x"], int) and isinstance(p["y"], int)
    }
    x_ints = list(range(x_lo, x_hi + 1))
    y_ints = list(range(y_lo, y_hi + 1))
    frames = max(len(x_ints), len(y_ints), 1)
    hist_x: list[float] = []
    hist_y: list[float] = []

    root = tk.Tk()
    root.title(f"MDC dos trenes — n={n}")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(cx, cy, color="#4a90d9", lw=1.5, label="diofántica")
    for p in collisions:
        if isinstance(p["x"], int) and isinstance(p["y"], int):
            ax.scatter([p["x"]], [p["y"]], s=220, marker="*", c="#e74c3c", zorder=5)
            ax.annotate(f"{p['s']}×{p['t']}", (p["x"], p["y"]), fontsize=8, xytext=(5, 5),
                        textcoords="offset points")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("x  (tren X)")
    ax.set_ylabel("y  (tren Y)")
    ax.set_title(f"(2x+3)(2y+3)={n}")
    ax.grid(True, alpha=0.3)
    train_x_dot, = ax.plot([], [], "o", color="#2ecc71", ms=11, label="Tren X")
    train_y_dot, = ax.plot([], [], "s", color="#9b59b6", ms=9, label="Tren Y")
    trail_x, = ax.plot([], [], "-", color="#2ecc71", alpha=0.4)
    status = ax.text(0.02, 0.98, "", transform=ax.transAxes, va="top", fontsize=9)
    ax.legend(loc="upper right", fontsize=8)

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def update(frame: int) -> None:
        if x_ints:
            xi = x_ints[min(frame, len(x_ints) - 1)]
            yf = diofantic_y(n, Fraction(xi))
            if yf is not None:
                yfl = float(yf)
                train_x_dot.set_data([xi], [yfl])
                hist_x.append(xi)
                hist_y.append(yfl)
                trail_x.set_data(hist_x, hist_y)
        if y_ints:
            yi = y_ints[min(frame, len(y_ints) - 1)]
            xf = diofantic_x(n, Fraction(yi))
            if xf is not None:
                train_y_dot.set_data([float(xf)], [yi])
        hit = ""
        if x_ints:
            xi = x_ints[min(frame, len(x_ints) - 1)]
            y = diofantic_y(n, Fraction(xi))
            if y and y.denominator == 1 and (xi, int(y)) in coll_xy:
                pair = factors_from_xy(Fraction(xi), Fraction(int(y)), n)
                if pair:
                    hit = f"  COLISIÓN {pair[0]}×{pair[1]}"
        status.set_text(f"t={frame + 1}/{frames}{hit}")
        canvas.draw()

    def tick() -> None:
        update(tick.frame % frames)
        tick.frame += 1
        root.after(max(50, 1000 // fps), tick)

    tick.frame = 0
    tick()
    root.mainloop()