#!/usr/bin/env python3
"""
Restos + hipotenusa + animación MDC dos trenes.
Reexporta visualización desde AntiPC mdc_lib.
"""

from __future__ import annotations

import argparse
import os
import sys

_ANTIPC = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "ideas-para-gpt-antipc", "antipc", "src", "antipc")
)
if os.path.isdir(_ANTIPC):
    sys.path.insert(0, _ANTIPC)

try:
    from mdc_lib.mdc_visual import animate_trains, launch_gui
except ImportError:
    import tkinter as tk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    import matplotlib.pyplot as plt
    import numpy as np

    def plot_restos(n: int):
        i_vals = np.arange(1, min(n, 120))
        y_vals = n % i_vals
        fig, ax = plt.subplots()
        ax.plot(i_vals, y_vals, label="n mod i")
        ax.plot([0, i_vals[-1]], [0, i_vals[-1]], "k--", label="Hipotenusa")
        ax.set_title(f"Restos de n={n}")
        ax.legend()
        return fig

    def launch_gui(n: int, proper_only: bool = True) -> None:
        root = tk.Tk()
        root.title("Visualización de restos")
        fig = plot_restos(n)
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack()
        root.mainloop()

    def animate_trains(n, **kwargs):
        raise SystemExit("Instala AntiPC mdc_lib o matplotlib completo")


def main() -> None:
    p = argparse.ArgumentParser(description="MDC visual — restos + trenes")
    p.add_argument("n", type=int, nargs="?", default=97)
    p.add_argument("-o", "--output", help="GIF de salida")
    p.add_argument("--gui", action="store_true", help="Solo ventana Tk")
    p.add_argument("--proper", action="store_true")
    p.add_argument("--fps", type=int, default=8)
    args = p.parse_args()

    if args.gui:
        launch_gui(args.n, proper_only=args.proper)
        return

    out = args.output or f"mdc_trains_{args.n}.gif"
    path = animate_trains(
        args.n,
        proper_only=args.proper,
        save_path=out,
        show=False,
        fps=args.fps,
    )
    print(f"OK: {path}")


if __name__ == "__main__":
    main()