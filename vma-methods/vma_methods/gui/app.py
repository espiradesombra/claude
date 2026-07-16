#!/usr/bin/env python3
"""VMA-METHODS GUI — Tkinter + matplotlib."""

from __future__ import annotations

import io
import sys
import time
import tkinter as tk
from tkinter import messagebox, scrolledtext, ttk

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from vma_methods import classic, cribas, criva, newton


class VmaMethodsApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        root.title("VMA-METHODS — Teoría de números")
        root.geometry("920x640")
        root.minsize(800, 520)

        header = ttk.Label(
            root,
            text="VMA-METHODS  ·  Víctor Manzanares Alberola  ·  cribas · Criva · Newton Rápido",
            font=("Segoe UI", 10, "bold"),
        )
        header.pack(fill=tk.X, padx=8, pady=6)

        nb = ttk.Notebook(root)
        nb.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        self._tab_cribas(nb)
        self._tab_criva(nb)
        self._tab_newton(nb)
        self._tab_classic(nb)

        status = ttk.Label(root, text="Listo", relief=tk.SUNKEN, anchor=tk.W)
        status.pack(fill=tk.X, side=tk.BOTTOM)
        self.status = status

    def _log(self, widget: scrolledtext.ScrolledText, text: str) -> None:
        widget.delete("1.0", tk.END)
        widget.insert(tk.END, text)
        self.status.config(text="OK")

    def _tab_cribas(self, nb: ttk.Notebook) -> None:
        frame = ttk.Frame(nb)
        nb.add(frame, text="Cribas")

        top = ttk.Frame(frame)
        top.pack(fill=tk.X, padx=8, pady=8)

        ttk.Label(top, text="Límite:").grid(row=0, column=0, sticky=tk.W)
        var_limit = tk.StringVar(value="5000")
        ttk.Entry(top, textvariable=var_limit, width=12).grid(row=0, column=1, padx=4)

        ttk.Label(top, text="Tipo:").grid(row=0, column=2, padx=(12, 0))
        var_tipo = tk.StringVar(value="6k")
        ttk.Combobox(
            top,
            textvariable=var_tipo,
            values=["desmemoriada", "6k", "hibrida"],
            width=14,
            state="readonly",
        ).grid(row=0, column=3, padx=4)

        out = scrolledtext.ScrolledText(frame, height=12, font=("Consolas", 9))
        out.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        fig = Figure(figsize=(7, 3), dpi=100)
        ax = fig.add_subplot(111)
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        def run_one() -> None:
            try:
                limit = int(var_limit.get())
            except ValueError:
                messagebox.showerror("Error", "Límite inválido")
                return
            mapping = {
                "desmemoriada": cribas.CribaDesmemoriada,
                "6k": cribas.CribaModular6k,
                "hibrida": cribas.CribaHibrida,
            }
            t0 = time.perf_counter()
            primes = mapping[var_tipo.get()](limit).run()
            ms = (time.perf_counter() - t0) * 1000
            self._log(
                out,
                f"Criba {var_tipo.get()} hasta {limit}\n"
                f"Primos: {len(primes)}\n"
                f"Tiempo: {ms:.2f} ms\n"
                f"Último: {primes[-1] if primes else '—'}\n"
                f"Muestra: {primes[:15]}",
            )

        def run_compare() -> None:
            try:
                limit = int(var_limit.get())
            except ValueError:
                messagebox.showerror("Error", "Límite inválido")
                return
            buf = io.StringIO()
            old = sys.stdout
            sys.stdout = buf
            try:
                cribas.comparar_cribas(limit=limit, verbose=True)
            finally:
                sys.stdout = old
            self._log(out, buf.getvalue())

            names, times_ms, counts = [], [], []
            ref = classic.pi_count(limit)
            for name, cls in [
                ("Desmemoriada", cribas.CribaDesmemoriada),
                ("6k±1", cribas.CribaModular6k),
                ("Híbrida", cribas.CribaHibrida),
            ]:
                t0 = time.perf_counter()
                n = len(cls(limit).run())
                times_ms.append((time.perf_counter() - t0) * 1000)
                names.append(name)
                counts.append(n)
            ax.clear()
            colors = ["#3498db", "#2ecc71", "#9b59b6"]
            bars = ax.bar(names, times_ms, color=colors, edgecolor="k")
            ax.set_ylabel("ms")
            ax.set_title(f"Benchmark cribas (π({limit})={ref})")
            for bar, c in zip(bars, counts):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height(),
                    f"{int(c)}p",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )
            fig.tight_layout()
            canvas.draw()

        btn_row = ttk.Frame(top)
        btn_row.grid(row=1, column=0, columnspan=4, pady=8, sticky=tk.W)
        ttk.Button(btn_row, text="Ejecutar criba", command=run_one).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_row, text="Comparar 3 cribas", command=run_compare).pack(side=tk.LEFT, padx=2)

    def _tab_criva(self, nb: ttk.Notebook) -> None:
        frame = ttk.Frame(nb)
        nb.add(frame, text="Criva")

        top = ttk.Frame(frame)
        top.pack(fill=tk.X, padx=8, pady=8)
        ttk.Label(top, text="x:").pack(side=tk.LEFT)
        var_x = tk.StringVar(value="100,500,1000,5000,10000")
        ttk.Entry(top, textvariable=var_x, width=40).pack(side=tk.LEFT, padx=4)

        out = scrolledtext.ScrolledText(frame, height=10, font=("Consolas", 9))
        out.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        fig = Figure(figsize=(7, 3.5), dpi=100)
        ax = fig.add_subplot(111)
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        def run() -> None:
            try:
                xs = [int(v.strip()) for v in var_x.get().split(",") if v.strip()]
            except ValueError:
                messagebox.showerror("Error", "Lista x inválida (separar por comas)")
                return
            buf = io.StringIO()
            old = sys.stdout
            sys.stdout = buf
            try:
                rows = criva.compare_criva_vs_pnt(xs, verbose=True)
            finally:
                sys.stdout = old
            self._log(out, buf.getvalue())

            ax.clear()
            x_plot = [r[0] for r in rows]
            criva_est = [r[1] for r in rows]
            pnt_est = [r[2] for r in rows]
            exact = [r[3] for r in rows]
            ax.plot(x_plot, exact, "o-", label="π(x) exacto", color="#2c3e50")
            ax.plot(x_plot, criva_est, "s--", label="Criva·x", color="#2ecc71")
            ax.plot(x_plot, pnt_est, "^:", label="PNT x/ln(x)", color="#e74c3c")
            ax.set_xlabel("x")
            ax.set_ylabel("conteo")
            ax.set_title("Criva vs PNT vs π(x)")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
            fig.tight_layout()
            canvas.draw()

        ttk.Button(top, text="Comparar", command=run).pack(side=tk.LEFT, padx=8)

    def _tab_newton(self, nb: ttk.Notebook) -> None:
        frame = ttk.Frame(nb)
        nb.add(frame, text="Newton Rápido")

        form = ttk.Frame(frame)
        form.pack(fill=tk.X, padx=8, pady=8)

        ttk.Label(form, text="E:").grid(row=0, column=0)
        var_e = tk.StringVar(value="121")
        ttk.Entry(form, textvariable=var_e, width=16).grid(row=0, column=1, padx=4)

        ttk.Label(form, text="base b:").grid(row=0, column=2, padx=(8, 0))
        var_b = tk.StringVar(value="10")
        ttk.Entry(form, textvariable=var_b, width=8).grid(row=0, column=3, padx=4)

        ttk.Label(form, text="Oráculo:").grid(row=0, column=4, padx=(8, 0))
        var_fam = tk.StringVar(value="cuadrados")
        ttk.Combobox(
            form,
            textvariable=var_fam,
            values=["general", "cuadrados", "cubos", "potencia", "kp", "mersenne", "none"],
            width=12,
            state="readonly",
        ).grid(row=0, column=5, padx=4)

        out = scrolledtext.ScrolledText(frame, height=18, font=("Consolas", 9))
        out.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        fig = Figure(figsize=(7, 2.5), dpi=100)
        ax = fig.add_subplot(111)
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        def run() -> None:
            try:
                E = float(var_e.get())
                b = float(var_b.get())
            except ValueError:
                messagebox.showerror("Error", "E o b inválidos")
                return
            fam = var_fam.get()
            if fam == "none":
                r = newton.newton_rapido(E, b=b)
            else:
                r = newton.log_con_oraculo(E, b=b, familia=fam)
            lines = [
                f"Newton Rápido — E={E}  b={b}",
                f"log_b(E) exacto : {r['j_exacto']:.14f}",
                f"log_b(E) VMA    : {r['j']:.14f}",
                f"iteraciones     : {r['iteraciones']}",
                f"error abs       : {r['error']:.3e}",
            ]
            if "j0" in r:
                lines.append(f"j0 oráculo      : {r['j0']:.8f}  ({r.get('familia', '')})")
            lines.append(f"\nHistorial j ({len(r['historial'])} pts):")
            for i, jv in enumerate(r["historial"][:20]):
                lines.append(f"  [{i:2d}] {jv:.12f}")
            if len(r["historial"]) > 20:
                lines.append("  ...")
            self._log(out, "\n".join(lines))

            ax.clear()
            hist = r["historial"]
            ax.plot(range(len(hist)), hist, "o-", color="#3498db", label="j iter")
            ax.axhline(r["j_exacto"], color="#e74c3c", ls="--", label="exacto")
            ax.set_xlabel("iteración")
            ax.set_ylabel("j")
            ax.set_title("Convergencia Newton Rápido")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
            fig.tight_layout()
            canvas.draw()

        ttk.Button(form, text="Calcular", command=run).grid(row=1, column=0, columnspan=6, pady=8)

    def _tab_classic(self, nb: ttk.Notebook) -> None:
        frame = ttk.Frame(nb)
        nb.add(frame, text="π(x) clásico")

        top = ttk.Frame(frame)
        top.pack(fill=tk.X, padx=8, pady=8)
        ttk.Label(top, text="x:").pack(side=tk.LEFT)
        var_x = tk.StringVar(value="10000")
        ttk.Entry(top, textvariable=var_x, width=12).pack(side=tk.LEFT, padx=4)

        out = scrolledtext.ScrolledText(frame, height=14, font=("Consolas", 9))
        out.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        def run() -> None:
            try:
                x = int(var_x.get())
            except ValueError:
                messagebox.showerror("Error", "x inválido")
                return
            t0 = time.perf_counter()
            n = classic.pi_count(x)
            ms = (time.perf_counter() - t0) * 1000
            self._log(
                out,
                f"π({x}) = {n}\n"
                f"densidad exacta  : {n / x:.8f}\n"
                f"PNT x/ln(x) ~    : {classic.pnt_estimate(x):.2f}\n"
                f"Criva π(x) ~     : {criva.criva(x) * x:.2f}\n"
                f"Eratóstenes      : {ms:.2f} ms",
            )

        ttk.Button(top, text="Contar primos", command=run).pack(side=tk.LEFT, padx=8)


def main() -> int:
    root = tk.Tk()
    try:
        ttk.Style().theme_use("vista")
    except tk.TclError:
        pass
    VmaMethodsApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())