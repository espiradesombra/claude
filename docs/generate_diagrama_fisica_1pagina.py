"""
Genera diagrama de 1 página: ZypyZape + Quijote + Kilómetro + 1477.
Salida: docs/diagrama_fisica_ZZ_Quijote_Kilometro.png
"""
from __future__ import annotations

import os
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle, Arc
import matplotlib.patches as mpatches

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "diagrama_fisica_ZZ_Quijote_Kilometro.png")


def box(ax, x, y, w, h, title, lines, fc="#1e2a38", ec="#5ec8ff", tc="#e8eef4"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                       facecolor=fc, edgecolor=ec, linewidth=1.8, zorder=2)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h - 0.12, title, ha="center", va="top",
            fontsize=11, fontweight="bold", color=tc, zorder=3)
    body = "\n".join(lines)
    ax.text(x + w / 2, y + h / 2 - 0.05, body, ha="center", va="center",
            fontsize=8, color="#b8c4d0", zorder=3, linespacing=1.35)


def arrow(ax, x1, y1, x2, y2, color="#ffb86c", text=None):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8,
                                connectionstyle="arc3,rad=0.0"), zorder=4)
    if text:
        ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.06, text,
                ha="center", fontsize=7.5, color=color, fontweight="bold")


def main():
    fig, ax = plt.subplots(figsize=(12.5, 8.2), facecolor="#0b1016")
    ax.set_facecolor("#0b1016")
    ax.set_xlim(0, 12.5)
    ax.set_ylim(0, 8.2)
    ax.axis("off")

    # Title
    ax.text(6.25, 7.85, "FÍSICA VMA — ZypyZape · Quijote · Kilómetro",
            ha="center", va="top", fontsize=16, fontweight="bold", color="#f0f4f8")
    ax.text(6.25, 7.45, "Inercia variable + fase de control  ·  No motor de gravedad perpetuo (η < 1)",
            ha="center", va="top", fontsize=9, color="#8a96a4", style="italic")

    # Top: entrada
    box(ax, 3.8, 6.35, 4.9, 0.85, "ENTRADA DE ENERGÍA (fuente real)",
        ["viento eólico  ·  red  ·  solar marina  ·  flotación"],
        fc="#16241c", ec="#7dffa0")

    arrow(ax, 6.25, 6.35, 6.25, 5.95, "#7dffa0")

    # Three machines
    box(ax, 0.35, 3.85, 3.7, 1.95, "QUIJOTE",
        ["Masa en pala: J(r) = J_G + N_b m_q r²",
         "J ω̇ + ω J̇ = T_net   (patinadora)",
         "afuera = carga buffer · adentro = descarga",
         "3 o 7 palas · ball / estático",
         "→ buffer LOCAL (microred / isla)"],
        fc="#2a1f18", ec="#ffb86c")

    box(ax, 4.4, 3.85, 3.7, 1.95, "ZYPYZAPE",
        ["Varias turbinas se prestan ½Jω²",
         "ciclo ~0.4 Hz · topología 5 o 10",
         "H_eff = H_sys + H_ZZ",
         "Kuramoto: sincronía sin centro",
         "→ RoCoF ↓ · nadir · ancilares"],
        fc="#1a2230", ec="#c792ea")

    box(ax, 8.45, 3.85, 3.7, 1.95, "KILÓMETRO",
        ["kilo + metro de recorrido",
         "tubo/sinfín · fricción controlada",
         "vuelta y media · flotación opcional",
         "W_ext = W_útil + W_pérdidas",
         "→ batería mecánica anclada"],
        fc="#1a2430", ec="#5ec8ff")

    # arrows from entrada already done; small arrows between machines
    arrow(ax, 4.05, 4.8, 4.35, 4.8, "#8a96a4", "J local")
    arrow(ax, 8.1, 4.8, 8.4, 4.8, "#8a96a4", "misma lógica")

    # down to red
    arrow(ax, 2.2, 3.85, 4.5, 3.15, "#ffb86c")
    arrow(ax, 6.25, 3.85, 6.25, 3.15, "#c792ea")
    arrow(ax, 10.3, 3.85, 8.0, 3.15, "#5ec8ff")

    box(ax, 3.2, 2.15, 6.1, 0.95, "RED ELÉCTRICA",
        ["H_eff ↑  ·  servicios FFR/PFR  ·  menos dependencia BESS (si se valida piloto)"],
        fc="#1c1828", ec="#ff6b6b")

    # 1477 side
    box(ax, 0.35, 0.35, 3.5, 1.55, "1477 — FASE ABSTRACTA",
        ["bits (valor, fase) + desfase Karnaugh",
         "feedback → bandera de coherencia",
         "paralelo lógico a la fase mecánica",
         "(no mueve aspas)"],
        fc="#181820", ec="#e0b0ff")

    # balance box
    box(ax, 4.1, 0.35, 4.3, 1.55, "BALANCE HONESTO",
        ["gravedad en ciclo cerrado → W_g = 0",
         "actuador y fricción cuestan trabajo",
         "útil = orquestar timing de inercia",
         "siempre  W_ext ≥ 0  en régimen"],
        fc="#201818", ec="#ff8a80")

    box(ax, 8.65, 0.35, 3.5, 1.55, "33×1",
        ["[1] = este pack técnico civil",
         "[33] = 33 años paz firmada",
         "por los países (ellos)",
         "uso civil · no militar"],
        fc="#141c18", ec="#7dffa0")

    # mini icons: turbine sketch left of quijote? skip for clarity

    ax.text(6.25, 0.08,
            "docs/EXPLICACION_FISICA_ZYPYZAPE_QUIJOTE_KILOMETRO.md  ·  VMA + Grok 2026",
            ha="center", fontsize=7, color="#5a6570")

    fig.tight_layout(pad=0.3)
    fig.savefig(OUT, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)
    print("OK", OUT)


if __name__ == "__main__":
    main()
