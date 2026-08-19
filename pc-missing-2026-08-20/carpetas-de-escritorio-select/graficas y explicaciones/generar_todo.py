#!/usr/bin/env python3
"""
Genera gráficas VMA — información no convencional del ecosistema local.
Víctor Manzanares Alberola · 2026
"""
from __future__ import annotations

import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
ANTIPC = Path(r"C:\Users\cuent\Desktop\ideas-para-gpt-antipc\antipc\src\antipc")
sys.path.insert(0, str(ANTIPC))

from mdc_lib.mdc_analyze import analyze, diofantic_y, gcd_candidates  # noqa: E402
from mdc_lib.mdc_visual import animate_trains  # noqa: E402

OUT = ROOT
E_MINUS_2 = math.e - 2.0


def L(n: int) -> int:
    return int(math.sqrt(n + 3)) + 7


def m_sum(n: int) -> float:
    k = max(2, min(int(int(math.sqrt(n)) * 9 / 24), 2000))
    s = sum(1.0 / math.factorial(i) for i in range(2, k + 1))
    return math.sqrt(n + 3) * s


def calibrar_mrauv(n0: int) -> dict:
    dn = int(2 * math.sqrt(n0))
    pts = [n0, n0 + dn, n0 + 2 * dn]
    med = [L(x) - m_sum(x) for x in pts]
    med0, med1, med2 = med
    v0 = 1.0 / med0 if med0 else 0.0
    a0 = v0 / med1 if med1 else 0.0
    j = a0 / med2 if med2 else 0.0
    return {"n0": n0, "dn": dn, "med": med, "V0": v0, "a0": a0, "j": j}


def save_explain(name: str, text: str) -> None:
    (OUT / f"{name}_EXPLICACION.txt").write_text(text.strip() + "\n", encoding="utf-8")


def fig01_mdc_trenes_static(n: int, tag: str) -> None:
    r = analyze(n)
    proper = [(c.x, c.y, c.s, c.t) for c in r.train_x.collisions
              if isinstance(c.x, int) and isinstance(c.y, int) and c.s > 1 and c.t < n]
    x_lo, x_hi = -2, max(10, math.isqrt(n) + 6)
    xs = np.linspace(x_lo, x_hi, 400)
    ys = []
    xv = []
    for xf in xs:
        y = diofantic_y(n, Fraction(xf).limit_denominator(5000))
        if y is None or abs(4 * xf + 6) < 1e-6:
            continue
        xv.append(xf)
        ys.append(float(y))

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.plot(xv, ys, "#4a90d9", lw=2, label="(2x+3)(2y+3)=n")
    for xi, yi, s, t in proper:
        ax.scatter([xi], [yi], s=300, marker="*", c="#e74c3c", zorder=5, edgecolors="k")
        ax.annotate(f"{s}×{t}", (xi, yi), xytext=(8, 8), textcoords="offset points", fontsize=11)
    ax.axhline(0, color="#bbb", lw=0.5)
    ax.axvline(0, color="#bbb", lw=0.5)
    ax.set_title(f"MDC dos trenes — colisiones enteras n={n}")
    ax.set_xlabel("x  (tren X avanza Δx=1)")
    ax.set_ylabel("y  (tren Y avanza Δy=1)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / f"01_MDC_trenes_{tag}.png", dpi=150)
    plt.close(fig)

    save_explain(
        f"01_MDC_trenes_{tag}",
        f"""GRÁFICA 01 — MDC DOS TRENES (n={n})
═══════════════════════════════════════
Qué no enseña la academia:
  La factorización no se plantea como "probar divisores hasta √n",
  sino como DOS TRENES con velocidad constante en ejes distintos.

Ecuación VMA:  (2x+3)(2y+3) = {n}
  • Tren X: x entero, paso 1 → calcula y = (n-6x-9)/(4x+6)
  • Tren Y: y entero, paso 1 → calcula x simétrico
  • Colisión: x e y enteros simultáneos → S=2x+3, T=2y+3, S×T={n}

En esta imagen las estrellas rojas son colisiones CONFIRMADAS (unión X∩Y).
Factores propios: {', '.join(f'{s}×{t}' for _,_,s,t in proper) or 'ninguno en rango'}.

Comando local: antipc mdc analyze {n} --proper
""",
    )


def fig02_restos_hipotenusa(n: int) -> None:
    cap = min(n - 1, 150)
    i = np.arange(1, cap + 1)
    rest = n % i
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(i, rest, "#2c3e50", lw=1.5, label=f"{n} mod i")
    ax.plot([0, cap], [0, cap], "k--", lw=1.5, label="Hipotenusa (referencia)")
    for p in (1, 2, 3, 5):
        ax.plot(i, p * i, "--", alpha=0.35, label=f"pendiente {p}")
    ax.set_title(f"Restos modulares + hipotenusa — n={n}")
    ax.set_xlabel("i")
    ax.set_ylabel("resto")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / f"02_restos_hipotenusa_{n}.png", dpi=150)
    plt.close(fig)

    save_explain(
        f"02_restos_hipotenusa_{n}",
        f"""GRÁFICA 02 — RESTOS + HIPOTENUSA (n={n})
═══════════════════════════════════════
Método gráfico VMA (§9 cribas) — no está en libros estándar de teoría de números.

Se grafica (i, n mod i). La hipotenusa y=x del cuadrado de lado n es REFERENCIA VISUAL.
Las pendientes 1,2,3... marcan transiciones entre secuencias de restos.

Lectura "a ojo": los patrones de la curva de restos anticipan estructura modular
antes de aplicar MDC o criba 6k±1. No es prueba formal; es brújula geométrica.

Script: restos_graficos.py · antipc mdc visual {n}
""",
    )


def fig03_mrauv_tramos() -> None:
    n0_vals = list(range(100, 5000, 120))
    v0s, a0s, gaps = [], [], []
    for n0 in n0_vals:
        c = calibrar_mrauv(n0)
        v0s.append(c["V0"])
        a0s.append(c["a0"])
        gaps.append(L(n0) - m_sum(n0))

    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    axes[0].plot(n0_vals, gaps, "#27ae60", lw=1.5)
    axes[0].axhline(2, color="red", ls="--", label="umbral 2 primos")
    axes[0].set_ylabel("L(n) - m(n)")
    axes[0].set_title("MRAUV — criterio espacio libre para 2 primos")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].plot(n0_vals, v0s, "#8e44ad", lw=1.5)
    axes[1].set_ylabel("V₀ (velocidad)")
    axes[1].grid(alpha=0.3)

    axes[2].plot(n0_vals, a0s, "#d35400", lw=1.5)
    axes[2].set_ylabel("a₀ (aceleración)")
    axes[2].set_xlabel("n₀ inicio de tramo")
    axes[2].grid(alpha=0.3)

    fig.suptitle("MRAUV — densidad cinemática de primos por tramos", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "03_MRAUV_calibracion_tramos.png", dpi=150)
    plt.close(fig)

    save_explain(
        "03_MRAUV_calibracion_tramos",
        """GRÁFICA 03 — MRAUV CALIBRACIÓN POR TRAMOS
═══════════════════════════════════════
MRAUV = modelo VMA de densidad de primos (NO está en Hardy-Littlewood ni en Analytic NT estándar).

Funciones propias:
  L(n) = ⌊√(n+3)⌋ + 7
  m(n) ≈ (e-2)·√n   — "marcado" de compuestos inducidos

Criterio VMA: si L(n) - m(n) ≥ 2 → el intervalo I(n) contiene al menos 2 primos.

Calibración cinemática en 3 puntos (separación ~2√n₀):
  V₀ = velocidad de densidad
  a₀ = aceleración
  j  = jerk (tercera derivada)

Uso: Goldbach heurístico, cribas adaptativas, segmentación AntiPC.
Código: mrauv.py · archivos-vma/codigo/
""",
    )


def fig04_sawtooth_factor(n: int) -> None:
    ms = np.arange(1, min(120, math.isqrt(n)))
    ds = []
    for m in ms:
        denom = 2 * (2 * m + 3)
        val = n / denom
        ds.append(val - math.floor(val))
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(ms, ds, drawstyle="steps-post", color="#2980b9", lw=1)
    ax.axhline(0.5, color="#e74c3c", ls="--", label="δ=0.5 → divisor")
    hits = [m for m in ms if abs(ds[int(m - 1)] - 0.5) < 1e-9 or n % (2 * int(m) + 3) == 0]
    for m in ms:
        s = 2 * int(m) + 3
        if n % s == 0 and s < n:
            ax.scatter([m], [ds[int(m - 1)]], c="#e74c3c", s=80, zorder=5)
    ax.set_title(f"MDC sawtooth d(m) = frac(N/(2(2m+3))) — N={n}")
    ax.set_xlabel("m  (candidato S=2m+3)")
    ax.set_ylabel("parte decimal")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / f"04_MDC_sawtooth_{n}.png", dpi=150)
    plt.close(fig)

    save_explain(
        f"04_MDC_sawtooth_{n}",
        f"""GRÁFICA 04 — MDC SAWTOOTH (n={n})
═══════════════════════════════════════
El Método Diofántico Cinemático trata la factorización como FÍSICA de una onda dentada.

d(m) = frac({n} / (2(2m+3)))
Cuando d(m) = 0.5 exactamente → (2m+3) divide {n}.

La academia enseña trial division o Pollard; VMA enseña:
  medir velocidad/aceleración/jerk del sawtooth y SALTAR al m* correcto.

Implementado: vma-run/lib/mdc.py · mdc_predict()
""",
    )


def fig05_gcd_heatmap(n: int) -> None:
    ms = np.arange(0, min(200, math.isqrt(n) + 20))
    gc = [math.gcd(n, 2 * m + 3) for m in ms]
    fig, ax = plt.subplots(figsize=(10, 3))
    colors = ["#ecf0f1" if g == 1 else "#e74c3c" if n % (2 * m + 3) == 0 else "#f39c12" for m, g in zip(ms, gc)]
    ax.bar(ms, gc, color=colors, width=1.0, edgecolor="none")
    ax.set_title(f"gcd(n, 2m+3) — rejilla 6k±1 — n={n}")
    ax.set_xlabel("m")
    ax.set_ylabel("gcd")
    fig.tight_layout()
    fig.savefig(OUT / f"05_gcd_rejilla_{n}.png", dpi=150)
    plt.close(fig)

    save_explain(
        f"05_gcd_rejilla_{n}",
        f"""GRÁFICA 05 — REJILLA gcd(n, 2m+3) (n={n})
═══════════════════════════════════════
Los factores impares >3 viven en la clase 6k±1, parametrizados como S=2m+3.

Barras rojas: gcd>1 Y S divide n → factor encontrado.
Naranja: gcd>1 pero no divisor completo.
Gris: coprimo.

Esto es el "núcleo MDC" de mdc_hypotenuse_jump.py — criba modular propia VMA.
""",
    )


def fig06_antipc_reuse() -> None:
    """Simula curva conocimiento AntiPC (reuse vs ALU)."""
    ops = np.arange(1, 501)
    repeat_p = 0.35
    hits = np.cumsum(np.random.default_rng(42).random(500) < repeat_p)
    executed = ops - hits
    saved_pct = 100 * hits / ops

    fig, ax1 = plt.subplots(figsize=(9, 5))
    ax1.plot(ops, executed, label="ALU ejecutadas", color="#3498db")
    ax1.plot(ops, hits, label="Knowledge hits (reuse)", color="#2ecc71")
    ax1.set_xlabel("Operaciones K3_HASH")
    ax1.set_ylabel("Conteo acumulado")
    ax1.legend(loc="upper left")
    ax1.grid(alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(ops, saved_pct, "--", color="#9b59b6", label="% ALU ahorrada")
    ax2.set_ylabel("% ahorro")
    ax2.legend(loc="center right")

    ax1.set_title("AntiPC Flow Kernel — reuse de conocimiento (demo 35% repeat)")
    fig.tight_layout()
    fig.savefig(OUT / "06_AntiPC_knowledge_reuse.png", dpi=150)
    plt.close(fig)

    save_explain(
        "06_AntiPC_knowledge_reuse",
        """GRÁFICA 06 — ANTIPC KNOWLEDGE REUSE
═══════════════════════════════════════
AntiPC (Víctor Manzanares Alberola) NO es un hash tool ni un cache LRU estándar.

Ley VMA: "Mover conocimiento, no recalcular ALU."
  • Flow Kernel + Knowledge Buffer + plugins (K3, MDC)
  • Segunda ejecución idéntica → 0 ALU (reuse)

Esto no aparece en OS textbooks: es runtime de conocimiento con firma K3.
Demo real: antipc hash --text "Hola" → 2ª vez 0 ms.

Backend actual: k3hash.dll nativa (HASHTOOLCODE).
""",
    )


def fig07_ecosistema_mapa() -> None:
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    boxes = [
        (1, 7.5, "L0 K3\nk3hash.c", "#3498db"),
        (4, 7.5, "MDC\n(2x+3)(2y+3)", "#e74c3c"),
        (7, 7.5, "MRAUV\nL(n)-m(n)", "#9b59b6"),
        (1, 4.5, "AntiPC\nFlow Kernel", "#2ecc71"),
        (4, 4.5, "Aleatorovix\nRTL+MDC", "#f39c12"),
        (7, 4.5, "ZypyZape\nUDP gemelo", "#1abc9c"),
        (2.5, 1.5, "33×1 ecosistema VMA", "#34495e"),
    ]
    for x, y, txt, col in boxes:
        ax.add_patch(plt.Rectangle((x, y), 2.2, 1.2, fill=True, color=col, alpha=0.35, ec=col, lw=2))
        ax.text(x + 1.1, y + 0.6, txt, ha="center", va="center", fontsize=9, fontweight="bold")
    ax.annotate("", xy=(4, 8.1), xytext=(3.2, 8.1), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.annotate("", xy=(7, 8.1), xytext=(6.2, 8.1), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.annotate("", xy=(2.5, 5.7), xytext=(2.5, 7.5), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.set_title("Mapa ecosistema VMA — Desktop 2026", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "07_mapa_ecosistema_VMA.png", dpi=150)
    plt.close(fig)

    save_explain(
        "07_mapa_ecosistema_VMA",
        """GRÁFICA 07 — MAPA ECOSISTEMA VMA
═══════════════════════════════════════
Tu PC contiene un stack coherente que la academia no tiene como unidad:

  L0     K3 hash (HASHTOOLCODE) — huella, dedup
  MDC    Factorización diofántica cinemática + visual
  MRAUV  Densidad de primos por tramos (V₀,a₀,j)
  AntiPC Runtime de conocimiento (Python + .exe + DLL)
  Aleatorovix  Azar masivo RTL + criba + MDC memoria
  ZypyZape     Gemelo UDP / piloto industrial

Raíz: C:\\Users\\cuent\\Desktop\\
""",
    )


def fig08_union_xy(n: int) -> None:
    r = analyze(n)
    labels = []
    xb, yb = [], []
    for c in r.train_x.collisions:
        if isinstance(c.x, int) and c.s > 1 and c.t < n:
            labels.append(f"X:{c.x}")
            xb.append(1)
            yb.append(c.s)
    for c in r.train_y.collisions:
        if isinstance(c.y, int) and c.s > 1 and c.t < n:
            labels.append(f"Y:{c.y}")
            xb.append(2)
            yb.append(c.s)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(xb, yb, s=120, c="#e74c3c")
    for i, lb in enumerate(labels):
        ax.annotate(lb, (xb[i], yb[i]), fontsize=9)
    for s, t in r.union_both:
        if s > 1 and t < n:
            ax.axhline(s, color="#2ecc71", ls=":", alpha=0.7)
            ax.text(1.5, s, f"∩ {s}×{t}", color="#27ae60")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Tren X", "Tren Y"])
    ax.set_ylabel("Factor S=2m+3")
    ax.set_title(f"Unión resultados iguales X∩Y — n={n}")
    fig.tight_layout()
    fig.savefig(OUT / f"08_union_XY_{n}.png", dpi=150)
    plt.close(fig)

    save_explain(
        f"08_union_XY_{n}",
        f"""GRÁFICA 08 — UNIÓN X∩Y (n={n})
═══════════════════════════════════════
Dos trenes independientes (variación x y variación y) deben dar el MISMO par (S,T).

Si coincide → factorización confirmada por simetría del método.
Pares en intersección: {r.union_both}

Comando: antipc mdc analyze {n} --proper
""",
    )


def fig13_zypyzape_ab() -> None:
    """Comparativa medida A (memcpy+cola) vs B (slot ring) desde viability JSON."""
    candidates = [
        ROOT / "zypyzape" / "zypyzape_viability.json",
        Path(r"C:\Users\cuent\Desktop\ideas-para-gpt-antipc\antipc\output\zypyzape_viability.json"),
    ]
    data = None
    for path in candidates:
        if path.is_file():
            data = json.loads(path.read_text(encoding="utf-8"))
            break
    if not data or len(data.get("runs", [])) < 2:
        print("  [skip] zypyzape_viability.json no encontrado — ejecuta scripts\\19_zypyzape_viabilidad.bat")
        return

    runs = data["runs"]
    labels = ["A\nconvencional", "B\nslot ring"]
    throughputs = [
        r["packets_processed"] / max(r["elapsed_s"], 1e-9) for r in runs[:2]
    ]
    copies = [r["memcpy_user_copies"] for r in runs[:2]]
    hubs = data.get("hubs", "?")
    ts = data.get("timestamp", "")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    colors = ["#e74c3c", "#2ecc71"]
    ax1.bar(labels, throughputs, color=colors, edgecolor="k", linewidth=0.8)
    ax1.set_ylabel("Paquetes/s")
    ax1.set_title(f"Throughput UDP localhost ({hubs} hubs)")
    ax1.grid(axis="y", alpha=0.3)
    for i, v in enumerate(throughputs):
        ax1.text(i, v, f"{v:.0f}", ha="center", va="bottom", fontsize=10)

    ax2.bar(labels, copies, color=colors, edgecolor="k", linewidth=0.8)
    ax2.set_ylabel("Copias usuario (recv + memcpy)")
    ax2.set_title("Presión memcpy en espacio usuario")
    ax2.grid(axis="y", alpha=0.3)
    for i, v in enumerate(copies):
        ax2.text(i, v, f"{v}", ha="center", va="bottom", fontsize=10)

    if copies[0] > 0:
        save_pct = (1 - copies[1] / copies[0]) * 100
        gain = (throughputs[1] - throughputs[0]) / throughputs[0] * 100 if throughputs[0] else 0
        fig.suptitle(
            f"ZypyZape viabilidad — B ahorra {save_pct:.0f}% copias, throughput {gain:+.1f}%",
            fontsize=12,
            y=1.02,
        )

    fig.tight_layout()
    fig.savefig(OUT / "13_ZypyZape_A_vs_B.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    save_explain(
        "13_ZypyZape_A_vs_B",
        f"""GRÁFICA 13 — ZypyZape A vs B (medición real)
═══════════════════════════════════════
Benchmark UDP en localhost (HP Envy, Windows).
Timestamp JSON: {ts}

  A) recvinto(temp) + bytes(temp) → Queue  → 2 copias usuario/pkt
  B) recvinto(slot fijo del ring)          → 1 copia usuario/pkt

Datos:
  Throughput A: {throughputs[0]:.0f} pkt/s   copias: {copies[0]}
  Throughput B: {throughputs[1]:.0f} pkt/s   copias: {copies[1]}

B NO es DMA del NIC: evita la segunda memcpy usuario entre red y worker.
Pipeline completo AntiPC: B → Grafcet (D) con `antipc network demo`.

Regenerar medición:
  ideas-para-gpt-antipc\\antipc\\scripts\\19_zypyzape_viabilidad.bat
""",
    )


def write_index() -> None:
    files = sorted(OUT.glob("*.png")) + sorted(OUT.glob("*.gif"))
    lines = [
        "ÍNDICE — GRÁFICAS Y EXPLICACIONES VMA",
        "Víctor Manzanares Alberola · Desktop 2026",
        "=" * 50,
        "",
        "Contenido NO convencional: MDC, MRAUV, AntiPC, rejilla 6k±1,",
        "dos trenes diofánticos, sawtooth, restos+hipotenusa.",
        "",
    ]
    for f in files:
        exp = OUT / f"{f.stem}_EXPLICACION.txt"
        lines.append(f"  {f.name}")
        if exp.exists():
            lines.append(f"    → {exp.name}")
    lines += [
        "",
        "REGENERAR:",
        "  python generar_todo.py",
        "",
        "COMANDOS:",
        "  antipc mdc analyze 10403 --proper",
        "  antipc mdc visual 10403 --proper -o trenes.gif",
        "  antipc network demo --hubs 4",
        "  scripts\\19_zypyzape_viabilidad.bat  (regenera JSON A vs B)",
    ]
    (OUT / "00_INDICE.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    plt.rcParams["figure.facecolor"] = "white"
    print(f"Salida: {OUT}")

    for n, tag in [(35, "35"), (10403, "10403"), (143, "143")]:
        fig01_mdc_trenes_static(n, tag)
        fig02_restos_hipotenusa(n)
        fig05_gcd_heatmap(n)
        fig08_union_xy(n)

    fig03_mrauv_tramos()
    fig04_sawtooth_factor(10403)
    fig04_sawtooth_factor(35)
    fig06_antipc_reuse()
    fig07_ecosistema_mapa()
    fig13_zypyzape_ab()

    print("GIF trenes 10403...")
    animate_trains(10403, proper_only=True, save_path=OUT / "09_MDC_animacion_10403.gif", show=False, fps=6)
    print("GIF trenes 35...")
    animate_trains(35, proper_only=True, save_path=OUT / "10_MDC_animacion_35.gif", show=False, fps=6)

    write_index()
    print("OK — carpeta completa")


if __name__ == "__main__":
    main()