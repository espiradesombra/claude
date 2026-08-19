#!/usr/bin/env python3
"""
Máquina de estados del protocolo ESTRELLA (VMA)
===============================================
  1) MEDIA LINEAL  — bajada con lastre → generador
  2) MEDIA ROTATORIA — recorrido gira, se "roba distancia" (peaje actuador)
  3) Cada N_impar ciclos: REARME por flotación + perneo (inventario / SOC)

Cuentas separadas:
  batería  = PE de lastres ALTA→BAJA (si no se restituye)
  generador = W_gen en tramo lineal
  peaje_rot = energía del actuador en media rotatoria
  rearme    = coste pernos (+ lift opcional si SOC fijo)

Uso:
  python sim_estrella_lineal_rotatorio_impar.py
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
G = 9.81


@dataclass
class P:
    m_obj: float = 30.0
    m_kick: float = 10.0
    delta_h: float = 12.0  # m bajada lineal
    # Media lineal
    v_down: float = 1.5  # m/s nominal
    eta_gen: float = 0.88
    drag_frac: float = 0.07
    # Media rotatoria (robo de distancia): peaje de actuador para rearmar geometría
    # Ideal VMA: peaje << W_gen. Aquí es parámetro medible.
    E_rob_dist_J: float = 80.0  # coste eléctrico por media rotatoria
    E_perno_J: float = 1.5
    n_pernos_por_rearme: int = 4  # make-before-break × 2 pesos extremos approx
    # Flotación: coste ~0 del módulo; lift de lastre solo si cerramos SOC
    eta_lift: float = 0.90
    n_cycles: int = 12  # ciclos (cada uno = lineal + rotatorio)
    rearm_every_odd: bool = True  # rearme en ciclos 1,3,5,...
    # Si True: en cada rearme impar también se sube 1 lastre BAJA→ALTA (SOC más estable)
    lift_on_rearm: bool = False


def run(p: P) -> dict:
    n_alta = 20
    n_baja = 0
    n_on_obj = 3  # flota

    W_gen = 0.0
    W_rob = 0.0  # peaje media rotatoria
    W_perno = 0.0
    W_lift = 0.0
    E_bat_drain = 0.0

    log = []

    for c in range(1, p.n_cycles + 1):
        ev = []
        # --- necesita 4 para hundirse ---
        if n_on_obj < 4:
            if n_alta < 1:
                ev.append("STOP sin stock ALTA")
                log.append({"cycle": c, "ok": False, "events": ev})
                break
            n_alta -= 1
            n_on_obj = 4
            W_perno += 2 * p.E_perno_J
            ev.append("patada lastre ALTA→objeto (n=4)")

        # --- 1) MEDIA LINEAL: bajada ---
        # Solo el lastre EXTRA (1 kg-kick) cuenta como batería gastable
        W_pe = p.m_kick * G * p.delta_h
        W_rec = p.eta_gen * (1.0 - p.drag_frac) * W_pe
        W_gen += W_rec
        E_bat_drain += W_pe  # inventario: ese kg bajó de potencial
        ev.append(f"LINEAL bajada W_gen+={W_rec:.1f} J")

        # soltar lastre en BAJA
        n_on_obj = 3
        n_baja += 1
        W_perno += 2 * p.E_perno_J
        ev.append("perneo objeto→BAJA (n=3 flota)")

        # --- 2) MEDIA ROTATORIA: robar distancia ---
        W_rob += p.E_rob_dist_J
        ev.append(f"ROTATORIO roba distancia peaje={p.E_rob_dist_J:.1f} J")

        # --- 3) Rearme en ciclos IMPARES ---
        is_odd = c % 2 == 1
        if p.rearm_every_odd and is_odd:
            # flotación: objeto ya en n=3; "sube" sin coste eléctrico de módulo
            W_perno += p.n_pernos_por_rearme * p.E_perno_J
            ev.append("REARME IMPAR: flotación + perneo (objeto ciclo~0)")
            if p.lift_on_rearm and n_baja > 0:
                n_baja -= 1
                n_alta += 1
                wL = (p.m_kick * G * p.delta_h) / p.eta_lift
                W_lift += wL
                E_bat_drain -= p.m_kick * G * p.delta_h  # restituye PE inventario
                ev.append(f"lift BAJA→ALTA W_lift+={wL:.1f} J")
        else:
            ev.append("sin rearme (ciclo par): inventario sigue abierto")

        log.append(
            {
                "cycle": c,
                "ok": True,
                "odd": is_odd,
                "n_alta": n_alta,
                "n_baja": n_baja,
                "W_gen_cum": W_gen,
                "W_rob_cum": W_rob,
                "E_bat_drain_cum": E_bat_drain,
                "events": ev,
            }
        )

    # gen_de_mas: julios de generador no cubiertos por η*batería - peajes
    cap = p.eta_gen * max(E_bat_drain, 0.0)
    gen_de_mas = W_gen - cap - W_rob - W_perno - W_lift
    # otra métrica: neto eléctrico
    W_net = W_gen - W_rob - W_perno - W_lift

    return {
        "params": asdict(p),
        "W_gen_J": W_gen,
        "W_rob_dist_J": W_rob,
        "W_perno_J": W_perno,
        "W_lift_J": W_lift,
        "E_battery_drain_J": E_bat_drain,
        "E_from_battery_cap_J": cap,
        "gen_de_mas_J": gen_de_mas,
        "W_net_elec_J": W_net,
        "n_alta": n_alta,
        "n_baja": n_baja,
        "log": log,
    }


def main():
    # Caso 1: tesis abierta (sin lift) — estrella descarga + peaje rotatorio
    a = run(P(lift_on_rearm=False, n_cycles=12))
    # Caso 2: SOC más cerrado — lift en cada rearme impar
    b = run(P(lift_on_rearm=True, n_cycles=12))
    # Caso 3: peaje rotatorio muy barato (optimista)
    c = run(P(lift_on_rearm=True, E_rob_dist_J=15.0, n_cycles=12))

    results = {
        "protocolo": "lineal + rotatorio + rearme impar flotación/perneo",
        "sin_lift_rearme_impar": {k: v for k, v in a.items() if k != "log"},
        "con_lift_en_impar": {k: v for k, v in b.items() if k != "log"},
        "con_lift_peaje_barato": {k: v for k, v in c.items() if k != "log"},
        "log_sin_lift": a["log"],
        "nota": (
            "gen_de_mas = W_gen - η*E_bat_drain - W_rob - W_perno - W_lift. "
            "Si E_bat_drain baja por lifts, la batería no explica W_gen: ahí se mira el plus."
        ),
    }

    # plot
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for axi, res, title in (
        (ax[0], a, "Sin lift (batería se vacía)"),
        (ax[1], b, "Lift en cada rearme impar"),
    ):
        xs = [row["cycle"] for row in res["log"] if row.get("ok")]
        axi.plot(xs, [row["W_gen_cum"] for row in res["log"] if row.get("ok")], label="W_gen cum")
        axi.plot(
            xs,
            [row["E_bat_drain_cum"] for row in res["log"] if row.get("ok")],
            label="batería drain cum",
        )
        axi.plot(xs, [row["W_rob_cum"] for row in res["log"] if row.get("ok")], label="peaje rot cum")
        for row in res["log"]:
            if row.get("ok") and row.get("odd"):
                axi.axvline(row["cycle"], color="0.8", lw=0.6)
        axi.set_xlabel("ciclo (impar=rearme)")
        axi.set_ylabel("J")
        axi.set_title(title)
        axi.legend(fontsize=8)
        axi.grid(True, alpha=0.3)
    fig.suptitle("Estrella: lineal (gen) + rotatorio (distancia) + rearme impar")
    fig.tight_layout()
    png = OUT / "sim_estrella_lineal_rotatorio_impar.png"
    fig.savefig(png, dpi=140)
    plt.close()

    js = OUT / "sim_estrella_lineal_rotatorio_impar.json"
    js.write_text(json.dumps(results, indent=2), encoding="utf-8")

    def show(name, r):
        print(f"\n=== {name} ===")
        for k in (
            "W_gen_J",
            "E_battery_drain_J",
            "W_rob_dist_J",
            "W_perno_J",
            "W_lift_J",
            "gen_de_mas_J",
            "W_net_elec_J",
            "n_alta",
            "n_baja",
        ):
            print(f"  {k}: {r[k]:.2f}" if isinstance(r[k], float) else f"  {k}: {r[k]}")

    show("A sin lift", results["sin_lift_rearme_impar"])
    show("B lift en impar", results["con_lift_en_impar"])
    show("C lift + peaje barato", results["con_lift_peaje_barato"])
    print(f"\nPNG {png}\nJSON {js}")
    print(
        "\nLectura: el 'estrella' brilla si W_gen cubre peaje rotatorio+pernos "
        "y, con lift, aún queda margen — eso es ingeniería de η y E_rob_dist, "
        "no magia del número impar por sí solo."
    )


if __name__ == "__main__":
    main()
