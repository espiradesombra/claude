#!/usr/bin/env python3
"""
Experimento: recorrido que ROTA + patadas de lastre PAR vs IMPAR
================================================================
Protocolo pedido (VMA): no confundir con "stock 5 → 5 ciclos".

  - El recorrido (tubo) gira: ángulo φ.
  - Hay N_patadas de lastre por vuelta (N par o impar).
  - Dos cuentas SEPARADAS:
      BATERÍA  = ΔPE del inventario de lastres (kg que cambian de cota)
      GENERADOR = trabajo eléctrico regenerativo del eje (η_gen)

  - Modo A "almacén abierto": las patadas pueden dejar lastre abajo
    → batería se descarga; el generador cobra parte de esa PE.
  - Modo B "SOC fijo": tras cada vuelta se RESTITUYE el inventario
    (lift contado en E_in). Aquí se ve si queda generador "de más"
    cuando la batería no se vacía.

Salidas: JSON + PNG en este directorio.
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
    m_obj: float = 40.0  # kg objeto sin lastre de patada
    m_kick: float = 10.0  # kg por patada de lastre
    R: float = 1.2  # m  brazo efectivo grav→par sobre el eje del tubo
    I0: float = 25.0  # kg·m² inercia del tubo + sinfín
    omega_ref: float = 2.0  # rad/s
    beta: float = 1.2  # fricción viscosa N·m·s/rad
    c_gen: float = 70.0  # freno regenerativo en fase favorable
    tau_mot_max: float = 500.0
    eta_gen: float = 0.90
    eta_mot: float = 0.90
    eta_lift: float = 0.90
    # Caída de cota efectiva cuando una patada "aparca" lastre abajo (modo A)
    delta_h_ballast: float = 8.0  # m
    dt: float = 5e-4
    n_turns: int = 6  # vueltas del recorrido a integrar


def kick_phases(n_kicks: int) -> np.ndarray:
    """Fases de ENGANCHE en [0, 2π). Soltar = engache + π/n (mitad del sector)."""
    if n_kicks <= 0:
        return np.array([])
    return (2.0 * np.pi) * np.arange(n_kicks) / n_kicks


def engaged_at(phi: float, n_kicks: int) -> bool:
    """
    Lastre enganchado en la primera mitad de cada sector 2π/N.
    Con N impar los sectores no se emparejan simétricamente con
    favorable/desfavorable de cos(φ) en el mismo patrón que con N par.
    """
    if n_kicks <= 0:
        return False
    s = phi % (2.0 * np.pi)
    sector = 2.0 * np.pi / n_kicks
    local = s % sector
    return local < 0.5 * sector


def run(
    n_kicks: int,
    *,
    soc_fixed: bool,
    p: P | None = None,
) -> dict:
    p = p or P()
    I_base = p.I0 + p.m_obj * p.R**2

    phi = 0.0
    omega = p.omega_ref
    phi_end = p.n_turns * 2.0 * np.pi

    W_mot = 0.0
    W_gen = 0.0
    W_grav = 0.0
    W_fric = 0.0
    W_lift = 0.0

    # Inventario: cada patada completada (engage→release) puede mover 1 kg ALTA→BAJA
    n_alta = 100  # stock grande para no cortar el ensayo
    n_baja = 0
    kicks_completed = 0
    # Seguimiento de si estamos dentro de un engage abierto en este sector
    was_engaged = False
    pending_ballast_drop = False  # release en modo A → 1 lastre a BAJA

    # Series
    t_hist = []
    phi_hist = []
    eng_hist = []
    gen_pow = []

    t = 0.0
    max_steps = int(1e7)
    stalled = False

    for _ in range(max_steps):
        eng = engaged_at(phi, n_kicks)
        m = p.m_obj + (p.m_kick if eng else 0.0)
        I = p.I0 + m * p.R**2

        # Detección borde engage→release (fin de patada)
        if was_engaged and not eng:
            kicks_completed += 1
            if not soc_fixed:
                # Modo batería: aparca lastre abajo
                if n_alta > 0:
                    n_alta -= 1
                    n_baja += 1
                    pending_ballast_drop = True
            # Modo SOC fijo: no cambia inventario aquí; el lift se paga al cerrar vuelta
        was_engaged = eng

        # Par gravitatorio (conservativo sobre el eje del tubo)
        # Convención: PE = m g R sin(φ); tau_g = -m g R cos(φ)
        tau_g = -m * G * p.R * np.cos(phi)
        favorable = np.cos(phi) < 0.0

        if favorable:
            tau_c = -p.c_gen * max(omega, 0.0)
            mode = "gen"
        else:
            err = p.omega_ref - omega
            tau_c = min(p.tau_mot_max, max(0.0, 100.0 * err + 0.35 * p.tau_mot_max))
            mode = "mot"

        tau_f = -p.beta * omega
        alpha = (tau_g + tau_c + tau_f) / I

        W_grav += tau_g * omega * p.dt
        W_fric += abs(tau_f * omega) * p.dt
        if mode == "mot":
            W_mot += (max(tau_c * omega, 0.0) / p.eta_mot) * p.dt
        else:
            W_gen += (max(-tau_c * omega, 0.0) * p.eta_gen) * p.dt

        omega_new = omega + alpha * p.dt
        if omega_new <= 0.0:
            stalled = True
            omega = 0.0
            break

        phi_next = phi + omega * p.dt

        # Al cruzar cada vuelta completa: si SOC fijo, restituir lastres bajados
        # (en modo A no hay que restituir; en B no bajamos inventario en release,
        #  pero igual forzamos balance de inventario=constante: W_lift de las
        #  patadas de la vuelta * m g Δh / η, como peaje de mantener SOC).
        turn0 = int(phi / (2.0 * np.pi))
        turn1 = int(phi_next / (2.0 * np.pi))
        if soc_fixed and turn1 > turn0:
            # N_patadas por vuelta → N restituciones teóricas de lastre
            # (mantenemos stock constante; peaje lift explícito)
            for _k in range(n_kicks):
                W_lift += (p.m_kick * G * p.delta_h_ballast) / p.eta_lift

        phi = phi_next
        omega = omega_new
        t += p.dt

        if _ % 40 == 0:
            t_hist.append(t)
            phi_hist.append(phi)
            eng_hist.append(1.0 if eng else 0.0)
            gen_pow.append(max(-tau_c * omega, 0.0) * p.eta_gen if mode == "gen" else 0.0)

        if phi >= phi_end:
            break

    # Contabilidad batería: PE de inventario movida ALTA→BAJA
    n_moved = n_baja  # en modo A; en SOC fijo queda 0
    E_battery_out = n_moved * p.m_kick * G * p.delta_h_ballast  # PE liberada del almacén
    # Lo que el generador "debería" como máximo de esa batería (ideal)
    E_from_battery_cap = p.eta_gen * E_battery_out

    # "Generador de más" (heurística del protocolo):
    #  - Modo A: W_gen - E_from_battery_cap  (si >0 parecería plus; suele ≤0 por fricción/motor)
    #  - Modo B: W_gen - W_mot - W_lift      (neto eléctrico a SOC fijo)
    if soc_fixed:
        gen_de_mas = W_gen - W_mot - W_lift
        battery_drain = 0.0
    else:
        gen_de_mas = W_gen - E_from_battery_cap - W_mot
        battery_drain = E_battery_out

    return {
        "n_kicks": n_kicks,
        "parity": "impar" if (n_kicks % 2 == 1) else "par",
        "soc_fixed": soc_fixed,
        "stalled": stalled,
        "turns": float(phi / (2.0 * np.pi)),
        "kicks_completed": kicks_completed,
        "n_alta": n_alta,
        "n_baja": n_baja,
        "W_gen_J": W_gen,
        "W_mot_J": W_mot,
        "W_lift_J": W_lift,
        "W_grav_J": W_grav,
        "W_fric_abs_J": W_fric,
        "E_battery_drain_J": battery_drain,
        "E_from_battery_cap_J": E_from_battery_cap,
        "gen_de_mas_J": gen_de_mas,
        "W_net_elec_J": W_gen - W_mot - W_lift,
        "t_hist": t_hist,
        "phi_hist": phi_hist,
        "eng_hist": eng_hist,
        "params": asdict(p),
    }


def main():
    kick_list = [1, 2, 3, 4, 5, 6]
    rows_open = []
    rows_soc = []

    for n in kick_list:
        r_a = run(n, soc_fixed=False)
        r_b = run(n, soc_fixed=True)
        # no guardar series enormes en tabla
        for r, bag in ((r_a, rows_open), (r_b, rows_soc)):
            bag.append({k: v for k, v in r.items() if not k.endswith("_hist")})

    # --- plot ---
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    def plot_side(ax, rows, title):
        ns = [r["n_kicks"] for r in rows]
        gen = [r["W_gen_J"] for r in rows]
        bat = [r["E_battery_drain_J"] for r in rows]
        plus = [r["gen_de_mas_J"] for r in rows]
        net = [r["W_net_elec_J"] for r in rows]
        x = np.arange(len(ns))
        w = 0.2
        ax.bar(x - 1.5 * w, gen, w, label="W_gen")
        ax.bar(x - 0.5 * w, bat, w, label="batería (drain PE)")
        ax.bar(x + 0.5 * w, plus, w, label="gen_de_mas")
        ax.bar(x + 1.5 * w, net, w, label="W_net_elec")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{n}\n({('impar' if n%2 else 'par')})" for n in ns])
        ax.set_title(title)
        ax.set_ylabel("J (en el ensayo)")
        ax.legend(fontsize=8)
        ax.grid(True, axis="y", alpha=0.3)

    plot_side(axes[0, 0], rows_open, "Modo A — almacén abierto (batería puede vaciarse)")
    plot_side(axes[0, 1], rows_soc, "Modo B — SOC fijo (lift paga restitución)")

    # detalle temporal N=3 vs N=4 modo B
    for ax, n, color in ((axes[1, 0], 3, "C0"), (axes[1, 1], 4, "C1")):
        r = run(n, soc_fixed=True)
        ph = np.array(r["phi_hist"])
        eng = np.array(r["eng_hist"])
        ax.plot(ph / (2 * np.pi), eng, color=color, lw=0.8)
        ax.set_xlabel("vueltas del recorrido")
        ax.set_ylabel("lastre enganchado (0/1)")
        ax.set_title(f"Patrón de patadas N={n} ({r['parity']}) — SOC fijo")
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Recorrido que rota + patadas par/impar\n"
        "Cuentas separadas: batería (inventario) vs generador; gen_de_mas = tesis VMA",
        fontsize=11,
    )
    fig.tight_layout()
    png = OUT / "sim_recorrido_rota_patadas_impar.png"
    fig.savefig(png, dpi=140)
    plt.close()

    # Series cortas de ejemplo no en JSON completo
    summary = {
        "idea": (
            "N patadas/vuelta sobre recorrido que rota. "
            "Modo A: release deja lastre abajo (batería). "
            "Modo B: SOC fijo con lift explícito. "
            "gen_de_mas_J: A→ W_gen - η*PE_drain - W_mot; "
            "B→ W_gen - W_mot - W_lift."
        ),
        "modo_A_almacen_abierto": rows_open,
        "modo_B_SOC_fijo": rows_soc,
        "lectura": {
            "si_gen_de_mas_B_positivo": "habría generador de más a SOC fijo (tesis)",
            "si_gen_de_mas_B_negativo": "a SOC fijo el peaje (motor+lift+η) gana; no hay plus",
            "si_solo_A_positivo_por_batería": "julios del vaciado de inventario, no plus de fase impar",
        },
    }
    js = OUT / "sim_recorrido_rota_patadas_impar.json"
    js.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("=== Modo A (almacén abierto) ===")
    print(f"{'N':>3} {'paridad':>7} {'W_gen':>10} {'batería':>10} {'gen_de_mas':>12} {'W_net':>10}")
    for r in rows_open:
        print(
            f"{r['n_kicks']:3d} {r['parity']:>7} {r['W_gen_J']:10.1f} "
            f"{r['E_battery_drain_J']:10.1f} {r['gen_de_mas_J']:12.1f} {r['W_net_elec_J']:10.1f}"
        )
    print("\n=== Modo B (SOC fijo) ===")
    print(f"{'N':>3} {'paridad':>7} {'W_gen':>10} {'W_lift':>10} {'gen_de_mas':>12} {'W_net':>10}")
    for r in rows_soc:
        print(
            f"{r['n_kicks']:3d} {r['parity']:>7} {r['W_gen_J']:10.1f} "
            f"{r['W_lift_J']:10.1f} {r['gen_de_mas_J']:12.1f} {r['W_net_elec_J']:10.1f}"
        )
    print(f"\nPNG: {png}")
    print(f"JSON: {js}")

    # Veredicto automático honesto
    best_b = max(rows_soc, key=lambda r: r["gen_de_mas_J"])
    print("\n--- Lectura automática (modo B = prueba de la tesis) ---")
    print(
        f"Mejor gen_de_mas a SOC fijo: N={best_b['n_kicks']} ({best_b['parity']}) "
        f"→ {best_b['gen_de_mas_J']:.1f} J"
    )
    if best_b["gen_de_mas_J"] > 0:
        print("RESULTADO: aparece plus a SOC fijo en este modelo (revisar supuestos).")
    else:
        print(
            "RESULTADO: en este modelo 1.ª-ley + lift, el impar NO saca generador "
            "de más a SOC fijo. El modo A puede lucir julios por vaciar batería."
        )


if __name__ == "__main__":
    main()
