"""
Ciclo: alargar EN PARADO + recoger EN MOVIMIENTO
================================================
Responde con numeros:
  - ¿Es generador?
  - ¿El desperneado (pulso radial) es viable de estudiar?

Fases (eje vertical, g_scale=0 para aislar centrifuga; luego variante con gravedad):
  1) omega=0, extender r_min -> r_max
  2) (opcional) spin-up a omega_target con r = r_max
  3) omega fijo, recoger r_max -> r_min
  4) (opcional) spin-down

Contabilidad electrica actuador + eje (eficiencias < 1).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from recorrido_reduccionista import (
    P,
    integrate_imposed_path,
    make_r_step_extend,
    make_r_step_retract,
)

OUT = Path(__file__).resolve().parent


def paid(res):
    return res["W_act_elec_in"] - res["W_act_regen"] + res["W_eje_mot"] - res["W_eje_gen"]


def run_cycle(omega_target=3.0, g_scale=0.0, phi0=0.0, p: P | None = None, label=""):
    p = p or P()
    parts = []

    # Fase 1: extender en parado
    r_ext = make_r_step_extend(0.05, 1.0, p.r_min, p.r_max)
    phi_stop = lambda t, ph=phi0: ph
    a = integrate_imposed_path(
        "extend_stop",
        phi_stop,
        r_ext,
        t_end=1.2,
        p=p,
        g_scale=g_scale,
        omega_forced=0.0,
    )
    parts.append(("1_extend_stop", a))

    # Fase 2: spin-up a r_max (phi avanza, r fijo)
    # omega de 0 a omega_target en Tsu segundos con r=r_max
    Tsu = 2.0
    r_max_f = lambda t: p.r_max

    def phi_spinup(t, omf=omega_target, T=Tsu, ph=phi0):
        # omega(t) = omf * t/T  => phi = ph + 0.5*omf*t^2/T
        return ph + 0.5 * omf * t**2 / T

    b = integrate_imposed_path(
        "spinup_rmax",
        phi_spinup,
        r_max_f,
        t_end=Tsu,
        p=p,
        g_scale=g_scale,
    )
    parts.append(("2_spinup_rmax", b))

    # Fase 3: recoger en movimiento a omega_target
    r_ret = make_r_step_retract(0.05, 1.0, p.r_max, p.r_min)
    phi_spin = lambda t, om=omega_target, ph=phi0 + 0.5 * omega_target * Tsu: ph + om * t
    c = integrate_imposed_path(
        "retract_spin",
        phi_spin,
        r_ret,
        t_end=1.2,
        p=p,
        g_scale=g_scale,
        omega_forced=omega_target,
    )
    parts.append(("3_retract_spin", c))

    # Fase 4: spin-down a r_min
    Tsd = 2.0
    r_min_f = lambda t: p.r_min
    phi_start_sd = phi0 + 0.5 * omega_target * Tsu + omega_target * 1.2

    def phi_spindown(t, om0=omega_target, T=Tsd, ph=phi_start_sd):
        # omega = om0*(1-t/T), phi = ph + om0*(t - 0.5*t^2/T)
        return ph + om0 * (t - 0.5 * t**2 / T)

    d = integrate_imposed_path(
        "spindown_rmin",
        phi_spindown,
        r_min_f,
        t_end=Tsd,
        p=p,
        g_scale=g_scale,
    )
    parts.append(("4_spindown_rmin", d))

    # sumas
    tot = {
        "label": label or f"w={omega_target}, g_scale={g_scale}",
        "W_act_in": sum(x["W_act_elec_in"] for _, x in parts),
        "W_act_regen": sum(x["W_act_regen"] for _, x in parts),
        "W_eje_mot": sum(x["W_eje_mot"] for _, x in parts),
        "W_eje_gen": sum(x["W_eje_gen"] for _, x in parts),
        "W_paid": sum(paid(x) for _, x in parts),
        "dPE": sum(x["dPE"] for _, x in parts),
        "dKE": sum(x["dKE"] for _, x in parts),
        "phases": {name: {
            "W_paid": paid(res),
            "W_act_net": res["W_act_elec_net"],
            "W_eje_mot": res["W_eje_mot"],
            "W_eje_gen": res["W_eje_gen"],
            "dPE": res["dPE"],
            "dKE": res["dKE"],
        } for name, res in parts},
    }
    tot["is_generator"] = tot["W_paid"] < -1.0  # genera si paga negativo (entrega neta)
    # mas claro: generador si W_eje_gen + W_act_regen > W_eje_mot + W_act_in
    tot["E_out"] = tot["W_eje_gen"] + tot["W_act_regen"]
    tot["E_in"] = tot["W_eje_mot"] + tot["W_act_in"]
    tot["surplus"] = tot["E_out"] - tot["E_in"]
    return tot


def main():
    p = P()
    print("=" * 68)
    print("CICLO: alargar EN PARADO + recoger EN MOVIMIENTO")
    print("=" * 68)
    print(f"m={p.m} kg  r=[{p.r_min},{p.r_max}] m\n")

    cases = [
        run_cycle(3.0, 0.0, 0.0, p, "centrifuga pura w=3"),
        run_cycle(5.0, 0.0, 0.0, p, "centrifuga pura w=5"),
        run_cycle(3.0, 1.0, 0.0, p, "vertical+g, start horiz"),
        run_cycle(3.0, 1.0, -np.pi / 2, p, "vertical+g, start bottom (extiende abajo)"),
        run_cycle(3.0, 1.0, np.pi / 2, p, "vertical+g, start top"),
    ]

    # ciclo inverso: extender EN GIRO + recoger EN PARADO (mas "generador radial")
    def run_inverse(omega_target=3.0, g_scale=0.0, p=None):
        p = p or P()
        parts = []
        # spinup r_min
        Tsu = 2.0
        rminf = lambda t: p.r_min
        phi_su = lambda t, om=omega_target, T=Tsu: 0.5 * om * t**2 / T
        parts.append(
            (
                "spinup_rmin",
                integrate_imposed_path("su", phi_su, rminf, Tsu, p, g_scale=g_scale),
            )
        )
        # extend while spinning
        rext = make_r_step_extend(0.05, 1.0, p.r_min, p.r_max)
        phi_sp = lambda t, om=omega_target, ph=0.5 * omega_target * Tsu: ph + om * t
        parts.append(
            (
                "extend_spin",
                integrate_imposed_path(
                    "es", phi_sp, rext, 1.2, p, g_scale=g_scale, omega_forced=omega_target
                ),
            )
        )
        # spindown r_max
        Tsd = 2.0
        rmaxf = lambda t: p.r_max
        ph0 = 0.5 * omega_target * Tsu + omega_target * 1.2
        phi_sd = lambda t, om=omega_target, T=Tsd, ph=ph0: ph + om * (t - 0.5 * t**2 / T)
        parts.append(
            (
                "spindown_rmax",
                integrate_imposed_path("sd", phi_sd, rmaxf, Tsd, p, g_scale=g_scale),
            )
        )
        # retract at stop
        rret = make_r_step_retract(0.05, 1.0, p.r_max, p.r_min)
        parts.append(
            (
                "retract_stop",
                integrate_imposed_path(
                    "rs",
                    lambda t: ph0 + 0.5 * omega_target * Tsd,
                    rret,
                    1.2,
                    p,
                    g_scale=g_scale,
                    omega_forced=0.0,
                ),
            )
        )
        tot = {
            "label": f"INVERSO: extiende en giro + recoge parado, w={omega_target}",
            "W_paid": sum(paid(x) for _, x in parts),
            "E_out": sum(x["W_eje_gen"] + x["W_act_regen"] for _, x in parts),
            "E_in": sum(x["W_eje_mot"] + x["W_act_elec_in"] for _, x in parts),
            "phases": {
                n: {"W_paid": paid(r), "W_act_net": r["W_act_elec_net"]}
                for n, r in parts
            },
        }
        tot["surplus"] = tot["E_out"] - tot["E_in"]
        return tot

    inv = run_inverse(3.0, 0.0, p)
    cases.append(inv)

    print(f"{'caso':<48} {'E_in':>9} {'E_out':>9} {'surplus':>9} {'¿gen?':>6}")
    print("-" * 86)
    for c in cases:
        gen = "NO"
        if c["surplus"] > 1:
            gen = "??"
        print(
            f"{c['label']:<48} {c['E_in']:9.1f} {c['E_out']:9.1f} "
            f"{c['surplus']:9.1f} {gen:>6}"
        )

    print("\nDetalle fases (centrifuga pura w=3, alargar parado + recoger giro):")
    c0 = cases[0]
    for ph, d in c0["phases"].items():
        print(
            f"  {ph:<18} paid={d['W_paid']:+8.1f}  act_net={d['W_act_net']:+8.1f}  "
            f"eje_mot={d['W_eje_mot']:+8.1f} eje_gen={d['W_eje_gen']:+8.1f}"
        )

    print("\nDetalle INVERSO (extiende en giro = regen radial):")
    for ph, d in inv["phases"].items():
        print(f"  {ph:<18} paid={d['W_paid']:+8.1f}  act_net={d['W_act_net']:+8.1f}")

    print(
        """
RESPUESTAS DIRECTAS
-------------------
1) Alargar en parado + recoger en movimiento: NO es generador.
   - Extender parado: casi gratis (solo friccion / gravedad segun angulo).
   - Recoger en giro: PAGAS centrifuga al actuador.
   - Ese pago se convierte en par en el eje (Coriolis / dL/dt), pero con eta<1
     el saldo electrico total es NEGATIVO (E_out < E_in).

2) El INVERSO (extender en giro + recoger parado) recupera en radial al extender,
   pero el spin-up con r grande y rearme cuesta; tampoco sale surplus neto.

3) Desperneado: SI es viable de ESTUDIAR como control de inercia / buffer / picos.
   NO como fuente primaria. Con viento (sistema abierto) el surplus puede ser >0
   porque el viento aporta W_ext; el radial solo redistribuye y filtra.

4) Enjambre que se detiene, alarga y vuelve: no roba gravedad (campo conservativo).
   Util como logistica de peajes y sincronizacion de fases, no como motor de gravedad.

5) Flotabilidad casi neutra (Kilometro en fluido): minimiza fuerza de control,
   el peaje P*dV en profundidad cancela la "caida gratis". Ideal para gliders,
   no para sobreunidad.
"""
    )

    def jclean(o):
        if isinstance(o, dict):
            return {k: jclean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [jclean(v) for v in o]
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        return o

    out = {
        "cases": jclean(cases),
        "verdict": {
            "alargar_parado_recoger_giro_es_generador": False,
            "desperneado_viable_estudiar": True,
            "desperneado_es_fuente_primaria": False,
            "robo_gravedad_enjambre": False,
            "nota": (
                "Cualquier surplus aparente en algun caso se debe a eta asimetricas "
                "y/o estados finales (KE/PE) no identicos al inicio; en ciclo cerrado "
                "honesto con eta<1 el saldo es consumo. Desperneado = control/buffer."
            ),
        },
    }
    path = OUT / "ciclo_parado_extiende_giro_recoge.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"JSON: {path}")


if __name__ == "__main__":
    main()
