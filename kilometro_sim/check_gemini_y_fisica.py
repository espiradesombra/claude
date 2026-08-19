"""
Chequeo del script Gemini del Kilometro + modelo fisico honesto.

1) Gemini original (tal cual del chat)
2) Gemini parcheado (completa fases; sigue sin cerrar PE)
3) Modelo cerrado 1a ley (N vueltas enteras, PE vuelve al inicio)
4) Sweeps

Veredicto esperado:
  Gemini "gana" por contabilidad incompleta o por no rearmar la PE.
  En ciclo cerrado real: W_gen < W_motor (perdidas friccion + eficiencias).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
G = 9.81


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    return obj


# ---------------------------------------------------------------------------
# 1) Gemini original
# ---------------------------------------------------------------------------
def run_gemini_original(
    m=50.0,
    r=1.0,
    eta_gen=0.92,
    eta_mot=0.90,
    dt=0.001,
    t_max=5.0,
    omega_target=0.5,
    F_patada=1200.0,
    t_patada_dur=0.2,
    F_avance_gen=300.0,
    dist_avance=2.0,
    omega0=0.1,
):
    theta = 0.0
    omega = omega0
    x_avance = 0.0
    v_avance = 0.0
    E_gen_bajada = 0.0
    E_cons_patada = 0.0
    E_gen_avance = 0.0
    time = 0.0
    fase = "A_BAJADA"
    t_patada_start = 0.0
    steps = 0
    bugs = []

    while time < t_max:
        if fase == "A_BAJADA":
            tau_grav = m * G * r * np.sin(theta)
            # BUG control: si omega < target, tau_gen crece y frena mas
            tau_gen = max(0.0, tau_grav - 0.1 * (omega - omega_target))
            P_gen = tau_gen * omega * eta_gen
            E_gen_bajada += P_gen * dt
            alpha = (tau_grav - tau_gen) / (m * r**2)
            omega += alpha * dt
            theta += omega * dt
            if theta >= np.pi:
                fase = "B_PATADA"
                t_patada_start = time
        elif fase == "B_PATADA":
            if (time - t_patada_start) <= t_patada_dur:
                if v_avance > 0:
                    P_mot = (F_patada * v_avance) / eta_mot
                else:
                    P_mot = (F_patada * 1.0) / eta_mot  # BUG: v ficticia 1 m/s
                    if "fake_v_patada" not in bugs:
                        bugs.append("fake_v_patada")
                E_cons_patada += P_mot * dt
                a_patada = (F_patada - m * G * 0.5) / m
                v_avance += a_patada * dt
                x_avance += v_avance * dt
            else:
                fase = "C_AVANCE_REGENERATIVO"
        elif fase == "C_AVANCE_REGENERATIVO":
            F_grav_fav = m * G * 0.8  # BUG: fuerza gravitatoria constante inventada
            if "const_F_grav_fav" not in bugs:
                bugs.append("const_F_grav_fav")
            a_avance = (F_grav_fav - F_avance_gen) / m if F_grav_fav > F_avance_gen else 0.0
            v_avance += a_avance * dt
            x_avance += v_avance * dt
            E_gen_avance += F_avance_gen * v_avance * eta_gen * dt
            if x_avance >= dist_avance:
                break
        time += dt
        steps += 1

    E_total = E_gen_bajada + E_gen_avance
    completed = fase == "C_AVANCE_REGENERATIVO" and x_avance >= dist_avance
    if not completed and fase == "A_BAJADA":
        bugs.append("stuck_in_bajada")
    return {
        "model": "gemini_original",
        "completed_cycle": completed,
        "fase_final": fase,
        "t": time,
        "steps": steps,
        "theta": theta,
        "omega": omega,
        "x_avance": x_avance,
        "v_avance": v_avance,
        "E_gen_bajada": E_gen_bajada,
        "E_cons_patada": E_cons_patada,
        "E_gen_avance": E_gen_avance,
        "E_total_generada": E_total,
        "E_neto": E_total - E_cons_patada,
        "bugs": bugs,
    }


# ---------------------------------------------------------------------------
# 2) Gemini parcheado (dinamicamente funciona; fisicamente incompleto en PE)
# ---------------------------------------------------------------------------
def run_gemini_patched(
    m=50.0,
    r=1.0,
    eta_gen=0.92,
    eta_mot=0.90,
    dt=0.001,
    t_max=30.0,
    omega_target=0.5,
    F_patada=1200.0,
    t_patada_dur=0.2,
    F_avance_gen=300.0,
    dist_avance=2.0,
    omega0=0.5,
):
    theta = 0.0
    omega = omega0
    x_avance = 0.0
    v_avance = 0.0
    E_gen_bajada = 0.0
    E_cons_patada = 0.0
    E_gen_avance = 0.0
    W_grav_claimed = 0.0
    time = 0.0
    fase = "A_BAJADA"
    t_patada_start = 0.0
    I = m * r**2

    while time < t_max:
        if fase == "A_BAJADA":
            tau_grav = m * G * r * np.sin(theta)
            # regeneracion solo si hay exceso de velocidad
            tau_gen = max(0.0, 8.0 * (omega - omega_target))
            # no superar el torque gravitatorio disponible en exceso
            tau_gen = min(tau_gen, max(tau_grav, 0.0) + 5.0)
            E_gen_bajada += tau_gen * max(omega, 0.0) * eta_gen * dt
            W_grav_claimed += tau_grav * omega * dt
            alpha = (tau_grav - tau_gen) / I
            omega = max(0.05, omega + alpha * dt)  # suelo minimo para avanzar
            theta += omega * dt
            if theta >= np.pi:
                fase = "B_PATADA"
                t_patada_start = time
                v_avance = max(v_avance, omega * r * 0.25)
        elif fase == "B_PATADA":
            if (time - t_patada_start) <= t_patada_dur:
                P_mot = (F_patada * max(v_avance, 0.0)) / eta_mot
                E_cons_patada += P_mot * dt
                a_patada = (F_patada - m * G * 0.5) / m
                v_avance += a_patada * dt
                x_avance += max(v_avance, 0.0) * dt
            else:
                fase = "C_AVANCE_REGENERATIVO"
        elif fase == "C_AVANCE_REGENERATIVO":
            F_grav_fav = m * G * 0.8  # inventada: no cierra PE
            a_avance = (F_grav_fav - F_avance_gen) / m
            v_avance = max(0.0, v_avance + a_avance * dt)
            x_avance += v_avance * dt
            E_gen_avance += F_avance_gen * v_avance * eta_gen * dt
            W_grav_claimed += F_grav_fav * v_avance * dt
            if x_avance >= dist_avance:
                break
        time += dt

    E_total = E_gen_bajada + E_gen_avance
    return {
        "model": "gemini_patched",
        "completed_cycle": x_avance >= dist_avance,
        "fase_final": fase,
        "t": time,
        "theta": theta,
        "omega": omega,
        "x_avance": x_avance,
        "v_avance": v_avance,
        "E_gen_bajada": E_gen_bajada,
        "E_cons_patada": E_cons_patada,
        "E_gen_avance": E_gen_avance,
        "E_total_generada": E_total,
        "E_neto": E_total - E_cons_patada,
        "W_grav_claimed": W_grav_claimed,
        "delta_PE_drop_ref_J": m * G * (2 * r),
        "note": "E_neto>0 NO es sobreunidad: se cobra PE sin rearmar.",
    }


# ---------------------------------------------------------------------------
# 3) Modelo cerrado: masa en circulo, N vueltas enteras
#    PE = m g R sin(phi); tau_g = -m g R cos(phi)
#    Tras N*2pi: Delta PE = 0 => W_grav ~ 0
# ---------------------------------------------------------------------------
@dataclass
class ClosedParams:
    m: float = 50.0
    R: float = 1.0
    I_rotor: float = 20.0
    beta: float = 1.5
    omega_ref: float = 2.0
    # control por sectores
    tau_motor_max: float = 600.0  # "patada" / motor zona desfavorable
    c_gen: float = 80.0  # N·m·s/rad  freno regenerativo en bajada
    eta_gen: float = 0.92
    eta_mot: float = 0.90
    dt: float = 0.00025
    n_turns: float = 5.0


def run_closed_cycle(p: ClosedParams):
    I = p.I_rotor + p.m * p.R**2
    phi = 0.0  # phi=0: altura 0; phi=pi/2: max altura R
    omega = p.omega_ref
    PE0 = p.m * G * p.R * np.sin(phi)
    KE0 = 0.5 * I * omega**2

    W_motor = 0.0
    W_gen = 0.0
    W_fric = 0.0
    W_grav = 0.0
    t = 0.0
    phi_target = p.n_turns * 2.0 * np.pi
    max_steps = int(200.0 / p.dt)  # tope seguridad

    stalled = False
    for _ in range(max_steps):
        # torque gravitatorio (conservativo)
        tau_g = -p.m * G * p.R * np.cos(phi)

        # sectores en [0, 2pi)
        s = phi % (2.0 * np.pi)
        # Subida desfavorable: s in (0, pi) si sin crece luego baja...
        # Con PE = m g R sin(phi):
        #   0->pi/2: sube (desfavorable)  pi/2->3pi/2: baja (favorable? sin baja de 1 a -1)
        #   3pi/2->2pi: sube de -1 a 0
        # Simplificamos:
        #   Favorable (gravedad acelera sentido +omega): dPE/dt < 0 con omega>0
        #   dPE/dt = m g R cos(phi) * omega  => favorable cuando cos(phi)<0 => phi in (pi/2, 3pi/2)
        favorable = np.cos(phi) < 0.0

        if favorable:
            # regeneracion proporcional a omega (cosecha)
            tau_c = -p.c_gen * max(omega, 0.0)
            mode = "gen"
        else:
            # motor mantiene omega_ref en zona dura ("patada" distribuida)
            err = p.omega_ref - omega
            tau_c = min(p.tau_motor_max, max(0.0, 120.0 * err + 0.5 * p.tau_motor_max))
            mode = "motor"

        tau_fric = -p.beta * omega
        alpha = (tau_g + tau_c + tau_fric) / I

        # trabajos con omega actual (antes de actualizar)
        W_grav += tau_g * omega * p.dt
        W_fric += abs(tau_fric * omega * p.dt)
        if mode == "motor":
            W_motor += (max(tau_c * omega, 0.0) / p.eta_mot) * p.dt
        else:
            W_gen += (max(-tau_c * omega, 0.0) * p.eta_gen) * p.dt

        omega_new = omega + alpha * p.dt
        # Si se para o invierte: fin del ensayo (sin suelo magico de energia)
        if omega_new <= 0.0:
            stalled = True
            omega = 0.0
            t += p.dt
            break

        phi = phi + omega * p.dt
        omega = omega_new
        t += p.dt
        if phi >= phi_target:
            break

    PE1 = p.m * G * p.R * np.sin(phi)
    KE1 = 0.5 * I * omega**2
    dPE = PE1 - PE0
    dKE = KE1 - KE0
    turns = phi / (2.0 * np.pi)

    # 1a ley mecanica+electrica:
    # W_motor_elec * eta paths ya estan en W_motor (electrica in) y W_gen (electrica out)
    # Balance mecanico: dKE + dPE = W_tau_c_mec + W_fric_signed
    # Contabilidad electrica reportada:
    #   entradas: W_motor (+ caida PE si dPE<0 cuenta como "source")
    # Verificamos residual energetico total:
    # W_motor + (-dPE if extracted) ... mejor:
    # Sistema: d(KE+PE) = P_motor_mec + P_gen_mec + P_fric
    # P_motor_mec = tau_c*omega >0, P_gen_mec = tau_c*omega <0, P_fric = tau_fric*omega <0
    # W_motor_elec = integral P_motor_mec/eta
    # W_gen_elec = integral |P_gen_mec|*eta
    # Residual check (mec): dKE + dPE - W_c_mec - W_fric_signed ~= 0
    # Usamos residual 1a ley en forma:
    # W_motor_elec * (approx) no cierra exacto por eta; residual mecanico puro:

    # Residual que usamos en el repo: W_motor + W_grav ~= W_gen + W_fric + dKE
    # (W_grav = integral tau_g omega = -dPE en limite exacto)
    lhs = W_motor + W_grav
    rhs = W_gen + W_fric + dKE
    # Nota: W_motor es electrico (mayor que mecanico), W_gen electrico (menor que mecanico)
    # por eso residual no es cero si se mezclan dominios. Mejor residual mecanico:

    # Recomputamos residual mecanico aproximado:
    # W_motor_mec ~ W_motor * eta_mot, W_gen_mec ~ W_gen / eta_gen
    W_motor_mec = W_motor * p.eta_mot
    W_gen_mec = W_gen / p.eta_gen if p.eta_gen > 0 else 0.0
    # dKE + dPE = W_motor_mec - W_gen_mec + W_fric_signed
    # W_fric_signed = -W_fric
    resid_mec = abs((dKE + dPE) - (W_motor_mec - W_gen_mec - W_fric))
    scale = max(abs(dKE) + abs(dPE) + W_motor_mec + W_gen_mec + W_fric, 1e-9)

    W_net = W_gen - W_motor  # electrico neto (debe ser < 0)

    return {
        "model": "closed_cycle_1st_law",
        "params": asdict(p),
        "turns": turns,
        "t": t,
        "omega_end": omega,
        "stalled_floor_hit": stalled,
        "W_motor_elec": W_motor,
        "W_gen_elec": W_gen,
        "W_fric": W_fric,
        "W_grav": W_grav,
        "dPE": dPE,
        "dKE": dKE,
        "W_net_elec": W_net,
        "W_grav_plus_dPE": W_grav + dPE,  # ~0 si integracion buena
        "mech_residual_J": resid_mec,
        "mech_residual_pct": 100.0 * resid_mec / scale,
        "overunity_elec": W_net > 0.0,
        "eta_paid": (W_gen / W_motor) if W_motor > 1e-9 else float("nan"),
    }


def sweep_closed():
    rows = []
    for omega_ref in [1.0, 2.0, 3.0]:
        for tau_m in [200.0, 600.0, 1200.0]:
            for c_gen in [20.0, 80.0, 200.0]:
                for beta in [0.5, 1.5, 4.0]:
                    p = ClosedParams(
                        omega_ref=omega_ref,
                        tau_motor_max=tau_m,
                        c_gen=c_gen,
                        beta=beta,
                        n_turns=5.0,
                    )
                    r = run_closed_cycle(p)
                    # Correccion por drenaje de bateria cinetica inicial
                    ke_drain = max(0.0, -r["dKE"])
                    surplus_after_ke = r["W_net_elec"] - p.eta_gen * ke_drain
                    valid = (
                        (not r["stalled_floor_hit"])
                        and r["turns"] >= p.n_turns * 0.99
                        and r["mech_residual_pct"] < 2.0
                        and abs(r["W_grav"] + r["dPE"]) < 50.0
                    )
                    rows.append(
                        {
                            "omega_ref": omega_ref,
                            "tau_motor_max": tau_m,
                            "c_gen": c_gen,
                            "beta": beta,
                            "W_net_elec": r["W_net_elec"],
                            "surplus_after_ke": surplus_after_ke,
                            "W_gen": r["W_gen_elec"],
                            "W_motor": r["W_motor_elec"],
                            "W_grav": r["W_grav"],
                            "dPE": r["dPE"],
                            "dKE": r["dKE"],
                            "turns": r["turns"],
                            "eta_paid": r["eta_paid"],
                            "resid_pct": r["mech_residual_pct"],
                            "overunity_raw": r["overunity_elec"],
                            "overunity_after_ke": surplus_after_ke > 1.0,
                            "valid_closed": valid,
                            "stalled": r["stalled_floor_hit"],
                        }
                    )
    valid_rows = [x for x in rows if x["valid_closed"]]
    rows_sorted = sorted(valid_rows or rows, key=lambda x: x["surplus_after_ke"], reverse=True)
    return {
        "n_combos": len(rows),
        "n_valid": len(valid_rows),
        "n_positive_raw": sum(1 for x in rows if x["overunity_raw"]),
        "n_positive_after_ke_valid": sum(
            1 for x in valid_rows if x["overunity_after_ke"]
        ),
        "best5": rows_sorted[:5],
        "worst5": rows_sorted[-5:],
        "max_surplus_after_ke": rows_sorted[0]["surplus_after_ke"] if rows_sorted else None,
        "min_surplus_after_ke": rows_sorted[-1]["surplus_after_ke"] if rows_sorted else None,
        "best_eta_paid_valid": max(
            (x["eta_paid"] for x in valid_rows if x["eta_paid"] == x["eta_paid"]),
            default=float("nan"),
        ),
    }


def sweep_gemini_patched():
    rows = []
    for F_patada in [200.0, 600.0, 1200.0, 2500.0]:
        for omega_target in [0.3, 0.5, 1.0]:
            for F_avance_gen in [100.0, 300.0, 500.0]:
                r = run_gemini_patched(
                    F_patada=F_patada,
                    omega_target=omega_target,
                    F_avance_gen=F_avance_gen,
                    t_max=40.0,
                )
                rows.append(
                    {
                        "F_patada": F_patada,
                        "omega_target": omega_target,
                        "F_avance_gen": F_avance_gen,
                        "completed": r["completed_cycle"],
                        "E_neto": r["E_neto"],
                        "E_gen": r["E_total_generada"],
                        "E_patada": r["E_cons_patada"],
                        "t": r["t"],
                    }
                )
    rows_sorted = sorted(rows, key=lambda x: x["E_neto"], reverse=True)
    return {
        "n_combos": len(rows),
        "n_positive": sum(1 for x in rows if x["E_neto"] > 0),
        "n_completed": sum(1 for x in rows if x["completed"]),
        "best5": rows_sorted[:5],
        "note": "Positivos Gemini = PE no rearmada, no motor perpetuo.",
    }


def main():
    print("=" * 64)
    print("1) GEMINI ORIGINAL (parametros del chat)")
    print("=" * 64)
    r0 = run_gemini_original()
    for k, v in r0.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")

    print("\n" + "=" * 64)
    print("2) GEMINI PARCHADO (completa fases; PE incompleta)")
    print("=" * 64)
    r1 = run_gemini_patched()
    for k, v in r1.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")

    print("\n" + "=" * 64)
    print("3) MODELO CERRADO 1a LEY (5 vueltas enteras)")
    print("=" * 64)
    r2 = run_closed_cycle(ClosedParams())
    for k, v in r2.items():
        if k == "params":
            continue
        print(f"  {k}: {v:.6f}" if isinstance(v, float) else f"  {k}: {v}")

    print("\n" + "=" * 64)
    print("4) SWEEP modelo cerrado")
    print("=" * 64)
    sw = sweep_closed()
    print(f"  combos: {sw['n_combos']}  validos: {sw['n_valid']}")
    print(f"  W_net raw > 0: {sw['n_positive_raw']}")
    print(f"  surplus_after_ke > 0 (validos): {sw['n_positive_after_ke_valid']}")
    print(f"  max surplus_after_ke: {sw['max_surplus_after_ke']:.4f} J")
    print(f"  min surplus_after_ke: {sw['min_surplus_after_ke']:.4f} J")
    print(f"  best eta_paid validos: {sw['best_eta_paid_valid']:.4f}")
    print("  best5 validos (surplus tras descontar drenaje KE):")
    for row in sw["best5"]:
        print(
            f"    w={row['omega_ref']}, tau_m={row['tau_motor_max']}, "
            f"c_gen={row['c_gen']}, beta={row['beta']} -> "
            f"W_net={row['W_net_elec']:.2f} J, after_KE={row['surplus_after_ke']:.2f}, "
            f"eta_paid={row['eta_paid']:.3f}, turns={row['turns']:.2f}, "
            f"valid={row['valid_closed']}"
        )

    print("\n" + "=" * 64)
    print("5) SWEEP Gemini parcheado")
    print("=" * 64)
    sg = sweep_gemini_patched()
    print(f"  combos: {sg['n_combos']}")
    print(f"  E_neto>0: {sg['n_positive']}  completed: {sg['n_completed']}")
    print(f"  note: {sg['note']}")
    for row in sg["best5"]:
        print(
            f"    Fpat={row['F_patada']}, w={row['omega_target']}, "
            f"Fgen={row['F_avance_gen']} -> E_neto={row['E_neto']:.1f} J "
            f"done={row['completed']}"
        )

    m, R = 50.0, 1.0
    dPE = m * G * (2 * R)
    print("\n" + "=" * 64)
    print("REFERENCIA + VEREDICTO")
    print("=" * 64)
    print(f"  Delta PE caida diametro 2R (m=50,R=1): {dPE:.1f} J  (una sola vez)")
    print(
        "  Gemini original: NO completa el ciclo; 'exito' falso por bajada parcial."
    )
    print(
        "  Gemini parcheado: completa fases pero inventa F_grav_fav y no rearma PE."
    )
    print(
        "  Cerrado 1a ley: W_net electrico debe ser <= 0; eta_paid < 1."
    )
    print(
        "  Doble-km por perneo: acoplamiento inercial / supercondensador mecanico;"
    )
    print(
        "  reparte picos al despernear, no crea energia de la gravedad."
    )

    report = {
        "gemini_original": r0,
        "gemini_patched": r1,
        "closed_cycle": {k: v for k, v in r2.items() if k != "params"},
        "closed_params": r2["params"],
        "sweep_closed": sw,
        "sweep_gemini_patched": sg,
        "delta_PE_drop_J": dPE,
        "verdict_es": (
            "El script de Gemini no valida un excedente fisico real. "
            "Con parametros del chat ni siquiera sale de la fase A. "
            "Parcheado da E_neto>0 porque cobra energia potencial sin devolverla. "
            "En modelo de ciclo cerrado (N vueltas), el balance electrico neto "
            "es negativo: Kilometro = bateria/convertidor, no motor perpetuo. "
            "Doble kilometro por perneo = descarga inercial acoplada (util como "
            "buffer de potencia), no generacion neta de julios."
        ),
    }
    out_path = OUT / "resultados_chequeo.json"
    out_path.write_text(
        json.dumps(_jsonable(report), indent=2), encoding="utf-8"
    )
    print(f"\nJSON: {out_path}")


if __name__ == "__main__":
    main()
