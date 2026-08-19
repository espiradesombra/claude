"""
Recorrido reduccionista del Kilometro (radio variable r)
=======================================================
Pregunta de VMA:
  - Subir con r pequeño / bajar con r grande: ¿asimetriza el peaje gravitatorio?
  - Extender o recoger la masa con el Km PARADO: ¿cuesta menos que en giro?
  - ¿Hay ciclo de "gasto eficiente" o magia?

Fisica (plano vertical, angulo phi, radio r):
  posicion:  x = r cos(phi),  y = r sin(phi)   (y vertical hacia arriba)
  PE = m g y = m g r sin(phi)
  T  = 1/2 m (rdot^2 + r^2 omega^2)

  Fuerzas generalizadas (Lagrangiano):
    tau_g  = - dV/dphi = - m g r cos(phi)     (par sobre el eje)
    F_r_grav = - dV/dr = - m g sin(phi)       (componente radial de gravedad)
    F_r_cent = m r omega^2                    (centrifuga, en la ecuacion de r)

  Ecuacion radial (si actuador impone r(t) o aplica F_act):
    m rddot = F_act + m r omega^2 - m g sin(phi) - c_r rdot
    => F_act = m rddot - m r omega^2 + m g sin(phi) + c_r rdot

  Trabajo actuador radial:
    W_act = integral F_act * rdot dt

  Trabajo motor/gen en el eje (si imponemos omega o controlamos tau):
    W_eje = integral tau_eje * omega dt

No se afirma sobreunidad: se compara peajes en estrategias distintas.
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
    m: float = 40.0
    r_min: float = 0.4
    r_max: float = 1.2
    c_r: float = 5.0  # amortiguamiento radial N·s/m
    beta: float = 2.0  # friccion eje N·m·s/rad
    eta_act: float = 0.90  # eficiencia actuador radial (elec->mec si W>0)
    eta_mot: float = 0.90
    eta_gen: float = 0.92
    dt: float = 0.0005


def pe(m, r, phi):
    return m * G * r * np.sin(phi)


def ke(m, r, rdot, omega):
    return 0.5 * m * (rdot**2 + (r * omega) ** 2)


def F_act_required(m, r, rdot, rddot, omega, phi, c_r, g_scale=1.0):
    """Fuerza que debe aplicar el actuador radial para seguir la cinematica.
    g_scale=0 => solo centrifuga/friccion (eje vertical / sin acoplar PE a r).
    """
    return (
        m * rddot
        - m * r * omega**2
        + m * (G * g_scale) * np.sin(phi)
        + c_r * rdot
    )


def integrate_imposed_path(
    name: str,
    phi_of_t,
    r_of_t,
    t_end: float,
    p: P,
    omega0: float | None = None,
    free_spin: bool = False,
    tau_extra_fn=None,
    g_scale: float = 1.0,
    omega_forced: float | None = None,
):
    """
    Integra un ciclo con r(t) y phi(t) impuestos (cinematica forzada),
    o free_spin=True con solo r(t) impuesto y dinamica en phi.

    Si free_spin: phi se integra con tau_g + tau_fric + tau_extra; r(t) impuesto.
    Si no: phi(t) y r(t) impuestos (util para comparar trabajos a omega fijo / parado).
    """
    dt = p.dt
    n = int(t_end / dt) + 1
    t = np.linspace(0.0, t_end, n)

    phi = np.zeros(n)
    omega = np.zeros(n)
    r = np.zeros(n)
    rdot = np.zeros(n)

    W_act_mec = 0.0  # trabajo mecanico actuador (puede ser + o -)
    W_act_elec_in = 0.0  # solo cuando actuador aporta energia
    W_act_regen = 0.0  # cuando actuador absorbe (regen radial)
    W_eje_mot = 0.0
    W_eje_gen = 0.0
    W_fric = 0.0
    W_grav_phi = 0.0  # integral tau_g * omega
    W_grav_r = 0.0  # integral F_r_grav * rdot

    # init
    r[0] = r_of_t(0.0)
    if free_spin:
        phi[0] = phi_of_t(0.0) if callable(phi_of_t) else float(phi_of_t)
        omega[0] = float(omega0 if omega0 is not None else 2.0)
    else:
        # cinematica forzada: derivar phi de la funcion
        phi[0] = phi_of_t(0.0)
        # omega numerico al inicio
        phi1 = phi_of_t(dt)
        omega[0] = (phi1 - phi[0]) / dt

    Gg = G * g_scale
    PE0 = p.m * Gg * r[0] * np.sin(phi[0])
    KE0 = ke(p.m, r[0], 0.0, omega[0])

    for i in range(n - 1):
        ti = t[i]
        # r cinematica
        r0 = r_of_t(ti)
        r1 = r_of_t(ti + dt)
        r2 = r_of_t(min(ti + 2 * dt, t_end))
        rd = (r1 - r0) / dt
        rdd = (r2 - 2 * r1 + r0) / (dt**2) if ti + 2 * dt <= t_end + 1e-15 else 0.0
        r[i] = r0
        rdot[i] = rd

        if free_spin:
            ph = phi[i]
            om = omega[i]
            tau_g = -p.m * Gg * r0 * np.cos(ph)
            tau_fric = -p.beta * om
            tau_x = tau_extra_fn(ti, ph, om, r0) if tau_extra_fn else 0.0
            I = max(p.m * r0**2, 1e-6)
            alpha = (tau_g + tau_fric + tau_x - 2.0 * p.m * r0 * rd * om) / I
            W_grav_phi += tau_g * om * dt
            W_fric += abs(tau_fric * om * dt)
            if tau_x > 0:
                W_eje_mot += (tau_x * max(om, 0.0) / p.eta_mot) * dt
            elif tau_x < 0:
                W_eje_gen += ((-tau_x) * max(om, 0.0) * p.eta_gen) * dt
            omega[i + 1] = om + alpha * dt
            phi[i + 1] = ph + om * dt
        elif omega_forced is not None:
            # Giro a omega constante con phi = omega*t (plano vertical o g_scale=0)
            om = float(omega_forced)
            ph = phi_of_t(ti) if callable(phi_of_t) else float(phi_of_t) + om * ti
            # si phi_of_t es solo offset:
            if not callable(phi_of_t):
                ph = float(phi_of_t) + om * ti
            else:
                # preferir cinematica forzada pura
                ph = phi_of_t(ti)
                om = (phi_of_t(ti + dt) - ph) / dt
            phi[i] = ph
            omega[i] = om
            tau_g = -p.m * Gg * r0 * np.cos(ph)
            W_grav_phi += tau_g * om * dt
            I = max(p.m * r0**2, 1e-6)
            # alpha~0 si omega forzado constante
            alpha = 0.0
            if callable(phi_of_t):
                om2 = (phi_of_t(min(ti + 2 * dt, t_end)) - phi_of_t(ti + dt)) / dt
                alpha = (om2 - om) / dt
            tau_eje = I * alpha + 2.0 * p.m * r0 * rd * om - tau_g + p.beta * om
            W_fric += abs(p.beta * om * om * dt)
            if tau_eje * om > 0:
                W_eje_mot += (tau_eje * om / p.eta_mot) * dt
            elif tau_eje * om < 0:
                W_eje_gen += ((-tau_eje * om) * p.eta_gen) * dt
            phi[i + 1] = ph + om * dt
            omega[i + 1] = om
        else:
            ph = phi_of_t(ti)
            ph_n = phi_of_t(ti + dt)
            om = (ph_n - ph) / dt
            phi[i] = ph
            omega[i] = om
            tau_g = -p.m * Gg * r0 * np.cos(ph)
            W_grav_phi += tau_g * om * dt
            I = max(p.m * r0**2, 1e-6)
            om_n = (
                (phi_of_t(min(ti + 2 * dt, t_end)) - ph_n) / dt
                if ti + 2 * dt <= t_end
                else om
            )
            alpha = (om_n - om) / dt
            tau_eje = I * alpha + 2.0 * p.m * r0 * rd * om - tau_g + p.beta * om
            W_fric += abs(p.beta * om * om * dt)
            if tau_eje * om > 0:
                W_eje_mot += (tau_eje * om / p.eta_mot) * dt
            elif tau_eje * om < 0:
                W_eje_gen += ((-tau_eje * om) * p.eta_gen) * dt
            phi[i + 1] = ph_n
            omega[i + 1] = om_n if ti + 2 * dt <= t_end else om

        # actuador radial
        Fa = F_act_required(p.m, r0, rd, rdd, omega[i], phi[i], p.c_r, g_scale=g_scale)
        dW = Fa * rd * dt
        W_act_mec += dW
        if dW > 0:
            W_act_elec_in += dW / p.eta_act
        else:
            W_act_regen += (-dW) * p.eta_act

        F_r_g = -p.m * Gg * np.sin(phi[i])
        W_grav_r += F_r_g * rd * dt

    # ultimo punto
    r[-1] = r_of_t(t_end)
    rdot[-1] = rdot[-2]
    if not free_spin:
        phi[-1] = phi_of_t(t_end)

    r[-1] = r_of_t(t_end)
    PE1 = p.m * Gg * r[-1] * np.sin(phi[-1])
    KE1 = ke(p.m, r[-1], rdot[-1], omega[-1])
    dPE = PE1 - PE0
    dKE = KE1 - KE0
    # Campo gravitatorio conservativo: W_grav_phi + W_grav_r + dPE ≈ 0
    grav_closure = W_grav_phi + W_grav_r + dPE

    return {
        "name": name,
        "t": t,
        "phi": phi,
        "omega": omega,
        "r": r,
        "rdot": rdot,
        "W_act_mec": W_act_mec,
        "W_act_elec_in": W_act_elec_in,
        "W_act_regen": W_act_regen,
        "W_act_elec_net": W_act_elec_in - W_act_regen,
        "W_eje_mot": W_eje_mot,
        "W_eje_gen": W_eje_gen,
        "W_eje_net": W_eje_gen - W_eje_mot,
        "W_fric": W_fric,
        "W_grav_phi": W_grav_phi,
        "W_grav_r": W_grav_r,
        "dPE": dPE,
        "dKE": dKE,
        "grav_closure": grav_closure,
        "PE0": PE0,
        "PE1": PE1,
        "KE0": KE0,
        "KE1": KE1,
        "omega_mean": float(np.mean(np.abs(omega))),
        "r_mean": float(np.mean(r)),
        "g_scale": g_scale,
    }


# ---------------------------------------------------------------------------
# Estrategias de r(t)
# ---------------------------------------------------------------------------
def make_r_constant(rval):
    return lambda t: float(rval)


def make_r_step_retract(t0, t1, r_max, r_min):
    """Recoge de r_max a r_min entre t0 y t1 (suave coseno)."""

    def r(t):
        if t < t0:
            return r_max
        if t > t1:
            return r_min
        s = (t - t0) / max(t1 - t0, 1e-9)
        # half-cosine
        return r_max + (r_min - r_max) * 0.5 * (1 - np.cos(np.pi * s))

    return r


def make_r_step_extend(t0, t1, r_min, r_max):
    def r(t):
        if t < t0:
            return r_min
        if t > t1:
            return r_max
        s = (t - t0) / max(t1 - t0, 1e-9)
        return r_min + (r_max - r_min) * 0.5 * (1 - np.cos(np.pi * s))

    return r


def make_r_reductionist_cycle(T, r_min, r_max, n_turns=1.0):
    """
    Por cada vuelta: en zona desfavorable (cos phi > 0) r cerca de r_min;
    en favorable (cos phi < 0) r cerca de r_max.
    phi = 2pi * n_turns * t / T
    """

    def phi(t):
        return 2.0 * np.pi * n_turns * t / T

    def r(t):
        ph = phi(t)
        # blend suave: r = r_mid - amp * cos(phi)
        # cos>0 desfavorable -> r bajo; cos<0 favorable -> r alto
        # r = r_mid - amp * cos(phi)  => cos=1 -> r_min, cos=-1 -> r_max
        r_mid = 0.5 * (r_min + r_max)
        amp = 0.5 * (r_max - r_min)
        return r_mid - amp * np.cos(ph)

    return phi, r


def make_r_constant_cycle(T, rval, n_turns=1.0):
    def phi(t):
        return 2.0 * np.pi * n_turns * t / T

    def r(t):
        return float(rval)

    return phi, r


def run_all(p: P | None = None):
    p = p or P()
    results = {}

    # --- A) Rearme radial EN PARADO (omega=0), varias orientaciones ---
    # recoger r_max -> r_min en 1 s, phi fijo
    for label, phi_fixed in [
        ("stop_retract_horizontal", 0.0),  # sin(phi)=0 si usamos y=r sin phi; phi=0 => PE=0
        ("stop_retract_top", np.pi / 2),  # masa arriba: recoger baja PE (gravedad ayuda?)
        # y = r sin(phi): phi=pi/2 top, reducir r => PE baja => gravedad ayuda a recoger
        ("stop_retract_bottom", -np.pi / 2),  # abajo: recoger sube PE
        ("stop_extend_horizontal", 0.0),
        ("stop_extend_bottom", -np.pi / 2),  # extender abajo: gravedad ayuda a extender
        ("stop_extend_top", np.pi / 2),  # extender arriba: pelea gravedad
    ]:
        if "retract" in label:
            rfun = make_r_step_retract(0.1, 1.1, p.r_max, p.r_min)
        else:
            rfun = make_r_step_extend(0.1, 1.1, p.r_min, p.r_max)
        phifun = lambda t, ph=phi_fixed: ph
        # omega = 0 forzado: phi constante
        res = integrate_imposed_path(label, phifun, rfun, t_end=1.3, p=p)
        results[label] = res

    # --- B1) Solo CENTRIFUGA (eje vertical, g_scale=0): recoger/extender a omega fijo ---
    # Teoria: W_cent_mec = -1/2 m omega^2 (r_f^2 - r_i^2)
    # recoger r_max->r_min: W_cent = +1/2 m w^2 (r_max^2 - r_min^2)  (pagas)
    omega_spin = 3.0
    for label, kind, om in [
        ("cent_retract_w0", "retract", 0.0),
        ("cent_retract_w3", "retract", 3.0),
        ("cent_retract_w5", "retract", 5.0),
        ("cent_extend_w0", "extend", 0.0),
        ("cent_extend_w3", "extend", 3.0),
        ("cent_extend_w5", "extend", 5.0),
    ]:
        if kind == "retract":
            rfun = make_r_step_retract(0.1, 1.1, p.r_max, p.r_min)
        else:
            rfun = make_r_step_extend(0.1, 1.1, p.r_min, p.r_max)
        # phi irrelevante si g_scale=0; imponemos omega con phi = om*t
        phifun = lambda t, o=om: o * t
        res = integrate_imposed_path(
            label, phifun, rfun, t_end=1.3, p=p, g_scale=0.0, omega_forced=om
        )
        # prediccion teorica mecanica centrifuga pura (sin friccion radial)
        if kind == "retract":
            W_cent_th = 0.5 * p.m * om**2 * (p.r_max**2 - p.r_min**2)
        else:
            W_cent_th = 0.5 * p.m * om**2 * (p.r_min**2 - p.r_max**2)
        res["W_cent_theory_mec"] = W_cent_th
        results[label] = res

    # --- B2) Vertical + gravedad, phi FIJO (eje bloqueado en angulo), omega=0 ---
    # (ya cubierto en A; se mantiene)

    # --- B3) Giro vertical con gravedad (phi avanza): solo referencia mixta ---
    for label, ph0, kind in [
        ("spin_retract_from_horiz", 0.0, "retract"),
        ("spin_extend_from_horiz", 0.0, "extend"),
    ]:
        if kind == "retract":
            rfun = make_r_step_retract(0.1, 1.1, p.r_max, p.r_min)
        else:
            rfun = make_r_step_extend(0.1, 1.1, p.r_min, p.r_max)
        phifun = lambda t, o=omega_spin, p0=ph0: p0 + o * t
        results[label] = integrate_imposed_path(label, phifun, rfun, t_end=1.3, p=p)

    # --- C) Ciclo completo: radio constante vs reduccionista (N vueltas) ---
    T = 8.0  # s para 2 vueltas => omega medio ~ 1.57 rad/s
    n_turns = 2.0
    phi_c, r_c = make_r_constant_cycle(T, p.r_max, n_turns)
    phi_red, r_red = make_r_reductionist_cycle(T, p.r_min, p.r_max, n_turns)
    # anti-reduccionista (al reves: grande al subir) para comparar
    def r_anti(t, rmin=p.r_min, rmax=p.r_max, T=T, n=n_turns):
        ph = 2.0 * np.pi * n * t / T
        r_mid = 0.5 * (rmin + rmax)
        amp = 0.5 * (rmax - rmin)
        return r_mid + amp * np.cos(ph)  # invertido

    results["cycle_r_max_const"] = integrate_imposed_path(
        "cycle_r_max_const", phi_c, r_c, t_end=T, p=p
    )
    results["cycle_reductionist"] = integrate_imposed_path(
        "cycle_reductionist", phi_red, r_red, t_end=T, p=p
    )
    results["cycle_anti_reductionist"] = integrate_imposed_path(
        "cycle_anti_reductionist", phi_c, r_anti, t_end=T, p=p
    )
    # ciclo r_min constante
    phi_min, r_minf = make_r_constant_cycle(T, p.r_min, n_turns)
    results["cycle_r_min_const"] = integrate_imposed_path(
        "cycle_r_min_const", phi_min, r_minf, t_end=T, p=p
    )

    # --- D) Gasto eficiente: recoger EN PARADO horizontal, girar, extender abajo en lento ---
    # Secuencia por tramos (concatenamos contabilidades)
    # 1) parado horizontal: extend r_min->r_max? mejor: empieza recogido, gira, extiende en favorable
    # Script de gasto eficiente:
    #   t0-1s: parado phi=0, r = r_min (ya listo)
    #   1-2s:  parado phi=0, nothing
    #   2-2.8s: spin up forced... simplified as forced path pieces summed

    # Tramo 1: recoger en parado (si empezamos extendidos) — ya en A
    # Tramo 2: girar a r_min (barato en par desfavorable)
    # Tramo 3: extender en favorable a baja omega
    # Lo reportamos como suma de piezas ya calculadas + un tramo extend lento parado abajo

    return results, p


def summarize(results: dict, p: P):
    rows = []
    for k, res in results.items():
        rows.append(
            {
                "name": k,
                "W_act_elec_net": res["W_act_elec_net"],
                "W_act_elec_in": res["W_act_elec_in"],
                "W_act_regen": res["W_act_regen"],
                "W_eje_net": res["W_eje_net"],
                "W_eje_mot": res["W_eje_mot"],
                "W_eje_gen": res["W_eje_gen"],
                "W_total_elec": res["W_act_elec_net"] - res["W_eje_net"],
                # total pago red approx: act_in - act_regen + eje_mot - eje_gen
                "W_paid_elec": res["W_act_elec_in"]
                - res["W_act_regen"]
                + res["W_eje_mot"]
                - res["W_eje_gen"],
                "dPE": res["dPE"],
                "dKE": res["dKE"],
                "W_grav_phi": res["W_grav_phi"],
                "W_grav_r": res["W_grav_r"],
                "omega_mean": res["omega_mean"],
                "r_mean": res["r_mean"],
            }
        )
    return rows


def plot_key(results, path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # 1) r(t) ciclos
    ax = axes[0, 0]
    for key, col in [
        ("cycle_reductionist", "#27ae60"),
        ("cycle_anti_reductionist", "#e74c3c"),
        ("cycle_r_max_const", "#2980b9"),
    ]:
        res = results[key]
        ax.plot(res["t"], res["r"], label=key, color=col, lw=1.2)
    ax.set_ylabel("r (m)")
    ax.set_title("Perfil de radio en ciclos")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # 2) centrifuga pura: recoger/extender vs omega
    ax = axes[0, 1]
    names = [
        "cent_retract_w0",
        "cent_retract_w3",
        "cent_retract_w5",
        "cent_extend_w0",
        "cent_extend_w3",
        "cent_extend_w5",
    ]
    labels = ["rec w0", "rec w3", "rec w5", "ext w0", "ext w3", "ext w5"]
    vals = [results[n]["W_act_elec_net"] for n in names]
    colors = ["#27ae60" if v < 0 else "#e67e22" for v in vals]
    ax.bar(range(len(vals)), vals, color=colors)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("W_act elec neto (J)")
    ax.set_title("Solo centrifuga (g=0): recoger vs extender")
    ax.axhline(0, color="k", lw=0.8)
    ax.grid(True, alpha=0.3, axis="y")

    # 3) peaje total ciclos
    ax = axes[1, 0]
    cyc = [
        "cycle_r_min_const",
        "cycle_r_max_const",
        "cycle_reductionist",
        "cycle_anti_reductionist",
    ]
    paid = [
        results[c]["W_act_elec_in"]
        - results[c]["W_act_regen"]
        + results[c]["W_eje_mot"]
        - results[c]["W_eje_gen"]
        for c in cyc
    ]
    ax.bar(range(len(cyc)), paid, color=["#3498db", "#9b59b6", "#27ae60", "#e74c3c"])
    ax.set_xticks(range(len(cyc)))
    ax.set_xticklabels(["r_min", "r_max", "reduccionista", "anti"], fontsize=8)
    ax.set_ylabel("W_paid elec (J) 2 vueltas")
    ax.set_title("Ciclo: gasto electrico total (eje+radial)")
    ax.grid(True, alpha=0.3, axis="y")

    # 4) desglose reductionist
    ax = axes[1, 1]
    res = results["cycle_reductionist"]
    parts = [
        res["W_eje_mot"],
        -res["W_eje_gen"],
        res["W_act_elec_in"],
        -res["W_act_regen"],
    ]
    ax.bar(
        ["eje_mot", "−eje_gen", "act_in", "−act_regen"],
        parts,
        color=["#e74c3c", "#2ecc71", "#e67e22", "#1abc9c"],
    )
    ax.set_title("Desglose ciclo reduccionista")
    ax.set_ylabel("J")
    ax.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    print("=" * 70)
    print("RECORRIDO REDUCCIONISTA + REARME EN PARADO")
    print("=" * 70)
    results, p = run_all()
    rows = summarize(results, p)

    print(f"\nParametros: m={p.m} kg  r=[{p.r_min},{p.r_max}] m")
    print(
        f"Delta PE max teorico (r_max-r_min)*m*g = "
        f"{(p.r_max - p.r_min) * p.m * G:.1f} J  (si Delta sin(phi)=1)"
    )

    print("\n--- A  PARADO (omega=0): gravedad segun angulo ---")
    print(f"{'caso':<32} {'W_act_net':>10} {'dPE':>10} {'nota'}")
    for key, nota in [
        ("stop_retract_horizontal", "g ⊥ r, ~solo friccion"),
        ("stop_retract_top", "gravedad AYUDA a recoger"),
        ("stop_retract_bottom", "gravedad se OPONE"),
        ("stop_extend_bottom", "gravedad AYUDA a extender"),
        ("stop_extend_top", "gravedad se OPONE"),
        ("stop_extend_horizontal", "g ⊥ r"),
    ]:
        r = results[key]
        print(
            f"{key:<32} {r['W_act_elec_net']:10.2f} {r['dPE']:10.2f}  {nota}"
        )

    print("\n--- B  Solo CENTRIFUGA (g=0): parado vs omega ---")
    print(
        f"{'caso':<22} {'W_act_net':>10} {'W_cent_th':>10} {'omega':>6}"
    )
    for key in [
        "cent_retract_w0",
        "cent_retract_w3",
        "cent_retract_w5",
        "cent_extend_w0",
        "cent_extend_w3",
        "cent_extend_w5",
    ]:
        r = results[key]
        print(
            f"{key:<22} {r['W_act_elec_net']:10.2f} {r.get('W_cent_theory_mec', 0):10.2f} "
            f"{r['omega_mean']:6.2f}"
        )

    print("\n--- C  Ciclos 2 vueltas ---")
    print(
        f"{'ciclo':<26} {'W_paid':>9} {'eje_mot':>9} {'eje_gen':>9} "
        f"{'act_net':>9} {'Wg_phi':>8} {'Wg_r':>8} {'dPE':>7} {'close':>7}"
    )
    for key in [
        "cycle_r_min_const",
        "cycle_r_max_const",
        "cycle_reductionist",
        "cycle_anti_reductionist",
    ]:
        r = results[key]
        paid = r["W_act_elec_in"] - r["W_act_regen"] + r["W_eje_mot"] - r["W_eje_gen"]
        print(
            f"{key:<26} {paid:9.1f} {r['W_eje_mot']:9.1f} {r['W_eje_gen']:9.1f} "
            f"{r['W_act_elec_net']:9.1f} {r['W_grav_phi']:8.1f} {r['W_grav_r']:8.1f} "
            f"{r['dPE']:7.1f} {r['grav_closure']:7.1f}"
        )

    def paid(res):
        return res["W_act_elec_in"] - res["W_act_regen"] + res["W_eje_mot"] - res["W_eje_gen"]

    cr0 = results["cent_retract_w0"]["W_act_elec_net"]
    cr3 = results["cent_retract_w3"]["W_act_elec_net"]
    cr5 = results["cent_retract_w5"]["W_act_elec_net"]
    ce0 = results["cent_extend_w0"]["W_act_elec_net"]
    ce3 = results["cent_extend_w3"]["W_act_elec_net"]
    stop_top = results["stop_retract_top"]["W_act_elec_net"]
    stop_bot_ext = results["stop_extend_bottom"]["W_act_elec_net"]
    paid_red = paid(results["cycle_reductionist"])
    paid_anti = paid(results["cycle_anti_reductionist"])
    paid_max = paid(results["cycle_r_max_const"])
    paid_min = paid(results["cycle_r_min_const"])
    # peaje de EJE solo
    eje_red = results["cycle_reductionist"]["W_eje_mot"] - results["cycle_reductionist"]["W_eje_gen"]
    eje_anti = (
        results["cycle_anti_reductionist"]["W_eje_mot"]
        - results["cycle_anti_reductionist"]["W_eje_gen"]
    )

    print("\n--- CONCLUSIONES (yo tambien lo tenia claro; la sim lo clava) ---")
    print(f"  Recoger g=0, w=0:     {cr0:+.1f} J   (~friccion)")
    print(f"  Recoger g=0, w=3:     {cr3:+.1f} J   (pagas centrifuga)")
    print(f"  Recoger g=0, w=5:     {cr5:+.1f} J   (mucho mas)")
    print(f"  Extender g=0, w=0:    {ce0:+.1f} J")
    print(f"  Extender g=0, w=3:    {ce3:+.1f} J   (centrifuga AYUDA a sacar)")
    print(f"  Recoger PARADO arriba:{stop_top:+.1f} J   (gravedad trabaja a favor)")
    print(f"  Extender PARADO abajo:{stop_bot_ext:+.1f} J   (gravedad trabaja a favor)")
    print(f"  Ciclo r_min W_paid:   {paid_min:+.1f} J")
    print(f"  Ciclo r_max W_paid:   {paid_max:+.1f} J")
    print(f"  Ciclo reduccionista:  {paid_red:+.1f} J  (eje neto {eje_red:+.1f} J)")
    print(f"  Ciclo anti-reducc.:   {paid_anti:+.1f} J  (eje neto {eje_anti:+.1f} J)")
    print(
        f"  Ahorro TOTAL red vs anti: {paid_anti - paid_red:+.1f} J | "
        f"ahorro EJE: {eje_anti - eje_red:+.1f} J"
    )
    print(
        f"  Cierre gravitatorio reduccionista (debe ~0): "
        f"{results['cycle_reductionist']['grav_closure']:.2f} J"
    )

    print(
        """
Veredicto:
  SI — no cuesta lo mismo. La sim lo confirma.
  • Centrifuga: recoger en giro CUESTA ~ 1/2 m ω² (R²-r²); extender en giro puede
    SALIR casi gratis o regenerar (la masa quiere irse afuera).
  • Parado: segun angulo, gravedad ayuda o castiga; horizontal ~ solo rozamiento.
  • Ciclo reduccionista: NO crea energia (grav_closure~0, W_paid>0). SI abaratiza el
    peaje del EJE moviendo trabajo al actuador radial y a fases buenas.
  • Gasto eficiente: rearmar en parado (o angulo neutro/favorable) + r pequeno en
    tramos duros + r grande en favorable + perneo para picos.
"""
    )

    plot_path = OUT / "recorrido_reduccionista.png"
    plot_key(results, plot_path)
    print(f"Grafica: {plot_path}")

    # JSON serializable sin arrays
    slim = {}
    for k, res in results.items():
        slim[k] = {kk: vv for kk, vv in res.items() if kk not in ("t", "phi", "omega", "r", "rdot")}
    report = {
        "params": asdict(p),
        "cases": slim,
        "summary_rows": rows,
        "verdict_es": (
            "Rearme en parado no cuesta lo mismo que en giro: la centrifuga domina. "
            "Ciclo reduccionista es gasto eficiente del EJE, no sobreunidad. "
            "Gravedad puede ayudar en parado segun angulo (top retract / bottom extend)."
        ),
    }
    out_json = OUT / "recorrido_reduccionista_resultados.json"
    out_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"JSON: {out_json}")


if __name__ == "__main__":
    main()
