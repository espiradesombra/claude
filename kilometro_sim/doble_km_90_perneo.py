"""
Doble Kilometro desfasado 90° + Perneo/Desperneo
================================================
Objetivo de ingenieria (NO sobreunidad):
  - Suavizar par / potencia electrica neta (buffer mecanico)
  - Regeneracion SOLO en fase gravitatoria favorable
  - En fase desfavorable: motor o libre (regen DESACTIVADA)
  - Perneo: Km2 desacoplado (trinquete inercial)
  - Desperneo: embrague acopla Km2 (absorbe pico / aporta inercia)

Metricas:
  - ripple par:  (tau_max - tau_min) / max(|tau_mean|, eps)
  - P_pico / P_media  (potencia electrica neta al bus)
  - eta_paid = W_gen / W_motor  (<1 en ciclo cerrado)
  - 1a ley residual

Perfiles de carga:
  - base: potencia media de red
  - pulso: demanda extra (kW, duracion) que el buffer debe absorber/aportar
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Literal

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
G = 9.81


# ---------------------------------------------------------------------------
# Parametros
# ---------------------------------------------------------------------------
@dataclass
class Phys:
    # Escala de banco de pruebas (par gravitatorio ~ m g R gestionable por actuadores)
    # m=40, R=1.2 => m g R ≈ 471 N·m  (pico por modulo)
    m: float = 40.0  # kg por modulo
    R: float = 1.2  # m
    I_rotor: float = 25.0  # kg·m² por eje de modulo (sin masa)
    beta: float = 3.0  # N·m·s/rad friccion por modulo
    eta_gen: float = 0.92
    eta_mot: float = 0.90


@dataclass
class Control:
    omega_ref: float = 2.0  # rad/s
    # Motor en zona DESFAVORABLE (regen OFF ahi)
    # bias ~ 0.55 * m g R para vencer el punto duro sin regenerar en contra
    Kp_motor: float = 400.0
    tau_motor_bias: float = 280.0  # N·m
    # Regeneracion SOLO en zona FAVORABLE
    c_gen: float = 35.0  # N·m·s/rad
    free_band: float = 0.12  # rad/s


@dataclass
class LoadPulse:
    """Demanda electrica positiva = red pide potencia (el sistema genera).
    Negativa = red inyecta / sobra (el sistema absorbe como freno).
    """
    P_base_W: float = 200.0  # W media pedida al generador (si puede)
    P_pulse_W: float = 2500.0  # W pico extra (2.5 kW)
    t_pulse_start: float = 8.0  # s
    t_pulse_dur: float = 0.8  # s descarga
    # Segundo pulso de absorcion (exceso de red / frenado)
    P_absorb_W: float = -1800.0
    t_absorb_start: float = 14.0
    t_absorb_dur: float = 0.6


@dataclass
class PerneoPolicy:
    enabled: bool = True
    # Perneo: Km2 bloqueado en reposo relativo (embrague abierto, I_extra=0 al eje)
    # Desperneo si:
    #   a) omega del activo > omega_hi (riesgo de embalamiento), o
    #   b) demanda de pico |P_load| > P_desperneo
    omega_hi: float = 2.6
    omega_lo: float = 1.2  # re-pernear si se frena demasiado tras acoplar
    P_desperneo_W: float = 1500.0
    # Tras desperneo, mantener acoplado al menos t_min
    t_hold: float = 1.2


ModeName = Literal["mono", "dual_rigid", "dual_perneo"]


# ---------------------------------------------------------------------------
# Fisica de fase
# ---------------------------------------------------------------------------
def tau_grav(phi: float, m: float, R: float) -> float:
    """Par gravitatorio. Convención v2: tau_g = -m g R cos(phi).
    Con omega>0: favorable (asiste) cuando tau_g > 0 => cos(phi) < 0.
    """
    return -m * G * R * np.cos(phi)


def is_favorable(phi: float) -> bool:
    """Gravedad asiste el sentido +omega."""
    return np.cos(phi) < 0.0


def load_power(t: float, load: LoadPulse) -> float:
    P = load.P_base_W
    if load.t_pulse_start <= t < load.t_pulse_start + load.t_pulse_dur:
        P += load.P_pulse_W
    if load.t_absorb_start <= t < load.t_absorb_start + load.t_absorb_dur:
        P += load.P_absorb_W
    return P


# ---------------------------------------------------------------------------
# Control por modulo: REGEN SOLO FAVORABLE
# ---------------------------------------------------------------------------
def control_module(
    phi: float,
    omega: float,
    ctrl: Control,
    phys: Phys,
    P_cmd_share: float,
) -> tuple[float, str, float, float]:
    """Devuelve (tau_c, mode, P_elec_motor_pos, P_elec_gen_pos).

    REGLA USUARIO:
      - Fase desfavorable: regeneracion DESACTIVADA.
        Solo motor (si hace falta) o libre.
      - Fase favorable: regeneracion permitida (y puede seguir carga).
    """
    fav = is_favorable(phi)
    err = ctrl.omega_ref - omega
    tau_g = tau_grav(phi, phys.m, phys.R)

    P_mot = 0.0
    P_gen = 0.0

    if not fav:
        # === DESFAVORABLE: SIN REGEN (ni freno electrico) ===
        # Feedforward cancela gravedad adversa + PI sobre omega
        if err > -ctrl.free_band:
            tau_ff = max(0.0, -tau_g)  # solo compensa si gravedad se opone
            tau_c = tau_ff + ctrl.Kp_motor * max(err, 0.0) + ctrl.tau_motor_bias * 0.25
            # Nunca negativo en desfavorable => no regenera
            tau_c = max(0.0, tau_c)
            mode = "motor_unfav"
            P_mot = max(tau_c * omega, 0.0) / phys.eta_mot
        else:
            # sobra velocidad: libre (NO regen en desfavorable)
            tau_c = 0.0
            mode = "free_unfav"
        return tau_c, mode, P_mot, P_gen

    # === FAVORABLE: REGEN ON (con techo: no matar omega) ===
    # Inercia sintetica: solo se cosecha EXCESO de velocidad / par gravitatorio.
    # No se persigue un P_cmd de kW si eso hunde el rotor (buffer, no planta).

    # 1) Cosecha por exceso de omega (primario)
    tau_speed = 0.0
    if omega > ctrl.omega_ref:
        tau_speed = -ctrl.c_gen * (omega - ctrl.omega_ref)

    # 2) Contribucion a carga: fraccion limitada del par gravitatorio favorable
    tau_load = 0.0
    if omega > ctrl.omega_ref * 0.85 and P_cmd_share > 0.0:
        tau_cmd = -P_cmd_share / (omega * phys.eta_gen)
        # techo: no freenes mas que ~70% del par gravitatorio asistente + holgura
        tau_cap = -0.70 * max(tau_g, 0.0) - 0.30 * ctrl.c_gen * max(omega - ctrl.omega_ref, 0.0)
        # tau_cmd y tau_cap son <= 0; el mas cercano a 0 es el menos agresivo
        tau_load = max(tau_cmd, tau_cap)  # max porque ambos negativos
    elif P_cmd_share < 0.0 and omega < ctrl.omega_ref * 1.15:
        # Red sobra energia: absorber como motor mecanico (carga el buffer)
        tau_load = min(
            (-P_cmd_share) * phys.eta_mot / max(omega, 0.05),
            ctrl.tau_motor_bias + ctrl.Kp_motor * max(err, 0.0),
        )

    tau_c = tau_speed + tau_load

    # 3) Proteccion de velocidad: si omega cae, corta regen (libera el eje)
    if omega < ctrl.omega_ref * 0.75 and tau_c < 0.0:
        tau_c = 0.0
        mode = "free_fav_protect"
        return tau_c, mode, P_mot, P_gen

    # 4) Rescue motor solo si deficit fuerte (raro en favorable)
    if err > 2 * ctrl.free_band and tau_c >= -1e-9:
        tau_c = ctrl.Kp_motor * err * 0.4
        mode = "motor_fav_rescue"
        P_mot = max(tau_c * omega, 0.0) / phys.eta_mot
        return tau_c, mode, P_mot, P_gen

    if tau_c < 0:
        mode = "gen_fav"
        P_gen = max(-tau_c * omega, 0.0) * phys.eta_gen
    elif tau_c > 0:
        mode = "motor_fav_absorb"
        P_mot = max(tau_c * omega, 0.0) / phys.eta_mot
    else:
        mode = "free_fav"

    return tau_c, mode, P_mot, P_gen


# ---------------------------------------------------------------------------
# Simulacion
# ---------------------------------------------------------------------------
@dataclass
class SimResult:
    mode: ModeName
    t: np.ndarray
    phi1: np.ndarray
    phi2: np.ndarray
    omega1: np.ndarray
    omega2: np.ndarray  # si perneo abierto, puede diferir; si rigid, igual
    omega_shaft: np.ndarray
    tau_g_total: np.ndarray
    tau_c_total: np.ndarray
    P_elec_net: np.ndarray  # + genera a bus, - consume del bus
    P_load: np.ndarray
    P_bus_error: np.ndarray  # P_elec_net - P_load (ideal 0)
    coupled: np.ndarray  # 1 si Km2 acoplado
    W_motor: float
    W_gen: float
    W_fric: float
    W_grav: float
    dKE: float
    metrics: dict = field(default_factory=dict)


def run_sim(
    mode: ModeName,
    phys: Phys | None = None,
    ctrl: Control | None = None,
    load: LoadPulse | None = None,
    perneo: PerneoPolicy | None = None,
    dt: float = 0.0005,
    t_max: float = 20.0,
) -> SimResult:
    phys = phys or Phys()
    ctrl = ctrl or Control()
    load = load or LoadPulse()
    perneo = perneo or PerneoPolicy(enabled=(mode == "dual_perneo"))

    n = int(t_max / dt)
    t = np.zeros(n)
    phi1 = np.zeros(n)
    phi2 = np.zeros(n)
    omega1 = np.zeros(n)
    omega2 = np.zeros(n)
    omega_shaft = np.zeros(n)
    tau_g_total = np.zeros(n)
    tau_c_total = np.zeros(n)
    P_elec_net = np.zeros(n)
    P_load_arr = np.zeros(n)
    P_bus_error = np.zeros(n)
    coupled_arr = np.zeros(n)

    I_m = phys.m * phys.R**2
    I1 = phys.I_rotor + I_m
    I2 = phys.I_rotor + I_m

    # Estado inicial
    p1 = 0.0
    # Dual 90°: desfase pi/2
    p2 = np.pi / 2 if mode != "mono" else 0.0
    o1 = ctrl.omega_ref
    o2 = ctrl.omega_ref if mode != "dual_perneo" else 0.0  # perneado: Km2 quieto
    coupled = mode != "dual_perneo"  # mono/rigid: siempre "acoplado" conceptual
    if mode == "mono":
        coupled = True
    t_coupled_since = -1e9

    W_motor = 0.0
    W_gen = 0.0
    W_fric = 0.0
    W_grav = 0.0
    KE0 = 0.5 * I1 * o1**2 + (0.5 * I2 * o2**2 if mode != "mono" else 0.0)

    for i in range(n):
        ti = i * dt
        t[i] = ti
        P_ld = load_power(ti, load)
        P_load_arr[i] = P_ld

        # --- logica perneo ---
        if mode == "dual_perneo" and perneo.enabled:
            in_hard_pulse = abs(P_ld) > perneo.P_desperneo_W
            if coupled:
                # No re-pernear en medio de un pulso duro (mantener I alta)
                if (
                    (not in_hard_pulse)
                    and (ti - t_coupled_since) > perneo.t_hold
                    and o1 < perneo.omega_lo
                ):
                    coupled = False
                    o2 = 0.0
            else:
                need = (o1 > perneo.omega_hi) or in_hard_pulse
                if need:
                    # Desperneo inelastico: conservacion de momento angular
                    # I1*o1 + I2*o2 = (I1+I2)*o_common  (o2~0 si estaba perneado)
                    o_common = (I1 * o1 + I2 * o2) / (I1 + I2)
                    o1 = o_common
                    o2 = o_common
                    p2 = p1 + np.pi / 2  # enganche a 90°
                    coupled = True
                    t_coupled_since = ti

        # --- torques gravitatorios ---
        tg1 = tau_grav(p1, phys.m, phys.R)
        if mode == "mono":
            tg2 = 0.0
        else:
            if mode == "dual_rigid" or coupled:
                p2 = p1 + np.pi / 2  # cigüeñal rígido 90°
            tg2 = tau_grav(p2, phys.m, phys.R)

        # reparto de consignas de carga entre modulos favorables activos
        fav1 = is_favorable(p1)
        km2_active = mode != "mono" and (mode == "dual_rigid" or coupled)
        fav2 = is_favorable(p2) if km2_active else False
        n_fav = int(fav1) + int(fav2)
        if n_fav == 0:
            share1 = share2 = 0.0
        elif n_fav == 1:
            share1 = P_ld if fav1 else 0.0
            share2 = P_ld if fav2 else 0.0
        else:
            share1 = share2 = 0.5 * P_ld

        # control
        tc1, _mode1, Pm1, Pg1 = control_module(p1, o1, ctrl, phys, share1)
        if mode == "mono":
            tc2, Pm2, Pg2 = 0.0, 0.0, 0.0
        elif km2_active:
            tc2, _mode2, Pm2, Pg2 = control_module(p2, o1, ctrl, phys, share2)
        else:
            # Km2 perneado (quieto): sin actuacion
            tc2, Pm2, Pg2 = 0.0, 0.0, 0.0

        tf1 = -phys.beta * o1
        if mode == "mono":
            tf2 = 0.0
        elif km2_active:
            tf2 = -phys.beta * o1
        else:
            tf2 = 0.0

        # --- dinamica (Euler: trabajos con omega actual, luego integrar) ---
        if mode == "mono":
            I_eff = I1
            tau_tot = tg1 + tc1 + tf1
            alpha = tau_tot / I_eff
            W_grav += tg1 * o1 * dt
            W_fric += abs(tf1 * o1 * dt)
            W_motor += Pm1 * dt
            W_gen += Pg1 * dt
            p1 = p1 + o1 * dt
            o1 = max(0.0, o1 + alpha * dt)
            o2 = o1
            p2 = p1
            o_shaft = o1
            tau_g_sum = tg1
            tau_c_sum = tc1
            P_net = Pg1 - Pm1
        elif mode == "dual_rigid" or coupled:
            I_eff = I1 + I2
            tau_tot = tg1 + tg2 + tc1 + tc2 + tf1 + tf2
            alpha = tau_tot / I_eff
            W_grav += (tg1 + tg2) * o1 * dt
            W_fric += abs((tf1 + tf2) * o1 * dt)
            W_motor += (Pm1 + Pm2) * dt
            W_gen += (Pg1 + Pg2) * dt
            p1 = p1 + o1 * dt
            o1 = max(0.0, o1 + alpha * dt)
            o2 = o1
            # cigüeñal rígido a 90°: desfase cinematico fijo
            p2 = p1 + np.pi / 2
            o_shaft = o1
            tau_g_sum = tg1 + tg2
            tau_c_sum = tc1 + tc2
            P_net = (Pg1 + Pg2) - (Pm1 + Pm2)
        else:
            # dual_perneo abierto: solo Km1 gira; Km2 bloqueado (phi fijo)
            I_eff = I1
            tau_tot = tg1 + tc1 + tf1
            alpha = tau_tot / I_eff
            W_grav += tg1 * o1 * dt
            W_fric += abs(tf1 * o1 * dt)
            W_motor += Pm1 * dt
            W_gen += Pg1 * dt
            p1 = p1 + o1 * dt
            o1 = max(0.0, o1 + alpha * dt)
            o2 = 0.0
            o_shaft = o1
            tau_g_sum = tg1
            tau_c_sum = tc1
            P_net = Pg1 - Pm1

        # log
        phi1[i] = p1
        phi2[i] = p2
        omega1[i] = o1
        omega2[i] = o2
        omega_shaft[i] = o_shaft
        tau_g_total[i] = tau_g_sum
        tau_c_total[i] = tau_c_sum
        P_elec_net[i] = P_net
        P_bus_error[i] = P_net - P_ld
        coupled_arr[i] = 1.0 if (mode != "dual_perneo" or coupled) else 0.0

    # recortar por si se usara break en el futuro (ahora corre t_max completo)
    KE1 = 0.5 * I1 * omega1[-1] ** 2
    if mode != "mono":
        KE1 += 0.5 * I2 * omega2[-1] ** 2
    dKE = KE1 - KE0

    metrics = compute_metrics(
        t, tau_g_total, P_elec_net, P_load_arr, W_motor, W_gen, W_fric, W_grav, dKE
    )
    metrics["mode"] = mode
    metrics["regen_unfavorable"] = False  # garantizado por control
    metrics["note"] = "regen SOLO en fase favorable (cos(phi)<0)"
    metrics["omega_end"] = float(omega_shaft[-1]) if len(omega_shaft) else 0.0
    metrics["omega_min"] = float(np.min(omega_shaft)) if len(omega_shaft) else 0.0
    metrics["omega_max"] = float(np.max(omega_shaft)) if len(omega_shaft) else 0.0
    metrics["mgR"] = float(phys.m * G * phys.R)

    # Metricas en ventana del pulso de descarga (utilidad buffer)
    pulse_mask = (t >= load.t_pulse_start) & (
        t < load.t_pulse_start + load.t_pulse_dur
    )
    if np.any(pulse_mask):
        metrics["omega_min_during_pulse"] = float(np.min(omega_shaft[pulse_mask]))
        metrics["P_err_rms_during_pulse_W"] = float(
            np.sqrt(np.mean((P_elec_net[pulse_mask] - P_load_arr[pulse_mask]) ** 2))
        )
        metrics["fraction_coupled_during_pulse"] = float(np.mean(coupled_arr[pulse_mask]))
        # Cobertura del pulso: energia electrica entregada / energia pedida
        E_req = float(np.trapz(np.maximum(P_load_arr[pulse_mask], 0.0), t[pulse_mask]))
        E_del = float(np.trapz(np.maximum(P_elec_net[pulse_mask], 0.0), t[pulse_mask]))
        metrics["E_pulse_requested_J"] = E_req
        metrics["E_pulse_delivered_J"] = E_del
        metrics["pulse_coverage"] = float(E_del / E_req) if E_req > 1e-9 else float("nan")
    else:
        metrics["omega_min_during_pulse"] = float("nan")
        metrics["P_err_rms_during_pulse_W"] = float("nan")
        metrics["fraction_coupled_during_pulse"] = float("nan")
        metrics["E_pulse_requested_J"] = float("nan")
        metrics["E_pulse_delivered_J"] = float("nan")
        metrics["pulse_coverage"] = float("nan")

    # Base load (sin pulsos): ondulación electrica
    base_mask = (t >= 2.0) & (t < load.t_pulse_start - 0.5)
    if np.any(base_mask):
        pb = P_elec_net[base_mask]
        metrics["P_base_std_W"] = float(np.std(pb))
        metrics["P_base_pp_W"] = float(np.max(pb) - np.min(pb))
    else:
        metrics["P_base_std_W"] = float("nan")
        metrics["P_base_pp_W"] = float("nan")

    return SimResult(
        mode=mode,
        t=t,
        phi1=phi1,
        phi2=phi2,
        omega1=omega1,
        omega2=omega2,
        omega_shaft=omega_shaft,
        tau_g_total=tau_g_total,
        tau_c_total=tau_c_total,
        P_elec_net=P_elec_net,
        P_load=P_load_arr,
        P_bus_error=P_bus_error,
        coupled=coupled_arr,
        W_motor=W_motor,
        W_gen=W_gen,
        W_fric=W_fric,
        W_grav=W_grav,
        dKE=dKE,
        metrics=metrics,
    )


def compute_metrics(t, tau_g, P_net, P_load, W_motor, W_gen, W_fric, W_grav, dKE):
    # tramo estable (descartar 1 s inicial si hay datos)
    if len(t) == 0:
        return {
            "tau_g_mean": 0.0,
            "tau_g_pp": 0.0,
            "ripple_tau_pp_over_rms": float("nan"),
            "P_mean_W": 0.0,
            "P_peak_abs_W": 0.0,
            "P_peak_over_P_mean": float("nan"),
            "P_track_err_rms_W": float("nan"),
            "W_motor_J": W_motor,
            "W_gen_J": W_gen,
            "W_fric_J": W_fric,
            "W_grav_J": W_grav,
            "dKE_J": dKE,
            "W_net_elec_J": W_gen - W_motor,
            "eta_paid": float("nan"),
            "omega_end": 0.0,
        }
    mask = t >= (t[0] + 1.0) if t[-1] > 2.0 else np.ones_like(t, dtype=bool)
    if not np.any(mask):
        mask = np.ones_like(t, dtype=bool)
    tg = tau_g[mask]
    pn = P_net[mask]
    pl = P_load[mask]
    tg_mean = float(np.mean(tg))
    tg_pp = float(np.max(tg) - np.min(tg))
    ripple_tau_pp_rms = tg_pp / max(float(np.sqrt(np.mean(tg**2))), 1e-6)

    P_mean = float(np.mean(pn))
    P_peak = float(np.max(np.abs(pn)))
    ratio = P_peak / max(abs(P_mean), 1e-6)
    err_rms = float(np.sqrt(np.mean((pn - pl) ** 2)))
    eta_paid = float(W_gen / W_motor) if W_motor > 1e-9 else float("nan")
    return {
        "tau_g_mean": tg_mean,
        "tau_g_pp": tg_pp,
        "ripple_tau_pp_over_rms": ripple_tau_pp_rms,
        "P_mean_W": P_mean,
        "P_peak_abs_W": P_peak,
        "P_peak_over_P_mean": ratio,
        "P_track_err_rms_W": err_rms,
        "W_motor_J": W_motor,
        "W_gen_J": W_gen,
        "W_fric_J": W_fric,
        "W_grav_J": W_grav,
        "dKE_J": dKE,
        "W_net_elec_J": W_gen - W_motor,
        "eta_paid": eta_paid,
    }


# ---------------------------------------------------------------------------
# Comparativa + plots
# ---------------------------------------------------------------------------
def plot_comparison(results: dict[str, SimResult], path: Path):
    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    colors = {"mono": "#e74c3c", "dual_rigid": "#2980b9", "dual_perneo": "#27ae60"}

    ax = axes[0]
    for name, r in results.items():
        ax.plot(r.t, r.tau_g_total, label=name, color=colors.get(name), lw=1.0)
    ax.set_ylabel("τ_grav total (N·m)")
    ax.set_title("Doble Km 90° + Perneo — regen SOLO en fase favorable")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for name, r in results.items():
        ax.plot(r.t, r.P_elec_net / 1000.0, label=name, color=colors.get(name), lw=1.0)
    # carga del primer result
    r0 = next(iter(results.values()))
    ax.plot(r0.t, r0.P_load / 1000.0, "k--", lw=1.2, label="P_load", alpha=0.7)
    ax.set_ylabel("P (kW)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)

    ax = axes[2]
    for name, r in results.items():
        ax.plot(r.t, r.omega_shaft, label=name, color=colors.get(name), lw=1.0)
    ax.set_ylabel("ω eje (rad/s)")
    ax.grid(True, alpha=0.3)

    ax = axes[3]
    if "dual_perneo" in results:
        r = results["dual_perneo"]
        ax.plot(r.t, r.coupled, color=colors["dual_perneo"], lw=1.5)
        ax.set_ylabel("perneo: 0=abierto 1=desperneado")
    else:
        ax.text(0.1, 0.5, "sin dual_perneo", transform=ax.transAxes)
    ax.set_xlabel("t (s)")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    print("=" * 68)
    print("DOBLE KILOMETRO 90° + PERNEO  |  regen SOLO fase favorable")
    print("=" * 68)
    print("Regla de control:")
    print("  cos(phi) >= 0  -> DESFAVORABLE: motor o libre; REGEN OFF")
    print("  cos(phi) <  0  -> FAVORABLE:    regeneracion ON (+ carga)")
    print()

    phys = Phys()
    ctrl = Control()
    # Perfil de pico: ~2.5 kW durante 0.8 s (ingenieria de buffer, no motor perpetuo)
    load = LoadPulse(
        P_base_W=80.0,  # W base pedida (buffer, no planta de generacion)
        P_pulse_W=2500.0,  # 2.5 kW pico / 0.8 s
        t_pulse_start=8.0,
        t_pulse_dur=0.8,
        P_absorb_W=-1800.0,  # -1.8 kW absorcion / 0.6 s
        t_absorb_start=14.0,
        t_absorb_dur=0.6,
    )
    perneo = PerneoPolicy(
        enabled=True,
        omega_hi=2.3,  # despernea si se embala
        omega_lo=1.0,  # re-pernea solo si ya bajo mucho (fuera de pulso)
        P_desperneo_W=800.0,  # ante demanda alta acopla Km2 (mas I)
        t_hold=1.5,
    )

    modes: list[ModeName] = ["mono", "dual_rigid", "dual_perneo"]
    results: dict[str, SimResult] = {}
    for m in modes:
        print(f"Simulando {m}...")
        results[m] = run_sim(
            m, phys=phys, ctrl=ctrl, load=load, perneo=perneo, dt=0.0005, t_max=20.0
        )

    # tabla
    print()
    print(
        f"{'modo':<14} {'τg_pp':>8} {'ωmin_p':>8} {'cobertura':>9} "
        f"{'E_del/E_req':>14} {'η_paid':>7} {'W_net':>9} {'ω_end':>7}"
    )
    print("-" * 90)
    table = {}
    for name, r in results.items():
        met = r.metrics
        print(
            f"{name:<14} {met['tau_g_pp']:8.1f} {met['omega_min_during_pulse']:8.3f} "
            f"{met['pulse_coverage']:9.2%} "
            f"{met['E_pulse_delivered_J']:.0f}/{met['E_pulse_requested_J']:.0f}J".rjust(14)
            + f" {met['eta_paid']:7.3f} {met['W_net_elec_J']:9.1f} {met['omega_end']:7.3f}"
        )
        table[name] = met
    print(f"\nm g R (por modulo) = {phys.m * G * phys.R:.1f} N·m")
    print(
        f"Pulso: {load.P_pulse_W/1000:.2f} kW x {load.t_pulse_dur:.2f} s "
        f"= {load.P_pulse_W * load.t_pulse_dur:.0f} J"
    )
    print(
        "Nota 90°: τ_g dual = -m g R √2 cos(φ+π/4) → τ_pp crece ×√2 vs mono; "
        "la ganancia es complementariedad motor/gen y perneo ante picos, "
        "no cancelacion del armónico (eso requiere 3@120°)."
    )

    plot_path = OUT / "doble_km_90_perneo.png"
    plot_comparison(results, plot_path)
    print(f"\nGrafica: {plot_path}")

    # --- demo explicita: un modulo nunca regenera en desfavorable ---
    # muestreo de modos en mono
    print("\nComprobacion regen desfavorable (muestreo mono, 2 s):")
    # re-sim corta instrumentada
    verify_regen_policy(phys, ctrl)

    report = {
        "phys": asdict(phys),
        "control": asdict(ctrl),
        "load_pulse": asdict(load),
        "perneo": asdict(perneo),
        "metrics": table,
        "policy": {
            "regen_in_unfavorable": False,
            "favorable_when": "cos(phi) < 0  (tau_g > 0 con omega>0)",
            "unfavorable_actions": ["motor", "free"],
            "favorable_actions": ["gen", "free", "motor_rescue"],
        },
        "engineering_verdict_es": (
            "Politica: regeneracion SOLO en fase favorable (cos phi < 0); "
            "en desfavorable solo motor o libre. "
            "Doble 90° no anula el 1er armonico de tau_g (amplitud ×√2); "
            "sirve para desfasar zonas motor/gen y para perneo/desperneo "
            "como supercondensador mecanico ante picos de demanda. "
            "eta_paid < 1 siempre: no es planta de generacion neta. "
            "Quijote 'gigante viable' vs molino terrestre: el cuello es "
            "par_aero / par_actuador y tiempo de ciclo, no magia gravitatoria. "
            "Kilometro aporta buffer e inercia sintetica, no sustituye al viento."
        ),
    }
    out_json = OUT / "doble_km_90_perneo_resultados.json"
    out_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"JSON: {out_json}")
    print("\n" + report["engineering_verdict_es"])


def verify_regen_policy(phys: Phys, ctrl: Control):
    """Barre un ciclo y cuenta intentos de regen en desfavorable (debe ser 0)."""
    n_unfav = n_unfav_gen = n_fav = n_fav_gen = 0
    for phi in np.linspace(0, 2 * np.pi, 720, endpoint=False):
        for omega in (0.5, 1.8, 2.5):
            tc, mode, Pm, Pg = control_module(phi, omega, ctrl, phys, P_cmd_share=500.0)
            if is_favorable(phi):
                n_fav += 1
                if Pg > 0 or "gen" in mode:
                    n_fav_gen += 1
            else:
                n_unfav += 1
                if Pg > 0 or "gen" in mode:
                    n_unfav_gen += 1
    print(f"  muestras desfavorables: {n_unfav}, con regen: {n_unfav_gen}  (debe ser 0)")
    print(f"  muestras favorables:    {n_fav}, con regen: {n_fav_gen}")
    assert n_unfav_gen == 0, "FALLO: regeneracion en fase desfavorable"
    print("  OK: regeneracion desactivada en desfavorable.")


if __name__ == "__main__":
    main()
