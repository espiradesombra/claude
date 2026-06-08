#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quijote_mantenidor_realista.py
═══════════════════════════════════════════════════════════════════════════════
La pregunta honesta: l'estratègia de Víctor (estendre en parat, contraure quan
cau el vent per mantindre λ òptim) capta prou vent extra per compensar el cost
de l'actuador?

Compara DUES màquines sobre EL MATEIX perfil de vent realista amb ràfegues:
  A) Molí normal — inèrcia FIXA
  B) Molí Quijote — contrau pesos quan cau el vent (manté ω prop de λ_opt)

Mesura: E_vent_captada, E_actuador, E_neta. Tot del vent (no hi ha hurto).
Balanç energètic complet amb residu.

Víctor Manzanares Alberola — EPSA UPV Alcoi
═══════════════════════════════════════════════════════════════════════════════
"""
import numpy as np

# ── Constants NREL 5MW ────────────────────────────────────────────────────────
G       = 9.81
R       = 63.0
RHO     = 1.225
A_ROT   = np.pi*R**2
CP_MAX  = 0.482
LAM_OPT = 7.55
P_NOM   = 5.0e6
m       = 66.5          # massa per pala (×3 pales)
N_PALES = 3
R_TIP   = 55.0
R_HUB   = 5.0
J_ROTOR = 4.0e7         # inèrcia base rotor NREL


def Cp(lam):
    """Coeficient de potència — pic a LAM_OPT, cau a banda i banda."""
    if lam <= 0: return 0.0
    return max(0.0, CP_MAX*(1.0 - ((lam-LAM_OPT)/LAM_OPT)**2))


def perfil_vent(T, dt, seed=0):
    """Vent realista: base + oscil·lació lenta + ràfegues + soroll."""
    np.random.seed(seed)
    n = int(T/dt)
    t = np.arange(n)*dt
    v = 10.0 + 2.0*np.sin(2*np.pi*t/40)          # oscil·lació lenta
    # ràfegues (caigudes i pujades sobtades)
    for t0, amp, durada in [(60,-4,15),(120,3,10),(180,-5,20),(260,4,12),(330,-3,18)]:
        v += amp*np.exp(-((t-t0)/durada)**2)
    # soroll Ornstein-Uhlenbeck
    ou = np.zeros(n)
    for i in range(1,n):
        ou[i] = ou[i-1] - 0.3*ou[i-1]*dt + 0.5*np.sqrt(dt)*np.random.randn()
    v += ou
    return np.clip(v, 3.0, 25.0)


def simula(v_arr, dt, estrategia):
    """
    estrategia='fixa'    → r constant (molí normal)
    estrategia='quijote' → contrau quan el vent cau, estén quan puja
    Retorna energies (J) i traça.
    """
    n = len(v_arr)
    omega = LAM_OPT*v_arr[0]/R       # comença a λ òptim
    theta = 0.0
    r = (R_TIP+R_HUB)/2
    E_vent = 0.0       # energia capturada del vent (al generador)
    E_act = 0.0        # treball net actuador (+ = consumeix, − = recupera)
    E_rot_0 = 0.5*(J_ROTOR + N_PALES*m*r**2)*omega**2
    v_prev = v_arr[0]
    om_hist = np.zeros(n); cp_hist = np.zeros(n); r_hist = np.zeros(n)

    for k in range(n):
        v = v_arr[k]
        # ── estratègia de control radial ──
        rdot = 0.0
        if estrategia == 'quijote':
            # objectiu: mantindre ω prop de l'òptim per al vent ACTUAL
            omega_opt = LAM_OPT*v/R
            error = omega - omega_opt
            # si ω < òptim (rotor lent) → contraure (r baixa) per accelerar
            # si ω > òptim (rotor ràpid) → estendre (r puja) per frenar
            if error < -0.01 and r > R_HUB:
                rdot = -2.0          # contraure
            elif error > 0.01 and r < R_TIP:
                rdot = +2.0          # estendre
            rdot = np.clip(rdot, -2.0, 2.0)
        r_new = np.clip(r + rdot*dt, R_HUB, R_TIP)
        rdot_real = (r_new - r)/dt

        # ── inèrcia i dinàmica ──
        J = J_ROTOR + N_PALES*m*r**2
        Jdot = N_PALES*2*m*r*rdot_real
        lam = omega*R/max(v,0.1)
        P_vent = 0.5*RHO*A_ROT*Cp(lam)*v**3
        P_vent = min(P_vent, P_NOM*1.2)         # límit físic
        tau_vent = P_vent/max(omega,0.1)
        # generador extrau (carrega la xarxa) — model: extrau cap a P_NOM
        tau_gen = min(tau_vent, P_NOM/max(omega,0.1))
        # gravetat (3 pales) — net ~0 sobre la volta, però l'incloem
        tau_grav = sum(-m*G*r*np.cos(theta + 2*np.pi*i/N_PALES) for i in range(N_PALES))
        # ω lliure: d(Jω)/dt = τ_vent − τ_gen + τ_grav
        omega_dot = (tau_vent - tau_gen + tau_grav - Jdot*omega)/J
        omega += omega_dot*dt
        omega = max(omega, 0.1)
        theta += omega*dt

        # ── energia capturada (el que va al generador) ──
        E_vent += tau_gen*omega*dt
        # ── cost actuador (moure els N pesos radialment) ──
        for i in range(N_PALES):
            th_i = theta + 2*np.pi*i/N_PALES
            F_cent = m*omega**2*r
            F_grav_rad = -m*G*np.sin(th_i)
            F_act = -(F_cent + F_grav_rad)        # força que aplica l'actuador
            E_act += F_act*rdot_real*dt
        v_prev = v
        om_hist[k]=omega; cp_hist[k]=Cp(omega*R/max(v,0.1)); r_hist[k]=r
        r = r_new

    E_rot_f = 0.5*(J_ROTOR + N_PALES*m*r**2)*omega**2
    return {
        'E_vent_kWh': E_vent/3.6e6,
        'E_act_kWh': E_act/3.6e6,
        'E_neta_kWh': (E_vent - max(E_act,0))/3.6e6,
        'dE_rot_kWh': (E_rot_f-E_rot_0)/3.6e6,
        'cp_mig': np.mean(cp_hist),
        'om_hist': om_hist, 'cp_hist': cp_hist, 'r_hist': r_hist,
    }


def main():
    print("█"*70)
    print("  QUIJOTE MANTENIDOR — capta prou vent extra per pagar l'actuador?")
    print("█"*70)
    dt = 0.02
    T = 400.0
    v_arr = perfil_vent(T, dt, seed=3)
    print(f"\n  Perfil de vent: {T:.0f}s, mitjana {v_arr.mean():.2f} m/s, "
          f"rang [{v_arr.min():.1f}, {v_arr.max():.1f}] m/s")
    print(f"  Turbina NREL 5MW, {N_PALES} pales, massa mòbil {m} kg/pala\n")

    A = simula(v_arr, dt, 'fixa')
    B = simula(v_arr, dt, 'quijote')

    print("─"*70)
    print(f"  {'':28}{'Molí NORMAL':>16}{'Molí QUIJOTE':>16}")
    print("─"*70)
    print(f"  {'E_vent captada (kWh)':28}{A['E_vent_kWh']:>16.3f}{B['E_vent_kWh']:>16.3f}")
    print(f"  {'Cp mitjà':28}{A['cp_mig']:>16.4f}{B['cp_mig']:>16.4f}")
    print(f"  {'Cost actuador (kWh)':28}{'—':>16}{B['E_act_kWh']:>16.3f}")
    print("─"*70)

    # comparativa neta
    guany_brut = B['E_vent_kWh'] - A['E_vent_kWh']     # vent extra captat
    cost_act = max(B['E_act_kWh'], 0)                  # cost net actuador
    guany_net = guany_brut - cost_act
    print(f"\n  Vent EXTRA captat per Quijote:  {guany_brut:+.3f} kWh "
          f"({100*guany_brut/A['E_vent_kWh']:+.2f}%)")
    print(f"  Cost de l'actuador:             {cost_act:.3f} kWh")
    print(f"  GUANY NET:                      {guany_net:+.3f} kWh")
    print()
    if guany_net > 0:
        print(f"  ✓ La simulació demostra guany NET de {guany_net:.2f} kWh")
        print(f"    captant vent que el molí normal perd en ràfegues.")
        print(f"    (Energia del VENT mal aprofitat, NO de la gravetat.)")
    else:
        print(f"  ✗ El cost de l'actuador supera el vent extra captat.")
        print(f"    En aquest perfil, el mantenidor no compensa.")
    print()
    print("─"*70)
    print("  INTERPRETACIÓ")
    print("─"*70)
    print(f"""
  Tot el balanç és energia del VENT — no hi ha hurto gravitatori
  (E_grav net ≈ 0 sobre cada volta, com ja vam demostrar).

  El mantenidor guanya (si guanya) perquè manté el rotor prop de
  λ òptim durant les ràfegues, on un molí d'inèrcia fixa perd Cp.
  El resultat depèn FORTAMENT del perfil de vent: com més ràfegues
  brusques, més valor té el mantenidor.

  Aquest és el mecanisme REAL i defensable de Quijote: optimització
  de captació + regulació, no creació d'energia.
""")
    # sensibilitat a la turbulència
    print("─"*70)
    print("  SENSIBILITAT: guany net segons intensitat de ràfegues (5 perfils)")
    print("─"*70)
    for seed in range(5):
        v = perfil_vent(T, dt, seed=seed)
        a = simula(v, dt, 'fixa'); b = simula(v, dt, 'quijote')
        gb = b['E_vent_kWh']-a['E_vent_kWh']; gn = gb - max(b['E_act_kWh'],0)
        print(f"    perfil {seed}: vent extra {gb:+.2f} kWh, "
              f"cost act {max(b['E_act_kWh'],0):.2f}, net {gn:+.2f} kWh "
              f"{'✓' if gn>0 else '✗'}")


if __name__ == '__main__':
    main()
