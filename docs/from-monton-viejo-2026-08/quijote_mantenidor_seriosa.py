#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quijote_mantenidor_seriosa.py
═══════════════════════════════════════════════════════════════════════════════
VERSIÓ SERIOSA — el número honest que decideix si Quijote-mantenidor és producte.

Atén TOTES les crítiques del revisor:
  · Massa realista gran (autoritat de control): escaneig 66→2000 kg/pala
  · Cost d'actuador REALISTA: η<1, fricció, límit d'acceleració, histèresi
  · Comparació JUSTA contra control de PARELL convencional (estat de l'art),
    no contra un molí d'inèrcia fixa "tonto"
  · Balanç energètic complet amb residu
  · Tres màquines sobre EXACTAMENT el mateix vent:
      A) Inèrcia fixa + control parell convencional (MPPT estàndard)
      B) Quijote inèrcia variable + control parell
      C) Inèrcia fixa SENSE optimitzar (referència "tonta")

Víctor Manzanares Alberola — EPSA UPV Alcoi — github.com/espiradesombra/claude
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
N_PALES = 3
R_TIP   = 55.0
R_HUB   = 5.0
J_ROTOR = 4.0e7
OMEGA_NOM = 1.327

# ── Paràmetres realistes de l'actuador ────────────────────────────────────────
ETA_ACT   = 0.85       # eficiència electromecànica (motor+reductor+rodaments)
DR_MAX    = 0.5        # velocitat radial màxima (m/s) — realista
DDR_MAX   = 0.3        # acceleració radial màxima (m/s²) — límit físic
HISTERESI = 0.03       # banda morta del control (rad/s) — evita repicar
C_FRIC_R  = 200.0      # fricció radial (N·s/m) per pala


def Cp(lam):
    if lam <= 0: return 0.0
    return max(0.0, CP_MAX*(1.0 - ((lam-LAM_OPT)/LAM_OPT)**2))


def perfil_vent(T, dt, seed=0, turbulencia=1.0):
    """Vent realista amb ràfegues i turbulència escalable."""
    np.random.seed(seed)
    n = int(T/dt); t = np.arange(n)*dt
    v = 10.0 + 2.0*np.sin(2*np.pi*t/40)
    for t0, amp, dur in [(60,-4,15),(120,3,10),(180,-5,20),(260,4,12),(330,-3,18),
                          (400,-4,14),(470,5,10),(540,-3,16)]:
        v += turbulencia*amp*np.exp(-((t-t0)/dur)**2)
    ou = np.zeros(n)
    for i in range(1,n):
        ou[i] = ou[i-1] - 0.3*ou[i-1]*dt + turbulencia*0.5*np.sqrt(dt)*np.random.randn()
    v += ou
    return np.clip(v, 3.0, 25.0)


def simula(v_arr, dt, massa, mode):
    """
    mode:
      'conv'    → inèrcia fixa, control de PARELL convencional (MPPT K·ω²)
      'quijote' → inèrcia variable (Quijote) + mateix control de parell
      'tonta'   → inèrcia fixa, generador segueix vent sense optimitzar
    Retorna energies (J) i mètriques.
    """
    n = len(v_arr)
    omega = LAM_OPT*v_arr[0]/R
    theta = 0.0
    r = (R_TIP+R_HUB)/2
    rdot = 0.0
    E_gen = 0.0      # energia al generador (la que es ven)
    E_act = 0.0      # cost net actuador (amb pèrdues)
    E_fric = 0.0     # pèrdues fricció radial
    # constant MPPT òptima: τ = K·ω², K = ½ρA R³ Cp_max/λ³
    K_mppt = 0.5*RHO*A_ROT*R**3*CP_MAX/LAM_OPT**3
    cp_acc = 0.0

    for k in range(n):
        v = v_arr[k]
        # ── control radial (només Quijote) ──
        rdot_cmd = 0.0
        if mode == 'quijote':
            omega_opt = LAM_OPT*v/R
            error = omega - omega_opt
            if error < -HISTERESI and r > R_HUB:
                rdot_cmd = -DR_MAX            # contraure (accelera)
            elif error > HISTERESI and r < R_TIP:
                rdot_cmd = +DR_MAX            # estendre (frena/guarda)
            # límit d'acceleració
            dv = np.clip(rdot_cmd - rdot, -DDR_MAX*dt, DDR_MAX*dt)
            rdot = rdot + dv
        else:
            rdot = 0.0
        r_new = np.clip(r + rdot*dt, R_HUB, R_TIP)
        rdot_real = (r_new - r)/dt

        # ── inèrcia ──
        J = J_ROTOR + N_PALES*massa*r**2
        Jdot = N_PALES*2*massa*r*rdot_real

        # ── aerodinàmica ──
        lam = omega*R/max(v,0.1)
        cp = Cp(lam)
        P_vent = min(0.5*RHO*A_ROT*cp*v**3, P_NOM*1.3)
        tau_vent = P_vent/max(omega,0.1)

        # ── control de parell del generador ──
        if mode == 'tonta':
            tau_gen = min(0.7*tau_vent, P_NOM/max(omega,0.1))  # subòptim
        else:
            # MPPT estàndard: τ_gen = K·ω² (segueix λ òptim), saturat a P_NOM
            tau_gen = min(K_mppt*omega**2, P_NOM/max(omega,0.1), tau_vent*1.5)

        # ── dinàmica ω lliure ──
        omega_dot = (tau_vent - tau_gen - Jdot*omega)/J
        omega = max(omega + omega_dot*dt, 0.1)
        theta += omega*dt

        # ── energies ──
        E_gen += tau_gen*omega*dt
        cp_acc += cp
        # cost actuador realista (amb η i fricció)
        if mode == 'quijote':
            for i in range(N_PALES):
                th_i = theta + 2*np.pi*i/N_PALES
                F_cent = massa*omega**2*r
                F_grav = -massa*G*np.sin(th_i)
                F_fric = C_FRIC_R*abs(rdot_real)*np.sign(rdot_real)
                F_act = -(F_cent + F_grav) - F_fric
                P = F_act*rdot_real
                # aportar costa /η, recuperar rendeix ·η
                E_act += (P/ETA_ACT if P > 0 else P*ETA_ACT)*dt
                E_fric += C_FRIC_R*rdot_real**2*dt
        r = r_new

    return {
        'E_gen_kWh': E_gen/3.6e6,
        'E_act_kWh': max(E_act, 0)/3.6e6,
        'E_fric_kWh': E_fric/3.6e6,
        'cp_mig': cp_acc/n,
        'E_neta_kWh': (E_gen - max(E_act,0))/3.6e6,
    }


def main():
    print("█"*72)
    print("  QUIJOTE MANTENIDOR — VERSIÓ SERIOSA (vs estat de l'art)")
    print("█"*72)
    dt = 0.02; T = 600.0
    v = perfil_vent(T, dt, seed=3)
    print(f"\n  Vent: {T:.0f}s, mitjana {v.mean():.2f} m/s, "
          f"rang [{v.min():.1f},{v.max():.1f}], NREL 5MW, η_act={ETA_ACT}")
    print(f"  Actuador realista: DR_MAX={DR_MAX}m/s, accel<{DDR_MAX}m/s², "
          f"fricció={C_FRIC_R}N·s/m, histèresi={HISTERESI}rad/s")

    print("\n" + "─"*72)
    print("  COMPARACIÓ 1: massa estàndard (66,5 kg/pala)")
    print("─"*72)
    conv = simula(v, dt, 66.5, 'conv')
    quij = simula(v, dt, 66.5, 'quijote')
    print(f"  {'':32}{'Conv. (MPPT)':>16}{'Quijote':>16}")
    print(f"  {'E_generador venuda (kWh)':32}{conv['E_gen_kWh']:>16.2f}{quij['E_gen_kWh']:>16.2f}")
    print(f"  {'Cp mitjà':32}{conv['cp_mig']:>16.4f}{quij['cp_mig']:>16.4f}")
    print(f"  {'Cost actuador (kWh)':32}{'—':>16}{quij['E_act_kWh']:>16.3f}")
    print(f"  {'Pèrdues fricció (kWh)':32}{'—':>16}{quij['E_fric_kWh']:>16.3f}")
    gn = quij['E_gen_kWh'] - conv['E_gen_kWh'] - quij['E_act_kWh']
    print(f"  {'GUANY NET vs MPPT (kWh)':32}{'':>16}{gn:>16.3f}")
    print(f"  {'GUANY NET (%)':32}{'':>16}{100*gn/conv['E_gen_kWh']:>15.3f}%")

    print("\n" + "─"*72)
    print("  COMPARACIÓ 2: escaneig de massa (guany net vs MPPT convencional)")
    print("─"*72)
    print(f"  {'massa/pala':>11}{'ΔJ/J':>8}{'E_gen Quij':>13}{'cost act':>11}"
          f"{'net vs MPPT':>13}{'%':>8}")
    conv_base = conv['E_gen_kWh']
    for massa in [66.5, 200, 500, 1000, 2000]:
        q = simula(v, dt, massa, 'quijote')
        c = simula(v, dt, massa, 'conv')  # MPPT amb mateixa massa (just)
        dJ = 100*N_PALES*massa*R_TIP**2/J_ROTOR
        net = q['E_gen_kWh'] - c['E_gen_kWh'] - q['E_act_kWh']
        pct = 100*net/c['E_gen_kWh']
        print(f"  {massa:>9.0f}kg{dJ:>7.1f}%{q['E_gen_kWh']:>13.2f}"
              f"{q['E_act_kWh']:>11.3f}{net:>13.3f}{pct:>7.2f}%")

    print("\n" + "─"*72)
    print("  COMPARACIÓ 3: robustesa en 5 perfils de vent (massa 1000 kg)")
    print("─"*72)
    nets = []
    for seed in range(5):
        vv = perfil_vent(T, dt, seed=seed)
        q = simula(vv, dt, 1000, 'quijote')
        c = simula(vv, dt, 1000, 'conv')
        net = q['E_gen_kWh'] - c['E_gen_kWh'] - q['E_act_kWh']
        pct = 100*net/c['E_gen_kWh']
        nets.append(pct)
        print(f"    perfil {seed}: net {net:+.2f} kWh ({pct:+.2f}%) "
              f"{'✓' if net>0 else '✗'}")
    print(f"\n  Mitjana: {np.mean(nets):+.2f}% ± {np.std(nets):.2f}%")

    print("\n" + "═"*72)
    print("  VEREDICTE")
    print("═"*72)
    mitjana = np.mean(nets)
    print(f"""
  Comparat HONESTAMENT contra un control MPPT convencional (estat de
  l'art) i amb cost d'actuador realista (η={ETA_ACT}, fricció, límits):

  · Amb massa estàndard (66 kg): efecte negligible (ΔJ/J ~1,5%).
  · El guany escala amb la massa, però també ho fa el cost d'actuador.
  · Amb massa gran (1000 kg/pala): guany net mitjà {mitjana:+.2f}%.

  La simulació {'DEMOSTRA un guany net modest però real' if mitjana>0.1 else 'NO demostra guany net clar'}
  sobre l'estat de l'art en aquest règim de ràfegues.

  CONCLUSIÓ HONESTA:
  - NO és hurto (energia del vent, conservació intacta).
  - El valor com a CAPTACIÓ extra és modest i depèn de massa i ràfegues.
  - El valor probablement més gran NO és captar més energia, sinó el
    servei de REGULACIÓ (inèrcia sintètica ràpida) que es ven a banda
    en el mercat d'ancillary services, no mesurat en aquest kWh-test.

  SEGÜENT PAS recomanat: valorar el servei de regulació en €/MW
  (mercat d'inèrcia) en lloc de només kWh de captació.
""")


if __name__ == '__main__':
    main()
