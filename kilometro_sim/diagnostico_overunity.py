from check_gemini_y_fisica import ClosedParams, run_closed_cycle

cases = [
    ClosedParams(omega_ref=3.0, tau_motor_max=200.0, c_gen=200.0, beta=0.5, n_turns=5.0),
    ClosedParams(omega_ref=2.0, tau_motor_max=600.0, c_gen=80.0, beta=1.5, n_turns=5.0),
    ClosedParams(omega_ref=2.0, tau_motor_max=600.0, c_gen=80.0, beta=1.5, n_turns=20.0),
    ClosedParams(omega_ref=3.0, tau_motor_max=200.0, c_gen=200.0, beta=0.5, n_turns=20.0),
]

for p in cases:
    r = run_closed_cycle(p)
    I = p.I_rotor + p.m * p.R**2
    KE0 = 0.5 * I * p.omega_ref**2
    KE1 = 0.5 * I * r["omega_end"] ** 2
    true_surplus2 = r["W_gen_elec"] - r["W_motor_elec"] - p.eta_gen * max(0.0, -r["dKE"])
    Wm = r["W_motor_elec"] * p.eta_mot
    Wg = r["W_gen_elec"] / p.eta_gen
    print("---")
    print(
        f"omega={p.omega_ref} tau_m={p.tau_motor_max} c_gen={p.c_gen} "
        f"beta={p.beta} turns={p.n_turns}"
    )
    print(
        f"W_motor={r['W_motor_elec']:.1f} W_gen={r['W_gen_elec']:.1f} "
        f"W_net={r['W_net_elec']:.1f}"
    )
    print(
        f"W_fric={r['W_fric']:.1f} W_grav={r['W_grav']:.1f} "
        f"dPE={r['dPE']:.1f} dKE={r['dKE']:.1f}"
    )
    print(f"KE0={KE0:.1f} KE1={KE1:.1f} omega_end={r['omega_end']:.3f}")
    print(f"eta_paid={r['eta_paid']:.3f} resid%={r['mech_residual_pct']:.4f}")
    print(f"surplus_after_eta*KEdrain={true_surplus2:.1f}")
    print(
        f"check dKE+dPE={r['dKE']+r['dPE']:.2f}  "
        f"Wm-Wg-Wfric={Wm - Wg - r['W_fric']:.2f}"
    )
