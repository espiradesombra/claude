import numpy as np

# --- PARAMETROS DEL SISTEMA (MODULO KILOMETRO) --- Gemini original
m = 50.0
g = 9.81
r = 1.0
eta_gen = 0.92
eta_mot = 0.90
dt = 0.001
t_max = 5.0
omega_target = 0.5
F_patada = 1200.0
t_patada_dur = 0.2
F_avance_gen = 300.0
dist_avance = 2.0

theta = 0.0
omega = 0.1
x_avance = 0.0
v_avance = 0.0
E_gen_bajada = 0.0
E_cons_patada = 0.0
E_gen_avance = 0.0
time = 0.0
fase = "A_BAJADA"
t_patada_start = 0.0
steps = 0

print("--- INICIANDO SIMULACION NUMERICA (script Gemini original) ---")

while time < t_max:
    if fase == "A_BAJADA":
        tau_grav = m * g * r * np.sin(theta)
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
            P_mot = (F_patada * v_avance) / eta_mot if v_avance > 0 else (F_patada * 1.0) / eta_mot
            E_cons_patada += P_mot * dt
            a_patada = (F_patada - m * g * 0.5) / m
            v_avance += a_patada * dt
            x_avance += v_avance * dt
        else:
            fase = "C_AVANCE_REGENERATIVO"
    elif fase == "C_AVANCE_REGENERATIVO":
        F_grav_fav = m * g * 0.8
        if F_grav_fav > F_avance_gen:
            a_avance = (F_grav_fav - F_avance_gen) / m
        else:
            a_avance = 0.0
        v_avance += a_avance * dt
        x_avance += v_avance * dt
        P_gen_lin = F_avance_gen * v_avance * eta_gen
        E_gen_avance += P_gen_lin * dt
        if x_avance >= dist_avance:
            break
    time += dt
    steps += 1

E_total_generada = E_gen_bajada + E_gen_avance
E_total_consumida = E_cons_patada
E_neto = E_total_generada - E_total_consumida

print(f"\nResultados del Ciclo Completo (Gemini original):")
print(f"  - Tiempo / pasos:                    {time:.3f} s / {steps}")
print(f"  - Fase final:                        {fase}")
print(f"  - theta final:                       {theta:.3f} rad")
print(f"  - omega final:                       {omega:.3f} rad/s")
print(f"  - x_avance / v_avance:                {x_avance:.3f} m / {v_avance:.3f} m/s")
print(f"  - Energia Generada en Bajada:        {E_gen_bajada:.2f} J")
print(f"  - Energia Consumida en Patada:       {E_cons_patada:.2f} J")
print(f"  - Energia Generada en Avance:        {E_gen_avance:.2f} J")
print(f"  --------------------------------------------------")
print(f"  - ENERGIA TOTAL GENERADA:            {E_total_generada:.2f} J")
print(f"  - BALANCE NETO POR CICLO:            {E_neto:.2f} J")
if E_neto > 0:
    print("\n [¡EXITO!] El ciclo produce un EXCEDENTE neto de energia (segun Gemini).")
else:
    print("\n [AJUSTE REQUERIDO] La patada consume mas de lo que compensa el ciclo.")
