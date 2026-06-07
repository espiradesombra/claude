"""
KILÒMETRE — Estudi de dimensions i bateria en cicles continus

Preguntes:
  1. Quina energia emmagatzema el sistema? (kJ, kWh?)
  2. Quina potència de pic pot donar?
  3. Quant tarda a carregar/descarregar?
  4. Quins paràmetres físics importen més?

Anàlisi dimensional:
  Energia cinètica: E_k = ½·I·ω²
  I = I_tub + I_objecte = ½·m_tub·R²_tub + m·R²
  → E_k ∝ m·R²·ω²

  Energia potencial disponible (si l'objecte baixa h):
  E_p = m·g·h

  Parell màxim generador via sinfín:
  τ_gen = k_gen·ω  → P_gen = k_gen·ω²

  Temps de descàrrega: t_desc = E_k / P_gen = I/(2·k_gen)
  → independent de ω! (com una bateria RC)

Variables de disseny:
  m   [kg]  → massa objecte
  R   [m]   → radi
  ω   [rad/s] → velocitat angular
  L   [m]   → longitud tub
  pas [m/rad] → pas sinfín
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

print("=" * 65)
print("KILÒMETRE — Estudi de dimensions")
print("Bateria cinètica via sinfín + tub rotant")
print("=" * 65)

# ── Paràmetres base ─────────────────────────────────────────
g    = 9.81   # m/s²

# Escombrat de dimensions
masses   = [1, 5, 10, 50, 100]       # kg
radis    = [0.05, 0.10, 0.15, 0.30, 0.50]  # m
omeges   = [1, 3, 6, 10, 20]         # rad/s
k_gens   = [0.5, 1.0, 2.0, 5.0]     # N·m·s/rad

print("\n=== TAULA: Energia cinètica E_k = ½·m·R²·ω² ===")
print(f"{'m(kg)':>6} | {'R(m)':>5} | {'ω(rad/s)':>9} | "
      f"{'I(kg·m²)':>9} | {'E_k(J)':>8} | {'E_k(Wh)':>8} | {'P_max(W)':>9}")
print("-" * 70)

casos_base = [
    # (m, R, omega, k_gen)  — casos representatius
    (10,   0.15, 3,  2.0),   # prototip petit
    (50,   0.30, 6,  5.0),   # mig
    (100,  0.50, 10, 10.0),  # gran
    (500,  1.00, 6,  20.0),  # industrial
    (1000, 1.50, 3,  50.0),  # molt gran
]

for m, R, om, k_gen in casos_base:
    m_tub  = m * 0.3   # tub = 30% massa objecte
    R_tub  = R * 1.2
    I      = 0.5*m_tub*R_tub**2 + m*R**2
    E_k    = 0.5 * I * om**2
    P_max  = k_gen * om**2   # potència pico generador
    t_desc = E_k / P_max     # temps descàrrega (s)
    print(f"{m:>6} | {R:>5.2f} | {om:>9.1f} | "
          f"{I:>9.3f} | {E_k:>8.1f} | {E_k/3600:>8.5f} | {P_max:>9.1f}  "
          f"t_desc={t_desc:.1f}s")

print()
print("=== ANÀLISI DIMENSIONAL: quins paràmetres importen? ===")
print("""
E_k = ½·(m_tub·R²_tub/2 + m·R²)·ω²  ≈  ½·m·R²·ω²

Dependències:
  E_k ∝ m    (lineal en massa)          → doblar m dobla E_k
  E_k ∝ R²   (quadràtica en radi)       → doblar R quadruplica E_k  ★
  E_k ∝ ω²   (quadràtica en velocitat)  → doblar ω quadruplica E_k  ★

Per maximitzar la bateria: R i ω són els paràmetres clau.

Potència pic: P = k_gen·ω²
  → també ∝ ω²: a més velocitat, més potència instantània

Temps de descàrrega: t = E_k/P = I/(2·k_gen)
  → independent de ω: sempre el mateix temps!
  → per augmentar t_desc cal augmentar I o reduir k_gen

Comparació amb bateria de liti:
  Li-ion: ~150-250 Wh/kg  (densitat energia alta, lenta)
  Volant cinètic: ~5-100 Wh/kg  (densitat menor, molt ràpid)
  Kilòmetre: depèn de R i ω
""")

# ── Cicle continu: càrrega i descàrrega ─────────────────────
print("=" * 65)
print("CICLE CONTINU: càrrega (motor) i descàrrega (gen)")
print("=" * 65)

# Paràmetres del cas base
m      = 10.0;  R = 0.15
m_tub  = 3.0;   R_tub = 0.18
I      = 0.5*m_tub*R_tub**2 + m*R**2
beta   = 0.02
eta_gen = 0.90;  eta_mot = 0.92

omega_max = 10.0   # rad/s  (carregat)
omega_min = 1.0    # rad/s  (descarregat)
k_gen_desc = 3.0   # N·m·s/rad  (frenada descàrrega)
tau_motor  = I * 3.0  # N·m  (acceleració càrrega)

dt = 0.001
t  = 0.0

# Fase 1: CÀRREGA (motor accelera de omega_min a omega_max)
print("\n--- Fase CÀRREGA ---")
om = omega_min
W_el_in_carrega = 0.0
t_carrega = 0.0
while om < omega_max:
    tau_fr = -beta * om
    alpha  = (tau_motor + tau_fr) / I
    om_new = om + alpha * dt
    dW_mot = tau_motor * om * dt
    W_el_in_carrega += max(dW_mot, 0) / eta_mot
    om = min(om_new, omega_max)
    t += dt
    t_carrega += dt

E_k_max = 0.5 * I * omega_max**2
E_k_min = 0.5 * I * omega_min**2
dE_k    = E_k_max - E_k_min

print(f"  Temps càrrega:     {t_carrega:.2f} s")
print(f"  ΔE_k emmagatzemat: {dE_k:.3f} J")
print(f"  W_elèctric consumit: {W_el_in_carrega:.3f} J")
print(f"  Eficiència càrrega: {100*dE_k/W_el_in_carrega:.1f}%")

# Fase 2: DESCÀRREGA (generador de omega_max a omega_min)
print("\n--- Fase DESCÀRREGA ---")
om = omega_max
W_el_out_desc = 0.0
t_desc_real   = 0.0
while om > omega_min:
    tau_fr  = -beta * om
    tau_gen = -k_gen_desc * om
    alpha   = (tau_gen + tau_fr) / I
    om_new  = om + alpha * dt
    dW_gen  = abs(tau_gen * om * dt)
    W_el_out_desc += dW_gen * eta_gen
    om = max(om_new, omega_min)
    t += dt
    t_desc_real += dt

print(f"  Temps descàrrega:   {t_desc_real:.2f} s")
print(f"  W_elèctric produït: {W_el_out_desc:.3f} J")
print(f"  Eficiència descàrrega: {100*W_el_out_desc/dE_k:.1f}%")
print(f"  Eficiència TOTAL (round-trip): "
      f"{100*W_el_out_desc/W_el_in_carrega:.1f}%")

# ── Simulació completa: N cicles ─────────────────────────────
print("\n=== SIMULACIÓ: 10 cicles càrrega/descàrrega ===")
dt   = 0.0001
N_ci = 10

t_arr   = []
om_arr  = []
mode_arr= []
Wei_arr = []
Weo_arr = []

om  = omega_min
Wei = 0.0; Weo = 0.0
t   = 0.0

for cicle in range(N_ci):
    # CÀRREGA
    while om < omega_max - 0.01:
        tau_fr = -beta * om
        alpha  = (tau_motor + tau_fr) / I
        om    += alpha * dt
        om     = min(om, omega_max)
        dWm    = max(tau_motor * om * dt, 0) / eta_mot
        Wei   += dWm
        t     += dt
        t_arr.append(t); om_arr.append(om)
        mode_arr.append(1); Wei_arr.append(Wei); Weo_arr.append(Weo)

    # DESCÀRREGA
    while om > omega_min + 0.01:
        tau_fr  = -beta * om
        tau_gen = -k_gen_desc * om
        alpha   = (tau_gen + tau_fr) / I
        om     += alpha * dt
        om      = max(om, omega_min)
        dWg     = abs(tau_gen * om * dt) * eta_gen
        Weo    += dWg
        t      += dt
        t_arr.append(t); om_arr.append(om)
        mode_arr.append(2); Wei_arr.append(Wei); Weo_arr.append(Weo)

t_arr    = np.array(t_arr)
om_arr   = np.array(om_arr)
mode_arr = np.array(mode_arr)
Wei_arr  = np.array(Wei_arr)
Weo_arr  = np.array(Weo_arr)

print(f"  Temps total: {t:.2f} s")
print(f"  W_el consumit total: {Wei:.3f} J")
print(f"  W_el produït total:  {Weo:.3f} J")
print(f"  Eficiència round-trip: {100*Weo/Wei:.1f}%")
print(f"  Energia per cicle: {Weo/N_ci:.3f} J  |  "
      f"Potència mitja: {Weo/t:.3f} W")

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 12))
fig.suptitle("KILÒMETRE — Bateria cinètica: dimensions i cicles\n"
             f"m={m}kg  R={R}m  ω∈[{omega_min},{omega_max}] rad/s  "
             f"η_total={100*Weo/Wei:.0f}%",
             fontsize=12, fontweight='bold')

gs = gridspec.GridSpec(3, 2, figure=fig)

# 1. E_k vs R per diverses masses
ax1 = fig.add_subplot(gs[0, 0])
Rs = np.linspace(0.05, 1.0, 100)
for m_plot in [10, 50, 100, 500]:
    I_plot = m_plot * Rs**2  # simplificat
    E_plot = 0.5 * I_plot * omega_max**2
    ax1.plot(Rs, E_plot, lw=2, label=f'm={m_plot}kg')
ax1.set_xlabel('Radi R [m]')
ax1.set_ylabel('E_k màxima [J]')
ax1.set_title(f'Energia vs radi  (ω={omega_max} rad/s)')
ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)

# 2. E_k vs ω per diversos radis
ax2 = fig.add_subplot(gs[0, 1])
oms = np.linspace(1, 30, 100)
for R_plot in [0.10, 0.20, 0.50, 1.00]:
    I_plot = m * R_plot**2
    E_plot = 0.5 * I_plot * oms**2
    ax2.plot(oms, E_plot, lw=2, label=f'R={R_plot}m')
ax2.set_xlabel('Velocitat angular ω [rad/s]')
ax2.set_ylabel('E_k màxima [J]')
ax2.set_title(f'Energia vs velocitat  (m={m}kg)')
ax2.legend(fontsize=9); ax2.grid(True, alpha=0.3)

# 3. Cicles: velocitat angular
ax3 = fig.add_subplot(gs[1, :])
mask_m = mode_arr == 1
mask_g = mode_arr == 2
ax3.fill_between(t_arr, om_arr, omega_min,
                 where=mask_m, color='red', alpha=0.25, label='Càrrega (motor)')
ax3.fill_between(t_arr, om_arr, omega_min,
                 where=mask_g, color='green', alpha=0.25, label='Descàrrega (gen)')
ax3.plot(t_arr, om_arr, 'steelblue', lw=1.5, label='ω(t)')
ax3.axhline(omega_max, color='red', lw=1, ls='--', alpha=0.5)
ax3.axhline(omega_min, color='green', lw=1, ls='--', alpha=0.5)
ax3.set_ylabel('ω [rad/s]'); ax3.set_xlabel('Temps [s]')
ax3.set_title(f'{N_ci} cicles càrrega/descàrrega')
ax3.legend(fontsize=9); ax3.grid(True, alpha=0.3)

# 4. Energies acumulades
ax4 = fig.add_subplot(gs[2, 0])
ax4.plot(t_arr, Wei_arr, 'red',   lw=1.5, label='W_el_in (consumit)')
ax4.plot(t_arr, Weo_arr, 'green', lw=1.5, label='W_el_out (produït)')
ax4.set_xlabel('Temps [s]'); ax4.set_ylabel('Energia [J]')
ax4.set_title('Energies acumulades')
ax4.legend(fontsize=9); ax4.grid(True, alpha=0.3)

# 5. Eficiència acumulada
ax5 = fig.add_subplot(gs[2, 1])
eta_cum = np.where(Wei_arr > 0.001, Weo_arr/Wei_arr*100, 0)
ax5.plot(t_arr, eta_cum, 'purple', lw=1.5)
ax5.axhline(eta_gen*eta_mot*100, color='gray', lw=1, ls='--',
            label=f'η_teòric={eta_gen*eta_mot*100:.0f}%')
ax5.set_ylim(0, 100)
ax5.set_xlabel('Temps [s]'); ax5.set_ylabel('Eficiència [%]')
ax5.set_title('Eficiència round-trip acumulada')
ax5.legend(fontsize=9); ax5.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/kilometre_dimensions.png', dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_dimensions.py  |  .png")
