"""
KILÒMETRE v8 — Flotació contínua (V_submergit variable)

Correcció del model v7: en lloc de flotació "tot o res" (escaló),
el volum submergit canvia contínuament quan l'objecte creua la superfície.

Geometria: objecte ESFÈRIC de radi r_obj
  - Centre de l'esfera a posició vertical: y_c = R·sin(φ) + h_eix
  - Volum submergit:
      y_c ≥  r_obj  → 0           (completament fora)
      y_c ≤ -r_obj  → V_total     (completament dins)
      -r_obj < y_c < r_obj → V_parcial (casquet esfèric)

  V_parcial(d) = π/3 · (r_obj - d)² · (2·r_obj + d)
    on d = y_c (distància centre a superfície, negatiu = dins)

Força de flotació: F_flot(φ) = ρ_agua · V_sub(φ) · g  [variable!]
Parell de flotació: τ_flot(φ) = F_flot(φ) · R · cos(φ)

Treball de flotació per volta:
  W_flot = ∮ τ_flot dφ
  Si W_flot ≠ 0 → font d'energia real (hidrostàtica)
  Si W_flot = 0 → la flotació és conservativa en aquest cicle

Verificació 1a Llei:
  ΔKE = 0 (omega constant)
  W_grav + W_flot + W_motor_mec + W_gen_mec - W_fric = 0
"""

import numpy as np
import matplotlib.pyplot as plt

# ════════════════════════════════════════════════════════════
#  PARÀMETRES
# ════════════════════════════════════════════════════════════
m        = 10.0    # kg  massa objecte
rho_obj  = 2000.0  # kg/m³
V_total  = m / rho_obj          # m³  volum total objecte
r_obj    = (3*V_total/(4*np.pi))**(1/3)  # radi esfera equivalent

rho_agua = 1000.0  # kg/m³
g        = 9.81    # m/s²
R        = 0.15    # m   radi tub (distància eix-objecte)

m_tub    = 2.0;  R_tub = 0.18
I        = 0.5*m_tub*R_tub**2 + m*R**2

beta     = 0.02    # N·m·s/rad
omega_ref = 3.0    # rad/s
alfa_gen  = 0.60
eta_gen   = 0.90
eta_mot   = 0.92

F_flot_max = rho_agua * V_total * g   # flotació màxima (completament submergit)

print(f"Objecte: m={m}kg  ρ={rho_obj}kg/m³  V={V_total*1e6:.2f}cm³  r={r_obj*100:.2f}cm")
print(f"Flotació màxima: {F_flot_max:.3f} N  ({100*F_flot_max/(m*g):.0f}% del pes)")
print(f"R_tub={R}m  |  r_obj={r_obj*100:.2f}cm  |  R/r_obj={R/r_obj:.1f}")
print()

# ════════════════════════════════════════════════════════════
#  FUNCIÓ: volum submergit d'una esfera
# ════════════════════════════════════════════════════════════
def V_submergit(y_centre, r):
    """
    Volum submergit d'una esfera de radi r amb centre a alçada y_centre.
    y_centre > 0  → centre per sobre l'aigua
    y_centre < 0  → centre per sota l'aigua
    """
    # d = y_centre (distància signada del centre a la superfície)
    d = y_centre
    if d >= r:
        return 0.0           # completament fora
    elif d <= -r:
        return 4/3*np.pi*r**3  # completament dins
    else:
        # Casquet esfèric: V = π/3·(r-d)²·(2r+d)
        return np.pi/3 * (r - d)**2 * (2*r + d)

V_sub_vec = np.vectorize(V_submergit)

# ════════════════════════════════════════════════════════════
#  SIMULACIÓ PER A UN h_eix DONAT
# ════════════════════════════════════════════════════════════
def simula(h_eix, alfa_gen=0.60, N=200000):
    phi  = np.linspace(0, 2*np.pi, N, endpoint=False)
    dphi = 2*np.pi / N
    om   = omega_ref

    # Alçada centre objecte respecte superfície
    y_c  = R * np.sin(phi) + h_eix

    # Volum submergit continu
    V_sub = V_sub_vec(y_c, r_obj)
    frac_sub_mig = np.mean(V_sub) / V_total   # fracció mitjana submergida

    # Forces i parells
    F_flot_v = rho_agua * V_sub * g            # variable
    tau_g    = -m * g * R * np.cos(phi)
    tau_f    = F_flot_v * R * np.cos(phi)      # parell flotació continu
    tau_fr   = -beta * om * np.ones(N)

    tau_gf   = tau_g + tau_f

    # Generador: actiu quan parell net (grav+flot) ajuda el gir
    tau_gen  = np.where(tau_gf > 0, -alfa_gen * tau_gf, 0.0)

    # Motor: compensa per mantindre omega constant
    tau_mot  = -(tau_gf + tau_gen + tau_fr)

    # Treball per element (regla única: τ·ω·dφ/ω = τ·dφ, amb factor ω pendent)
    # Nota: dW = τ·ω·dt = τ·ω·(dφ/ω) = τ·dφ  (per omega constant)
    dW_g   = tau_g   * dphi
    dW_f   = tau_f   * dphi
    dW_fr  = tau_fr  * dphi
    dW_gen = tau_gen * dphi
    dW_mot = tau_mot * dphi

    W_grav = np.sum(dW_g)
    W_flot = np.sum(dW_f)
    W_fric = np.sum(np.abs(dW_fr))

    # Verificació 1a Llei: W_grav+W_flot+W_mot+W_gen-W_fric = 0
    check  = W_grav + W_flot + np.sum(dW_mot) + np.sum(dW_gen) - W_fric
    residu = abs(check) / max(W_fric, 0.001) * 100

    # Elèctric
    W_el_gen = np.sum(np.where(dW_gen < 0, np.abs(dW_gen)*eta_gen, 0.0))
    W_el_in  = np.sum(np.where(dW_mot > 0, dW_mot/eta_mot, 0.0))
    W_el_rec = np.sum(np.where(dW_mot < 0, np.abs(dW_mot)*eta_gen, 0.0))
    W_el_out = W_el_gen + W_el_rec
    net_el   = W_el_out - W_el_in

    return dict(h_eix=h_eix, frac_sub=frac_sub_mig,
                W_grav=W_grav, W_flot=W_flot, W_fric=W_fric,
                W_el_in=W_el_in, W_el_out=W_el_out, net_el=net_el,
                residu=residu,
                phi=phi, tau_g=tau_g, tau_f=tau_f,
                tau_gen=tau_gen, tau_mot=tau_mot,
                V_sub=V_sub, y_c=y_c)

# ════════════════════════════════════════════════════════════
#  ESCOMBRAT % SUBMERGIT (30%→70%)
# ════════════════════════════════════════════════════════════
# Busquem h_eix per a cada fracció objectiu
# frac ≈ (π - arcsin(h_eix/R)) / (2π)  → h_eix = R·sin(π(1-2·frac))
fracs = np.array([0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70])
h_eixos = R * np.sin(np.pi * (1 - 2*fracs))

print("=" * 82)
print("ESCOMBRAT % SUBMERGIT — Flotació contínua (esfera)")
print(f"{'% sub':>6} | {'h_eix':>6} | {'W_grav':>9} | {'W_flot':>9} | "
      f"{'W_fric':>8} | {'el_in':>8} | {'el_out':>8} | {'net_el':>9} | {'res%':>6}")
print("-" * 82)

resultats = []
for frac, h in zip(fracs, h_eixos):
    r = simula(h)
    resultats.append(r)
    print(f"{r['frac_sub']*100:>5.1f}% | {h:>6.3f} | {r['W_grav']:>9.6f} | "
          f"{r['W_flot']:>9.6f} | {r['W_fric']:>8.6f} | "
          f"{r['W_el_in']:>8.5f} | {r['W_el_out']:>8.5f} | "
          f"{r['net_el']:>9.5f} | {r['residu']:>6.3f}%")

print("=" * 82)
print()

# Interpretació W_flot
print("=== Interpretació W_flot (flotació contínua) ===")
for r in resultats:
    signe = "↑ font d'energia hidrostàtica" if r['W_flot'] > 1e-6 else \
            ("↓ cost hidrostàtic" if r['W_flot'] < -1e-6 else "≈ 0 conservatiu")
    print(f"  {r['frac_sub']*100:.0f}% sub:  W_flot = {r['W_flot']:+.8f} J/volta  {signe}")

print()
nets = [r['net_el'] for r in resultats]
idx_max = np.argmax(nets)
r_best  = resultats[idx_max]
print(f"Millor balanç: {r_best['frac_sub']*100:.0f}% submergit → net = {r_best['net_el']:+.6f} J/volta")
print()
if r_best['net_el'] > 0.001:
    print(f"★ BALANÇ ELÈCTRIC POSITIU")
    print(f"  Font: W_flot = {r_best['W_flot']:+.6f} J (hidrostàtica)")
    print(f"  Però: cal verificar que W_flot > pèrdues (fricció + conversió)")
    print(f"  W_fric = {r_best['W_fric']:.6f} J  → W_flot cobreix {100*r_best['W_flot']/r_best['W_fric']:.2f}% de la fricció")
else:
    print(f"✓ Balanç negatiu: la flotació contínua és conservativa o quasi-conservativa")
    print(f"  Les pèrdues (fricció + conversió) superen qualsevol aportació hidrostàtica")

# ════════════════════════════════════════════════════════════
#  GRÀFICS
# ════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14, 16))
fig.suptitle("KILÒMETRE v8 — Flotació contínua (esfera)\n"
             "V_submergit(φ) calculat geomètricament, sense discontinuïtat",
             fontsize=12, fontweight='bold')

# 1. W_flot vs % submergit
ax1 = fig.add_subplot(3, 2, 1)
fracs_pct = [r['frac_sub']*100 for r in resultats]
W_flots   = [r['W_flot']*1e6 for r in resultats]  # en µJ
colors_f  = ['green' if w > 0 else 'red' for w in W_flots]
ax1.bar(fracs_pct, W_flots, width=4, color=colors_f, alpha=0.8)
ax1.axhline(0, color='black', lw=1)
ax1.set_xlabel('% submergit (mitjana per volta)')
ax1.set_ylabel('W_flot [µJ/volta]')
ax1.set_title('Treball de flotació per volta\n(continu vs escaló)')
ax1.grid(True, alpha=0.3)

# 2. Net elèctric vs % submergit
ax2 = fig.add_subplot(3, 2, 2)
nets_mJ = [r['net_el']*1000 for r in resultats]
colors_n = ['green' if n > 0 else 'red' for n in nets_mJ]
ax2.bar(fracs_pct, nets_mJ, width=4, color=colors_n, alpha=0.8)
ax2.axhline(0, color='black', lw=1)
ax2.set_xlabel('% submergit')
ax2.set_ylabel('Net elèctric [mJ/volta]')
ax2.set_title('Balanç net elèctric per volta')
ax2.grid(True, alpha=0.3)

# 3-6. Detall per 4 casos
casos_idx = [0, 2, 4, 6]
for ip, ir in enumerate(casos_idx):
    ax = fig.add_subplot(3, 2, 3+ip)
    r  = resultats[ir]
    ph_deg = np.degrees(r['phi'])

    # Volum submergit (normalitzat)
    ax2b = ax.twinx()
    ax2b.fill_between(ph_deg, r['V_sub']/V_total*100, 0,
                      color='royalblue', alpha=0.12, label='% vol. sub.')
    ax2b.set_ylabel('% vol. submergit', color='royalblue', fontsize=8)
    ax2b.tick_params(axis='y', labelcolor='royalblue', labelsize=7)
    ax2b.set_ylim(0, 120)

    ax.plot(ph_deg, r['tau_g'],   'firebrick',   lw=1.2, label='τ_grav')
    ax.plot(ph_deg, r['tau_f'],   'royalblue',   lw=1.5, label='τ_flot', alpha=0.9)
    ax.plot(ph_deg, r['tau_gen'], 'forestgreen', lw=1.2, label='τ_gen')
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlim(0, 360)
    ax.set_title(f"{r['frac_sub']*100:.0f}% sub  |  "
                 f"W_flot={r['W_flot']*1e6:+.3f}µJ  |  "
                 f"net={r['net_el']*1000:+.2f}mJ/v")
    ax.set_xlabel('Angle [°]')
    ax.set_ylabel('Parell [N·m]')
    if ip == 0:
        ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/kilometre_v8_flotacio_continua.png',
            dpi=150, bbox_inches='tight')
print("\nFitxers guardats: kilometre_v8_flotacio_continua.py  |  .png")
