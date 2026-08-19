"""
KILÒMETRE v7 — Bateria submarina amb flotació asimètrica
Eix del tub desplaçat respecte la superfície de l'aigua.
Percentatge submergit variable (40% a 60%).

Física nova:
  - Gravetat:  τ_g = -m·g·R·cos(φ)          → W_grav = 0 sempre
  - Flotació:  τ_f = +F_flot·R·cos(φ)        → W_flot ≠ 0 si asimètric
    actua quan: R·sin(φ) + h_eix < 0
    (objecte per sota la superfície)
  - La font d'energia real és la DIFERÈNCIA DE PRESSIÓ
    entre l'aire i l'aigua (energia potencial hidràulica)

h_eix = alçada de l'eix respecte la superfície:
  h_eix > 0 → eix per sobre → menys del 50% submergit
  h_eix < 0 → eix per sota  → més del 50% submergit
  h_eix = 0 → exactament 50% (W_flot = 0)
"""

import numpy as np
import matplotlib.pyplot as plt

# ════════════════════════════════════════════════════════════
#  PARÀMETRES
# ════════════════════════════════════════════════════════════
m        = 10.0    # kg
R        = 0.15    # m
g        = 9.81    # m/s²
rho_agua = 1000.0  # kg/m³
rho_obj  = 2000.0  # kg/m³   (objecte més dens que l'aigua)
V_obj    = m / rho_obj       # m³
F_flot   = rho_agua * V_obj * g   # N  (50% del pes)

m_tub    = 2.0;  R_tub = 0.18
I        = 0.5*m_tub*R_tub**2 + m*R**2

beta     = 0.02    # N·m·s/rad  fricció viscosa
omega_ref = 3.0    # rad/s  (mantingut pel motor)
alfa_gen  = 0.60   # fracció de frenada en fase favorable
eta_gen   = 0.90
eta_mot   = 0.92

print(f"Pes objecte:  {m*g:.2f} N")
print(f"Flotació màx: {F_flot:.2f} N  ({100*F_flot/(m*g):.0f}% del pes)")
print(f"R = {R} m  |  omega = {omega_ref} rad/s  |  I = {I:.4f} kg·m²")
print()

# ════════════════════════════════════════════════════════════
#  FUNCIÓ: simula una volta amb h_eix donat
# ════════════════════════════════════════════════════════════
def simula_volta(h_eix, alfa_gen=0.60, N=100000):
    """
    Retorna balanç energètic per una volta completa.
    h_eix: alçada eix sobre superfície [m]
           positiu → eix per sobre (menys submergit)
           negatiu → eix per sota  (més submergit)
    """
    phi  = np.linspace(0, 2*np.pi, N, endpoint=False)
    dphi = 2*np.pi / N
    om   = omega_ref

    # Posició vertical objecte (positiu = per sobre l'aigua)
    y_obj = R * np.sin(phi) + h_eix

    # Submergit quan y_obj < 0
    submergit = y_obj < 0
    frac_sub  = np.mean(submergit)

    # Parells
    tau_g  = -m * g * R * np.cos(phi)
    tau_f  = np.where(submergit, +F_flot * R * np.cos(phi), 0.0)
    tau_fr = -beta * om * np.ones(N)

    # Parell net gravitatori + flotació
    tau_gf = tau_g + tau_f

    # Generador: embriagat quan tau_gf > 0 (força neta ajuda el gir)
    # i quan objecte va prou de pressa
    tau_gen = np.where(tau_gf > 0, -alfa_gen * tau_gf, 0.0)

    # Motor: compensa tot per mantindre omega constant
    tau_mot = -(tau_gf + tau_gen + tau_fr)

    # Treball mecànic per element
    dW_g   = tau_g   * om * dphi
    dW_f   = tau_f   * om * dphi
    dW_fr  = tau_fr  * om * dphi
    dW_gen = tau_gen * om * dphi   # ≤ 0 (extreu mecànica)
    dW_mot = tau_mot * om * dphi   # + si aporta, - si recupera

    W_grav = np.sum(dW_g)
    W_flot = np.sum(dW_f)
    W_fric = np.sum(np.abs(dW_fr))

    # Elèctric generador
    W_el_gen = np.sum(np.where(dW_gen < 0, np.abs(dW_gen) * eta_gen, 0.0))

    # Elèctric motor (separat: consum i recuperació)
    W_el_in  = np.sum(np.where(dW_mot > 0, dW_mot / eta_mot, 0.0))
    W_el_rec = np.sum(np.where(dW_mot < 0, np.abs(dW_mot) * eta_gen, 0.0))
    W_el_out = W_el_gen + W_el_rec

    # Verificació 1a Llei (ΔKE=0 perquè omega constant)
    check = W_grav + W_flot + np.sum(dW_mot) + np.sum(dW_gen) - W_fric
    residu = abs(check) / max(W_fric, W_el_in, 0.001) * 100

    return {
        'h_eix':    h_eix,
        'frac_sub': frac_sub,
        'W_grav':   W_grav,
        'W_flot':   W_flot,
        'W_fric':   W_fric,
        'W_el_in':  W_el_in,
        'W_el_out': W_el_out,
        'net_el':   W_el_out - W_el_in,
        'residu%':  residu,
        'tau_g':    tau_g,
        'tau_f':    tau_f,
        'tau_gen':  tau_gen,
        'tau_mot':  tau_mot,
        'phi':      phi,
        'submergit':submergit,
    }

# ════════════════════════════════════════════════════════════
#  ESCOMBRAT DE h_eix (% submergit 30% → 70%)
# ════════════════════════════════════════════════════════════
# h_eix tal que frac_sub = target:
# frac_sub ≈ arccos(h_eix/R)/π  (aproximació geomètrica)
# h_eix = R·cos(π·frac_sub)
fracs_target = np.linspace(0.30, 0.70, 9)
h_eixos      = R * np.cos(np.pi * fracs_target)

print("=" * 78)
print("ESCOMBRAT % SUBMERGIT")
print(f"{'% sub':>6} | {'h_eix':>6} | {'W_grav':>8} | {'W_flot':>8} | "
      f"{'W_fric':>7} | {'el_in':>8} | {'el_out':>8} | {'net_el':>9} | {'res%':>5}")
print("-" * 78)

resultats = []
for frac, h in zip(fracs_target, h_eixos):
    r = simula_volta(h)
    resultats.append(r)
    print(f"{r['frac_sub']*100:>5.1f}% | {h:>6.3f} | {r['W_grav']:>8.5f} | "
          f"{r['W_flot']:>8.5f} | {r['W_fric']:>7.5f} | "
          f"{r['W_el_in']:>8.5f} | {r['W_el_out']:>8.5f} | "
          f"{r['net_el']:>9.5f} | {r['residu%']:>5.2f}%")

print("=" * 78)
print()

# Resultat màxim
nets    = [r['net_el'] for r in resultats]
idx_max = np.argmax(nets)
r_best  = resultats[idx_max]
print(f"Millor cas: {r_best['frac_sub']*100:.1f}% submergit  →  "
      f"net = {r_best['net_el']:+.5f} J/volta")
print(f"  W_flot = {r_best['W_flot']:+.5f} J  (font real d'energia: pressió hidrostàtica)")
print(f"  W_grav = {r_best['W_grav']:+.5f} J  (ha de ser ~0)")
print()
if r_best['net_el'] > 0:
    print(f"★ BALANÇ POSITIU: {r_best['net_el']:+.5f} J/volta")
    print(f"  Potència neta: {r_best['net_el']*omega_ref/(2*np.pi)*1000:+.2f} mW")
    print(f"  Font: diferència de pressió aire/aigua (NO la gravetat)")
else:
    print(f"✓ BALANÇ NEGATIU en tots els casos: la flotació asimètrica")
    print(f"  no compensa les pèrdues de conversió")

# ════════════════════════════════════════════════════════════
#  GRÀFICS
# ════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14, 16))
fig.suptitle("KILÒMETRE v7 — Flotació asimètrica\n"
             "Font d'energia: diferència de pressió hidrostàtica",
             fontsize=12, fontweight='bold')

# 1. Net elèctric vs % submergit
ax1 = fig.add_subplot(3, 2, 1)
fracs_pct = [r['frac_sub']*100 for r in resultats]
nets_mJ   = [r['net_el']*1000 for r in resultats]
colors_bar = ['green' if n > 0 else 'red' for n in nets_mJ]
ax1.bar(fracs_pct, nets_mJ, width=4, color=colors_bar, alpha=0.8)
ax1.axhline(0, color='black', lw=1)
ax1.set_xlabel('% submergit')
ax1.set_ylabel('Net elèctric [mJ/volta]')
ax1.set_title('Balanç net vs % submergit')
ax1.grid(True, alpha=0.3)

# 2. W_flot i W_fric vs % submergit
ax2 = fig.add_subplot(3, 2, 2)
W_flots = [r['W_flot']*1000 for r in resultats]
W_frics = [-r['W_fric']*1000 for r in resultats]
ax2.plot(fracs_pct, W_flots, 'b-o', lw=2, label='W_flot (mJ)')
ax2.plot(fracs_pct, W_frics, 'r-o', lw=2, label='-W_fric (mJ)')
ax2.axhline(0, color='gray', lw=0.5)
ax2.set_xlabel('% submergit')
ax2.set_ylabel('Energia [mJ/volta]')
ax2.set_title('Flotació vs Fricció')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# 3-6. Parells per 4 casos representatius
casos = [0, 2, 4, 6]   # 30%, 40%, 50%, 60%
for idx_plot, idx_r in enumerate(casos):
    ax = fig.add_subplot(3, 2, 3 + idx_plot)
    r  = resultats[idx_r]
    ph_deg = np.degrees(r['phi'])

    ax.plot(ph_deg, r['tau_g'],   'firebrick',   lw=1.2, label='τ_grav', alpha=0.8)
    ax.plot(ph_deg, r['tau_f'],   'royalblue',   lw=1.2, label='τ_flot', alpha=0.8)
    ax.plot(ph_deg, r['tau_gen'], 'forestgreen', lw=1.2, label='τ_gen',  alpha=0.8)
    ax.plot(ph_deg, r['tau_mot'], 'purple',      lw=1,   label='τ_mot',  alpha=0.6, ls='--')

    # Zona submergida
    sub = r['submergit']
    in_sub = False; x0 = 0
    for j in range(len(ph_deg)):
        if sub[j] and not in_sub:
            x0 = ph_deg[j]; in_sub = True
        elif not sub[j] and in_sub:
            ax.axvspan(x0, ph_deg[j], alpha=0.08, color='blue')
            in_sub = False

    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlim(0, 360)
    ax.set_title(f"{r['frac_sub']*100:.0f}% submergit  |  "
                 f"net={r['net_el']*1000:+.2f} mJ/volta")
    ax.set_xlabel('Angle [°]')
    ax.set_ylabel('Parell [N·m]')
    if idx_plot == 0:
        ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/kilometre_v7_flotacio.png', dpi=150, bbox_inches='tight')
print("Fitxers guardats: kilometre_v7_flotacio.py  |  kilometre_v7_flotacio.png")
