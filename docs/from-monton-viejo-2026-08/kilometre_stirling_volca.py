"""
KILÒMETRE STIRLING — Objecte amb 3 dipòsits + Volcà submarí

L'objecte té 3 cambres interconnectades:
  1. Cambra GAS (He o aire): s'expandeix en calor, es contrau en fred
  2. Cambra LÍQUID (oli/mercuri): fluid de transferència, incompressible
  3. Cambra BUIT PARCIAL: receptor del líquid quan el gas es contrau

Cicle:
  ZONA CALENTA (volcà):
    Gas s'expandeix → empeny líquid cap a cambra buit → V_total↑ → Fn↑
  ZONA FREDA (mar profund):
    Gas es contrau → líquid torna → V_total↓ → Fn↓

Geometria kilòmetre:
  - R gran (10-50m) → temps de trànsit llarg → intercanvi tèrmic complet
  - Tub mig submergit en zona volcànica calenta
  - Tub mig en aigües profundes fredes

Font d'energia: ΔT volcà/mar = 50-300°C
Eficiència màxima teòrica (Carnot): η = 1 - T_fred/T_cal
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

g     = 9.81
R_gas = 8.314   # J/(mol·K)

print("=" * 72)
print("KILÒMETRE STIRLING — Motor tèrmic amb objecte 3 dipòsits")
print("Font d'energia: gradient tèrmic volcà submarí")
print("=" * 72)

# ── Termodinàmica de l'objecte 3 dipòsits ────────────────────
print("\n=== TERMODINÀMICA: objecte 3 cambres ===")

# Gas ideal dins la cambra 1
# P·V = n·R·T → V ∝ T/P (a pressió quasi-constant en el cicle)
# Pressió ≈ constant (profunditat del tub ≈ constant per R gran)

T_fred_K = 278   # 5°C  (aigües profundes)
T_cal_K  = 373   # 100°C (zona volcànica suau — molt conservador)
T_cal_K2 = 573   # 300°C (zona volcànica activa)

# Volum gas a temperatura freda (referència)
V_gas_fred = 1.0  # m³ (normalitzat)

# Volum gas a temperatura calenta (gas ideal)
V_gas_cal  = V_gas_fred * T_cal_K  / T_fred_K
V_gas_cal2 = V_gas_fred * T_cal_K2 / T_fred_K

print(f"\nGas ideal (He o aire comprimit):")
print(f"  T_fred = {T_fred_K-273}°C ({T_fred_K}K)")
print(f"  T_cal  = {T_cal_K-273}°C ({T_cal_K}K)  → V_gas × {V_gas_cal:.2f}  (expansió {(V_gas_cal-1)*100:.0f}%)")
print(f"  T_cal  = {T_cal_K2-273}°C ({T_cal_K2}K) → V_gas × {V_gas_cal2:.2f}  (expansió {(V_gas_cal2-1)*100:.0f}%)")

# Eficiència de Carnot
eta_carnot  = 1 - T_fred_K/T_cal_K
eta_carnot2 = 1 - T_fred_K/T_cal_K2
print(f"\nEficiència Carnot màxima:")
print(f"  ΔT={T_cal_K-T_fred_K}°C: η_Carnot = {eta_carnot*100:.1f}%")
print(f"  ΔT={T_cal_K2-T_fred_K}°C: η_Carnot = {eta_carnot2*100:.1f}%")
print(f"  Eficiència real estimada: ~30-50% de Carnot")
print(f"  → {eta_carnot*0.4*100:.1f}% a {eta_carnot2*0.4*100:.1f}%  (realista)")

# ── Física del kilòmetre Stirling ────────────────────────────
print("\n=== FÍSICA: Fn variable per temperatura ===")

# Fluid exterior: aigüa marina (ρ=1025) + zona volcànica calenta (ρ<1025)
rho_mar_fred = 1025    # kg/m³ aigüa freda profunda
rho_mar_cal  = 960     # kg/m³ aigüa calenta volcànica (~80°C)

# Objecte: cambra de gas + dipòsit líquid
# V_obj = V_gas(T) + V_liquid + V_buit
# En zona freda: V_gas petit → V_obj mínim → Fn mínima (o negativa!)
# En zona calenta: V_gas gran → V_obj màxim → Fn màxima

# Massa objecte (carcassa + líquid fix)
m_carcassa = 500   # kg
V_liquid   = 0.3   # m³  (líquid de transferència, fix)
V_gas_min  = 0.5   # m³  (gas a T_fred)
rho_liquid = 800   # kg/m³ (oli)

m_liquid   = rho_liquid * V_liquid
m_gas      = 0.164 * V_gas_min * 0.1  # He comprimit (negligible)
m_obj_total = m_carcassa + m_liquid + m_gas

def V_obj_total(T_K, V_gas_ref=V_gas_min, T_ref=T_fred_K):
    """Volum total objecte en funció de temperatura."""
    V_gas_T = V_gas_ref * T_K / T_ref   # gas ideal
    return V_gas_T + V_liquid            # + cambra buit absorbeix diferència

def Fn_objecte(T_K, rho_fluid, V_gas_ref=V_gas_min):
    """Força neta (flotació - pes) de l'objecte."""
    V_tot = V_obj_total(T_K, V_gas_ref)
    F_flot = rho_fluid * V_tot * g
    F_pes  = m_obj_total * g
    return F_flot - F_pes

# Cicle complet
print(f"\nObjecte 3 cambres:")
print(f"  Massa total: {m_obj_total:.0f} kg")
print(f"  V_gas_fred:  {V_gas_min:.2f} m³")
print(f"  V_liquid:    {V_liquid:.2f} m³")

print(f"\n{'Zona':>10} | {'T(°C)':>6} | {'ρ_fluid':>8} | {'V_obj(m³)':>10} | "
      f"{'F_flot(kN)':>11} | {'F_pes(kN)':>10} | {'Fn(kN)':>8}")
print("-" * 80)

casos_T = [
    ('Fred profund', T_fred_K, rho_mar_fred),
    ('Temperat',     323,      1000),
    ('Calent 100°',  T_cal_K,  rho_mar_cal),
    ('Volcà 300°',   T_cal_K2, 920),
]

Fns = {}
for nom, T_K, rho_fl in casos_T:
    V_tot  = V_obj_total(T_K)
    F_fl   = rho_fl * V_tot * g / 1000
    F_pes  = m_obj_total * g / 1000
    Fn_kN  = F_fl - F_pes
    Fns[nom] = Fn_kN
    print(f"{nom:>10} | {T_K-273:>6.0f} | {rho_fl:>8.0f} | {V_tot:>10.3f} | "
          f"{F_fl:>11.2f} | {F_pes:>10.2f} | {Fn_kN:>8.2f}")

delta_Fn = Fns['Calent 100°'] - Fns['Fred profund']
delta_Fn2 = Fns['Volcà 300°'] - Fns['Fred profund']
print(f"\nΔFn (100°C cicle): {delta_Fn:.2f} kN")
print(f"ΔFn (300°C cicle): {delta_Fn2:.2f} kN")

# ── Kilòmetre Stirling: dimensions i potència ────────────────
print("\n=== KILÒMETRE STIRLING: POTÈNCIA vs DIMENSIONS ===")

# Per R gran: temps de trànsit llarg → intercanvi tèrmic complet
# Temps trànsit: t = L / v_axial = L / (pas·ω)
# Per R=50m, pas=0.5m/rad, ω=0.7 rad/s: t ≈ 200/(0.5×0.7) ≈ 570s ≈ 10min

pas_sinfin = 0.5   # m/rad

print(f"\nTemps de trànsit objecte pel tub (per intercanvi tèrmic complet):")
print(f"{'R(m)':>6} | {'L=4R(m)':>8} | {'ω_max(r/s)':>11} | "
      f"{'v_ax(m/s)':>10} | {'t_transit(min)':>15} | {'intercanvi?':>12}")
print("-" * 72)

for R in [5, 10, 20, 50, 100]:
    L       = R * 4
    om_max  = np.sqrt(250e6/(7800*R**2))
    v_ax    = pas_sinfin * om_max
    t_trans = (L/2) / v_ax / 60  # minuts (mig recorregut = mig gir)
    ok      = "✓ Complet" if t_trans > 2 else "⚠ Parcial"
    print(f"{R:>6} | {L:>8.0f} | {om_max:>11.3f} | "
          f"{v_ax:>10.3f} | {t_trans:>15.1f} | {ok:>12}")

# ── Potència neta del motor Stirling ─────────────────────────
print(f"\n=== POTÈNCIA NETA: W_net per volta ===")
print(f"W_net = ΔFn × R × ∫_zona_cal cos(φ)dφ  ← però ara ΔFn≠0!")
print(f"""
En el kilòmetre Stirling:
  - Mig gir en zona CALENTA (φ: 0→π): Fn = Fn_cal (gran)
  - Mig gir en zona FREDA   (φ: π→2π): Fn = Fn_fred (petit)

W_net_volta = ∫_0^π Fn_cal·R·cos(φ)dφ + ∫_π^2π Fn_fred·R·cos(φ)dφ
            = Fn_cal·R·[sin(π)-sin(0)] + Fn_fred·R·[sin(2π)-sin(π)]
            = Fn_cal·R·0 + Fn_fred·R·0
            = 0  ← SEGUEIX SENT ZERO!

ESPERA — el problema és el mateix que abans!
La integral de cos(φ) sobre qualsevol interval de 2π és zero.

PERÒ: en el kilòmetre, l'objecte NO segueix una trajectòria circular!
Segueix el tub (recte) mentre el tub gira.
L'altura de l'objecte és: y = R·sin(φ)  (posició radial × sin gir)

Si Fn varia amb la posició angular (calent dalt, fred baix):
  φ ∈ [0,π]   objecte a la meitat SUPERIOR (y>0) → zona calenta
  φ ∈ [π,2π]  objecte a la meitat INFERIOR (y<0) → zona freda

W = ∫_0^π Fn_cal·R·cos(φ)dφ + ∫_π^2π Fn_fred·R·cos(φ)dφ
  = (Fn_cal - Fn_fred)·R·∫_0^π cos(φ)dφ + Fn_fred·R·∮cos(φ)dφ
  = (Fn_cal - Fn_fred)·R·[sin(π)-sin(0)] + 0
  = (Fn_cal - Fn_fred)·R·0 = 0

PERÒ si la zona calenta és la meitat INFERIOR (fons volcànic):
  φ ∈ [0,π]   objecte a DALT → fred → Fn_fred
  φ ∈ [π,2π]  objecte a BAIX → calent → Fn_cal (gran, volcà!)

W = ∫_0^π Fn_fred·R·cos(φ)dφ + ∫_π^2π Fn_cal·R·cos(φ)dφ
  = Fn_fred·R·[0] + Fn_cal·R·[0] = 0

Segueix zero... La geometria circular mata l'asimetria.
""")

print(">>> LA SOLUCIÓ: el tub NO és horitzontal sinó VERTICAL! <<<")
print("""
Si el tub és VERTICAL (eix de rotació horitzontal):
  L'objecte puja en una meitat i baixa en l'altra.
  Zona calenta = part inferior (volcà al fons)
  Zona freda   = part superior (mar fred)

  En pujada (zona freda): Fn_fred → treball NEGATIU (costa pujar)
  En baixada (zona calenta): Fn_cal → treball POSITIU (guanya en baixar)

  ΔW = (Fn_cal - Fn_fred) × L_vertical
     = ΔFn × L

  ARA SÍ! W_net ≠ 0 perquè l'objecte puja per un costat
  i baixa per l'altre, i Fn és diferent a cada costat!
""")

# Càlcul potència motor Stirling vertical
print("=== MOTOR STIRLING VERTICAL: càlcul potència ===\n")
print(f"{'R(m)':>5} {'L(m)':>6} {'ΔFn(kN)':>9} {'W_cicle(kJ)':>12} "
      f"{'rpm':>7} {'P_neta(kW)':>11} {'P_Carnot(kW)':>13}")
print("-" * 68)

# Flux de calor disponible (volcà)
# Q_volca ≈ 100 MW per km² (volcà submarí actiu) → 100 W/m²
Q_flux = 10000   # W/m² (zona volcànica activa prop del cràter)

for R in [5, 10, 20, 50]:
    L         = R * 4
    om_max    = np.sqrt(250e6/(7800*R**2))
    rpm       = om_max * 60 / (2*np.pi)
    A_tub     = np.pi * R**2   # àrea secció

    # W per cicle (objecte puja fred, baixa calent)
    W_cicle   = abs(delta_Fn) * 1000 * L   # J (ΔFn en N × L en m)

    # Potència mecànica neta
    freq      = om_max / (2*np.pi)   # Hz (voltes/s)
    P_mec     = W_cicle * freq       # W

    # Potència limitada pel flux de calor del volcà
    A_intercanvi = 2 * np.pi * R * L   # m² superfície tub
    Q_disponible = Q_flux * A_intercanvi  # W
    P_util       = min(P_mec, Q_disponible * eta_carnot * 0.4)  # W

    print(f"{R:>5} {L:>6.0f} {abs(delta_Fn):>9.2f} {W_cicle/1000:>12.2f} "
          f"{rpm:>7.1f} {P_util/1000:>11.2f} "
          f"{Q_disponible*eta_carnot*0.4/1000:>13.2f}")

print(f"""
Notes:
  ΔFn = {abs(delta_Fn):.2f} kN  (objecte 3 cambres, ΔT=95°C)
  η_Carnot × 0.4 = {eta_carnot*0.4*100:.1f}%  (eficiència realista)
  Q_flux volcà actiu = {Q_flux} W/m²

Limitació real: el flux de calor del volcà és el factor limitant,
no les dimensions mecàniques del tub.
""")

# ── Volcà submarí: dades reals ───────────────────────────────
print("=== VOLCANS SUBMARINES: dades reals ===")
print("""
Tipus de volcà/font hidrotermal:
                                    T_fluid   Flux calor    Profunditat
  Font negra ("black smoker")       350°C     10-100 MW     2000-3000m
  Font blanca ("white smoker")      40-75°C   1-10 MW       500-2000m
  Volcà submarí actiu (flanc)       50-200°C  1-50 W/m²     200-1000m
  Font hidrotermal difusa           5-30°C    0.1-1 W/m²    qualsevol

Per al kilòmetre Stirling:
  → Fonts blanques (40-75°C): ΔT = 35-70°C  → η_Carnot = 11-20%
  → Flanc volcànic (50-200°C): ΔT = 45-195°C → η_Carnot = 14-40%
  → Font negra (350°C): ΔT = 345°C → η_Carnot = 53%  (si sobreviu el tub!)

Localitzacions reals viables:
  ✓ Açores (Mid-Atlantic Ridge): fonts 40-80°C a 200-800m
  ✓ Canàries (flanc volcànic): 30-60°C a 50-500m
  ✓ Mediterrani (Santorini, Panarèa): 50-100°C a 5-200m
  ✓ Japó/Pacífic: fonts actives a profunditats diverses
  ✓ Islàndia submarina: geotèrmica + oceà
""")

# ── Gràfics ──────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
fig.suptitle("KILÒMETRE STIRLING — Motor tèrmic oceànic\n"
             "Objecte 3 dipòsits + Volcà submarí com a font de calor",
             fontsize=13, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.38)

# 1. Volum objecte vs temperatura
ax = fig.add_subplot(gs[0, 0])
Ts = np.linspace(270, 600, 200)
Vs = [V_obj_total(T) for T in Ts]
ax.plot(Ts-273, Vs, 'firebrick', lw=2.5)
ax.axvline(T_fred_K-273, color='steelblue', lw=2, ls='--',
           label=f'T_fred={T_fred_K-273}°C')
ax.axvline(T_cal_K-273,  color='orange',    lw=2, ls='--',
           label=f'T_cal={T_cal_K-273}°C')
ax.axvline(T_cal_K2-273, color='red',       lw=2, ls='--',
           label=f'T_volcà={T_cal_K2-273}°C')
ax.fill_between(Ts-273, V_obj_total(T_fred_K),
                [V_obj_total(T) for T in Ts],
                where=[T>=T_fred_K for T in Ts],
                alpha=0.2, color='orange', label='ΔV útil')
ax.set_xlabel('Temperatura [°C]')
ax.set_ylabel('V_objecte [m³]')
ax.set_title('Volum objecte vs temperatura\n(cambra gas + liquid + buit)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# 2. Fn vs temperatura
ax = fig.add_subplot(gs[0, 1])
Fns_T_fred = [Fn_objecte(T, rho_mar_fred)/1000 for T in Ts]
Fns_T_cal  = [Fn_objecte(T, rho_mar_cal)/1000  for T in Ts]
ax.plot(Ts-273, Fns_T_fred, 'steelblue', lw=2, label=f'Fluid fred {rho_mar_fred}kg/m³')
ax.plot(Ts-273, Fns_T_cal,  'firebrick', lw=2, label=f'Fluid cal {rho_mar_cal}kg/m³')
ax.axhline(0, color='black', lw=1)
ax.axvline(T_fred_K-273, color='steelblue', lw=1, ls=':')
ax.axvline(T_cal_K-273,  color='orange',    lw=1, ls=':')
ax.fill_between(Ts-273,
                Fns_T_fred, Fns_T_cal,
                alpha=0.2, color='purple', label='ΔFn disponible')
ax.set_xlabel('T objecte [°C]')
ax.set_ylabel('Fn [kN]')
ax.set_title('Força neta vs temperatura\n(positiu = flota, negatiu = enfonsa)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# 3. Esquema motor (tub vertical)
ax = fig.add_subplot(gs[0, 2])
ax.set_xlim(0, 10); ax.set_ylim(0, 12); ax.axis('off')
ax.set_title('Esquema Motor Stirling Vertical', fontsize=10)

# Fons del mar (volcà)
ax.add_patch(plt.Rectangle((0,0), 10, 1.5,
             color='#ff6600', alpha=0.7))
ax.text(5, 0.75, '🌋 VOLCÀ SUBMARÍ\n(T=100-300°C)',
        ha='center', va='center', fontsize=9, fontweight='bold')

# Aigüa freda (dalt)
ax.add_patch(plt.Rectangle((0,9.5), 10, 2.5,
             color='#aaddff', alpha=0.6))
ax.text(5, 10.5, '❄️ AIGÜA FREDA PROFUNDA\n(T=5-10°C)',
        ha='center', va='center', fontsize=9)

# Tub vertical
ax.add_patch(plt.Rectangle((3.8,1.5), 0.2, 8,
             color='gray', alpha=0.8))  # paret esquerra
ax.add_patch(plt.Rectangle((6.0,1.5), 0.2, 8,
             color='gray', alpha=0.8))  # paret dreta

# Zona calenta (dins tub, baix)
ax.add_patch(plt.Rectangle((4.0,1.5), 2, 3.5,
             color='#ffaa44', alpha=0.4))
ax.text(5, 3.2, 'ZONA\nCALENTA\nFn↑ V↑',
        ha='center', va='center', fontsize=8, color='darkred')

# Zona freda (dins tub, dalt)
ax.add_patch(plt.Rectangle((4.0,5.0), 2, 4.5,
             color='#aaeeff', alpha=0.4))
ax.text(5, 7.2, 'ZONA\nFREDA\nFn↓ V↓',
        ha='center', va='center', fontsize=8, color='darkblue')

# Objecte en zona calenta (gran, flota molt)
ax.add_patch(mpatches.Ellipse((5, 2.8), 1.6, 0.7,
             color='#ffdd00', alpha=0.9, zorder=5))
ax.text(5, 2.8, 'Obj\ncalent\n(gran)',
        ha='center', va='center', fontsize=7, fontweight='bold')

# Objecte en zona freda (petit, poc empeny)
ax.add_patch(mpatches.Ellipse((5, 7.5), 1.0, 0.5,
             color='#88ccff', alpha=0.9, zorder=5))
ax.text(5, 7.5, 'Obj\nfred\n(petit)',
        ha='center', va='center', fontsize=7)

# Fletxes cicle
ax.annotate('', (4.2,6.8), (4.2,3.2),
            arrowprops=dict(arrowstyle='->', color='blue',
                           lw=2, connectionstyle='arc3,rad=-0.3'))
ax.annotate('', (5.8,3.5), (5.8,7.0),
            arrowprops=dict(arrowstyle='->', color='red',
                           lw=2, connectionstyle='arc3,rad=-0.3'))
ax.text(2.8, 5, 'Puja\nfred', ha='center', fontsize=7, color='blue')
ax.text(7.2, 5, 'Baixa\ncalent', ha='center', fontsize=7, color='red')

# Sinfín + generador
ax.text(5, 11.2, '⚡ GENERADOR (sinfín)',
        ha='center', fontsize=9, fontweight='bold', color='darkgreen')

# 4. Potència vs R (motor Stirling)
ax = fig.add_subplot(gs[1, 0])
Rs_p    = np.linspace(2, 100, 100)
Ps_p    = []
Qs_p    = []
for R_p in Rs_p:
    L_p      = R_p * 4
    om_p     = np.sqrt(250e6/(7800*R_p**2))
    freq_p   = om_p/(2*np.pi)
    W_c_p    = abs(delta_Fn)*1000*L_p
    P_mec_p  = W_c_p * freq_p
    A_p      = 2*np.pi*R_p*L_p
    Q_p      = Q_flux * A_p * eta_carnot * 0.4
    Ps_p.append(min(P_mec_p, Q_p)/1e6)
    Qs_p.append(Q_p/1e6)
ax.loglog(Rs_p, Ps_p, 'firebrick', lw=2.5, label='P_neta motor')
ax.loglog(Rs_p, Qs_p, 'orange', lw=2, ls='--', label='Límit flux volcà')
ax.axhline(1,   color='gray', lw=1, ls=':', label='1 MW')
ax.axhline(100, color='gray', lw=1, ls='--', label='100 MW')
ax.set_xlabel('R [m]'); ax.set_ylabel('Potència [MW]')
ax.set_title('Potència motor Stirling vs R\n(limitat pel flux del volcà)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which='both')

# 5. Eficiència Carnot vs ΔT
ax = fig.add_subplot(gs[1, 1])
DTs   = np.linspace(5, 350, 200)
T_f   = 278  # K (5°C)
etas  = 1 - T_f/(T_f + DTs)
etas_real = etas * 0.4
ax.plot(DTs, etas*100,      'steelblue', lw=2, label='η Carnot teòric')
ax.plot(DTs, etas_real*100, 'firebrick', lw=2, label='η real (40% Carnot)')
ax.axvline(20,  color='cyan',   lw=1, ls='--', label='OTEC ΔT=20°C')
ax.axvline(70,  color='orange', lw=1, ls='--', label='Font blanca ΔT=70°C')
ax.axvline(200, color='red',    lw=1, ls='--', label='Volcà ΔT=200°C')
ax.axvline(345, color='darkred',lw=1, ls='--', label='Font negra ΔT=345°C')
ax.set_xlabel('ΔT [°C]'); ax.set_ylabel('Eficiència [%]')
ax.set_title('Eficiència Carnot vs ΔT\n(T_fred=5°C referència)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

# 6. Mapa de volcans viables
ax = fig.add_subplot(gs[1, 2])
ax.axis('off')
ax.set_title('Fonts de calor oceàniques\nper al Kilòmetre Stirling',
             fontsize=10, fontweight='bold')

fonts = [
    ('OTEC (gradient tèrmic)',   '20°C',  '0.1 W/m²', '3-8%',  '✓ Global'),
    ('Font hidrotermal difusa',  '30°C',  '1 W/m²',   '5-10%', '✓ Disponible'),
    ('Flanc volcànic',           '70°C',  '10 W/m²',  '12-20%','✓ Açores,Canàries'),
    ('Font blanca (smoker)',     '75°C',  '100 W/m²', '15-22%','✓ Cuidat materials'),
    ('Volcà actiu (flanc)',      '200°C', '1000 W/m²','28-40%','⚠ Difícil accés'),
    ('Font negra (smoker)',      '350°C', '10000 W/m²','53%',  '✗ Massa agressiu'),
]

y = 0.92
ax.text(0.05, y, f"{'Font':^28} {'ΔT':^7} {'Flux':^11} {'η_real':^8} {'Viabilitat':^18}",
        transform=ax.transAxes, fontsize=7.5, fontweight='bold',
        fontfamily='monospace')
y -= 0.06
for f_nom, dT, flux, eta_r, viab in fonts:
    color = 'darkgreen' if '✓' in viab else ('orange' if '⚠' in viab else 'red')
    ax.text(0.05, y,
            f"{f_nom:<28} {dT:^7} {flux:^11} {eta_r:^8} {viab:<18}",
            transform=ax.transAxes, fontsize=7.5,
            fontfamily='monospace', color=color)
    y -= 0.10

ax.text(0.05, y-0.04,
        "★ RECOMANACIÓ: Flanc volcànic (Açores/Canàries)\n"
        "  ΔT=70°C, η≈15%, flux 10W/m², profunditat viable\n"
        "  Materials estàndard (acer inox), accés submergible",
        transform=ax.transAxes, fontsize=8,
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.savefig('/mnt/user-data/outputs/kilometre_stirling_volca.png',
            dpi=150, bbox_inches='tight')
print("Fitxers guardats: kilometre_stirling_volca.py  |  .png")
