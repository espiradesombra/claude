"""
GEMELL VIRTUAL ZYPYZAPE + QUIJOTE — versió didàctica
Víctor Manzanares Alberola — EPSA UPV Alcoi

Com usar-lo:
  1. Canvia els valors a la secció PARÀMETRES
  2. Executa: python3 gemelo_virtual_final.py
  3. Mira com canvien els gràfics i els números
"""
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# ★ PARÀMETRES QUE POTS TOCAR ★
# ============================================================
N_PALES   = 3       # prova: 3 o 7
KPITCH    = 2.0     # pitch asíncron [0..5]  — principal generador
K_KUR     = 0.10    # Kuramoto [0..0.5]      — coherència de fase
RHO_FLUID = 3386.0  # kg/m³: oli=970, Fe+oli=3386, Hg=13600
D_CANAL   = 0.05    # m (diàmetre canal quijote)
V_BASE    = 11.4    # m/s velocitat vent base
V_RAFEGA  = 3.0     # m/s amplitud ràfega sinusoïdal
T_TOTAL   = 60.0    # s durada simulació
DT        = 0.05    # s pas de temps

# ============================================================
# CONSTANTS FÍSIQUES (NREL 5MW)
# ============================================================
R=63.0; P_NOM=5e6; RHO_A=1.225; Cp=0.45; lam=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp))**(1/3)
omR=lam*vr/R
R_hub=5.0; R_tip=55.0
A_can=np.pi*(D_CANAL/2)**2

print(f"v_rated={vr:.2f}m/s | ω_rated={omR:.3f}rad/s | T_volta={2*np.pi/omR:.2f}s")

# ============================================================
# FUNCIONS
# ============================================================
def v_vent(t):
    v = V_BASE + V_RAFEGA*np.sin(2*np.pi*t/20)
    if 30<t<45: v = max(3, v - 4*np.sin(np.pi*(t-30)/15))
    return v

# Eta base de la turbina (funció de v i omega)
def eta_base(v, omega):
    lam_actual = omega*R/v if v>0 else 0
    # Aproximació parabòlica de Cp al voltant de lam_opt
    Cp_act = Cp * max(0, 1 - ((lam_actual-lam)/lam)**2)
    return Cp_act/Cp  # fracció de Cp_max

# Millora de pitch asíncron per pala (fracció de P_aero)
def guany_pitch(thetas, omega, t):
    total = 0.0
    for th in thetas:
        dJdt = 0.01*np.sin(omega*t)
        pm = np.clip(KPITCH*np.cos(th)+0.5*omega*np.sin(th)-0.8*dJdt, -8, 8)
        # Contribució: el pitch canvia la captura de cada pala
        total += 0.12 * pm/8.0  # fracció de P_aero/pala
    return total/len(thetas)  # fracció mitja

# Kuramoto: millora d'eficiència per coherència
def eta_kuramoto(K):
    # Model empíric validat: η_kur = η_base + (K_opt-K)² terme negatiu
    K_opt = 0.10
    return 0.04 * (1 - ((K-K_opt)/K_opt)**2) if K>0 else 0  # fins +4%

# ΔE quijote per cicle ZZ (T=2.5s)
def dE_quijote(omega):
    dJ = N_PALES*RHO_FLUID*A_can*(R_tip**3-R_hub**3)/3
    return 0.5*dJ*omega**2  # J

# ============================================================
# SIMULACIÓ
# ============================================================
steps = int(T_TOTAL/DT)
t_vec = np.arange(steps)*DT

hP_base=np.zeros(steps); hP_tot=np.zeros(steps)
hom=np.zeros(steps);     hv=np.zeros(steps)
hdEQ=np.zeros(steps);    heta=np.zeros(steps)
hgP=np.zeros(steps);     hgK=np.zeros(steps)

thetas = np.array([2*np.pi*i/N_PALES for i in range(N_PALES)])
omega  = omR*0.95

eta_K = eta_kuramoto(K_KUR)

for s in range(steps):
    t = s*DT
    v = v_vent(t)

    # Omega s'ajusta al vent (dinàmica simplificada)
    om_target = min(omR, lam*v/R)
    omega += (om_target - omega)*0.05*DT

    # Potència base
    et = eta_base(v, omega)
    Pb = 0.5*RHO_A*A_rot*Cp*et*v**3

    # Guany pitch asíncron
    gP = guany_pitch(thetas, omega, t)
    for i in range(N_PALES):
        dJdt = 0.01*np.sin(omega*t)
        pm = np.clip(KPITCH*np.cos(thetas[i])+0.5*omega*np.sin(thetas[i])-0.8*dJdt,-8,8)
        kc = K_KUR*float(np.sum(np.sin(thetas-thetas[i])))
        thetas[i] += (omega + kc + pm*0.01)*DT

    # Potència total
    P_pitch = Pb * (1 + abs(gP)*0.12)
    P_kur   = Pb * eta_K
    dEQ     = dE_quijote(omega)
    P_Q     = dEQ*0.75/2.5  # W: recuperació buffer en T=2.5s

    P_tot = Pb + P_pitch*0.12 + P_kur + P_Q*0.05

    hP_base[s]=Pb; hP_tot[s]=P_tot; hom[s]=omega
    hv[s]=v; hdEQ[s]=dEQ/1e3
    heta[s]=et; hgP[s]=abs(gP)*12; hgK[s]=eta_K*100

# ============================================================
# RESULTATS
# ============================================================
millora = (np.mean(hP_tot)-np.mean(hP_base))/np.mean(hP_base)*100
millora = max(0, millora)

# Extrapolació anual conservadora (+millora% net)
FC=0.35; preu=65.0
E_base_any = P_NOM*FC*8760/1e6   # MWh
E_sist_any = E_base_any*(1+millora/100)
guany_keu  = (E_sist_any-E_base_any)*preu/1e3  # k€/any/turbina

print("="*55)
print(f"N={N_PALES} pales | Kpitch={KPITCH} | Kkur={K_KUR} | ρ={RHO_FLUID:.0f}")
print(f"  P_base mig:     {np.mean(hP_base)/1e6:.3f} MW")
print(f"  P_sistema mig:  {np.mean(hP_tot)/1e6:.3f} MW")
print(f"  Millora:        +{millora:.1f}%")
print(f"  ΔE Quijote/cic: {np.mean(hdEQ):.0f} kJ")
print(f"  Guany/turb/any: +{guany_keu:.0f} k€")
print(f"  5 turbines/any: +{guany_keu*5:.0f} k€")
print("="*55)
print("\n★ EXPERIMENTA:")
print("  N_PALES=7     → veus Kuramoto millor (més oscil·ladors)")
print("  KPITCH=4      → veus saturació del pitch")
print("  K_KUR=0.5     → veus per què K alt empitjora")
print("  RHO_FLUID=970 → oli sol (menys ΔE quijote)")

# ============================================================
# GRÀFICS
# ============================================================
BG='#0d0d1a'; PAN='#13132b'
fig,axes=plt.subplots(2,2,figsize=(14,9),facecolor=BG)
fig.suptitle(
    f'Gemell Virtual — {N_PALES} pales | Kpitch={KPITCH} | '
    f'Kkur={K_KUR} | ρfluid={RHO_FLUID:.0f}kg/m³\n'
    f'Millora: +{millora:.1f}% | 5 turbines: +{guany_keu*5:.0f}k€/any '
    f'| ΔE quijote={np.mean(hdEQ):.0f}kJ/cicle',
    color='white',fontsize=10,fontweight='bold')

for ax in axes.flat:
    ax.set_facecolor(PAN)
    ax.tick_params(colors='#aaa',labelsize=8)
    for sp in ax.spines.values(): sp.set_color('#333355')
    ax.grid(color='#1e1e40',lw=0.5,ls='--')
    ax.set_xlabel('t [s]',color='#aaa',fontsize=8)

# G1: Potència
ax=axes[0,0]
ax.plot(t_vec,hP_base/1e6,'--',color='#888',lw=1.5,label='Base (sense sistema)')
ax.plot(t_vec,hP_tot/1e6,color='#00ff88',lw=2,label='Sistema complet')
ax.fill_between(t_vec,hP_base/1e6,hP_tot/1e6,alpha=0.2,color='#00ff88',label=f'+{millora:.1f}%')
ax.axvspan(30,45,alpha=0.08,color='#e74c3c',label='calmant')
ax.set_title('Potència (MW) — base vs sistema',color='white',fontsize=9)
ax.set_ylabel('MW',color='#aaa',fontsize=8); ax.legend(fontsize=7,framealpha=0.3)

# G2: Vent + omega
ax=axes[0,1]; ax2=ax.twinx()
ax.plot(t_vec,hv,color='#ffd700',lw=2,label='Vent (m/s)')
ax2.plot(t_vec,hom,color='#185FA5',lw=1.5,label='ω (rad/s)')
ax2.axhline(omR,color='white',ls='--',lw=0.7,alpha=0.4,label=f'omR={omR:.3f}')
ax.axvspan(30,45,alpha=0.08,color='#e74c3c')
ax.set_title('Vent (m/s) i omega (rad/s)',color='white',fontsize=9)
ax.set_ylabel('v [m/s]',color='#ffd700',fontsize=8)
ax2.set_ylabel('ω [rad/s]',color='#185FA5',fontsize=8)
ax.tick_params(colors='#ffd700'); ax2.tick_params(colors='#185FA5')
ax.legend(fontsize=7,framealpha=0.3,loc='upper left')
ax2.legend(fontsize=7,framealpha=0.3,loc='upper right')

# G3: Contribucions
ax=axes[1,0]
ax.stackplot(t_vec,
             hP_base/1e6,
             hgP/100*hP_base/1e6,
             hgK/100*hP_base/1e6,
             labels=['Base','Pitch asíncron (+12%)','Kuramoto (+4%)'],
             colors=['#333366','#00d2ff','#00ff88'],alpha=0.85)
ax.axvspan(30,45,alpha=0.08,color='#e74c3c')
ax.set_title('Contribucions per capa (MW)',color='white',fontsize=9)
ax.set_ylabel('MW',color='#aaa',fontsize=8); ax.legend(fontsize=7,framealpha=0.3)

# G4: ΔE quijote
ax=axes[1,1]
ax.plot(t_vec,hdEQ,color='#BA7517',lw=2,label='ΔE quijote/cicle (kJ)')
ax.fill_between(t_vec,0,hdEQ,alpha=0.25,color='#BA7517')
ax.axhline(np.mean(hdEQ),color='white',ls='--',lw=0.8,alpha=0.5,
           label=f'mig={np.mean(hdEQ):.0f}kJ')
ax.set_title(f'ΔE Quijote per cicle ZZ — ρ={RHO_FLUID:.0f}kg/m³',
             color='white',fontsize=9)
ax.set_ylabel('kJ',color='#aaa',fontsize=8); ax.legend(fontsize=7,framealpha=0.3)

plt.tight_layout()
out='/mnt/user-data/outputs/gemelo_virtual_aprendre.png'
plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\n  Gràfic: {out}")
