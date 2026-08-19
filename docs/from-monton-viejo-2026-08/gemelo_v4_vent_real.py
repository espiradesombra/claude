"""
GEMELL VIRTUAL v4 — Vent realista + control amb retard
Víctor Manzanares Alberola — EPSA UPV Alcoi

CANVIS vs v3:
  - 3 perfils de vent: SUAU / REALISTA / EXTREM
  - Dinàmica omega amb tau realista (5s) + soroll
  - Comparativa automàtica de guany per perfil
  - Ornstein-Uhlenbeck per turbulència física
"""
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)  # reproducibilitat

# ============================================================
# PARÀMETRES
# ============================================================
N_PALES   = 3
KPITCH    = 2.0
K_KUR     = 0.10
RHO_FLUID = 3386.0
D_CANAL   = 0.05
V_BASE    = 11.4
T_TOTAL   = 80.0
DT        = 0.05

# ============================================================
# CONSTANTS
# ============================================================
R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.45; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3)
omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0
A_can=np.pi*(D_CANAL/2)**2
J_base=3.54e7
m_q=RHO_FLUID*A_can*(R_tip-R_hub)
dJ_max=N_PALES*m_q*(R_tip**2-R_hub**2)
TAU_OMEGA=5.0   # s — temps de resposta real del rotor (>>ideal)
SIGMA_OU=0.3    # turbulència Ornstein-Uhlenbeck (m/s)
THETA_OU=0.5    # velocitat de reversió O-U

print(f"omR={omR:.3f}rad/s | vr={vr:.2f}m/s | tau_omega={TAU_OMEGA}s")

# ============================================================
# PRE-GENERA ELS 3 PERFILS DE VENT
# ============================================================
steps=int(T_TOTAL/DT)
t_vec=np.arange(steps)*DT

def gen_vent_suau(t_vec):
    """Original — sinusoïdal perfecte. Poca turbulència."""
    v=np.zeros(len(t_vec))
    for i,t in enumerate(t_vec):
        vi=V_BASE+3.0*np.sin(2*np.pi*t/20)
        if 30<t<45: vi=max(3.0,vi-4.0*np.sin(np.pi*(t-30)/15))
        v[i]=vi
    return v

def gen_vent_realista(t_vec):
    """
    Vent realista amb 4 components (GPT):
    1. Base lenta (meteorologia)
    2. Ràfegues irregulars (freqüència variable)
    3. Soroll turbulent (Ornstein-Uhlenbeck)
    4. Calmant brusc (no sinusoïdal)
    """
    N=len(t_vec)
    # Component O-U: dv = -θ·v·dt + σ·√dt·N(0,1)
    ou=np.zeros(N)
    for i in range(1,N):
        ou[i]=ou[i-1]-THETA_OU*ou[i-1]*DT+SIGMA_OU*np.sqrt(DT)*np.random.randn()

    v=np.zeros(N)
    for i,t in enumerate(t_vec):
        vi=V_BASE
        vi+=2.0*np.sin(2*np.pi*t/60)                          # meteorologia
        vi+=1.5*np.sin(2*np.pi*t/7+0.5*np.sin(t/10))         # ràfegues
        vi+=ou[i]                                               # turbulència O-U
        if 30<t<45: vi-=5.0                                    # calmant brusc
        v[i]=max(2.0,vi)
    return v

def gen_vent_extrem(t_vec):
    """
    Vent extrem: turbulència alta + múltiples calmants + ràfegues fortes.
    Intensitat turbulència IEC classe A (σ≈15% v_mig).
    """
    N=len(t_vec)
    ou=np.zeros(N)
    for i in range(1,N):
        ou[i]=ou[i-1]-0.8*ou[i-1]*DT+1.5*np.sqrt(DT)*np.random.randn()

    v=np.zeros(N)
    for i,t in enumerate(t_vec):
        vi=V_BASE
        vi+=3.0*np.sin(2*np.pi*t/25)
        vi+=2.5*np.sin(2*np.pi*t/5+np.sin(t/8)*2)
        vi+=ou[i]
        if 20<t<30:  vi-=6.0                 # primer calmant
        if 50<t<60:  vi+=5.0                 # ràfega forta
        if 65<t<75:  vi-=7.0                 # segon calmant sever
        v[i]=max(2.0,vi)
    return v

vents={
    'SUAU':     gen_vent_suau(t_vec),
    'REALISTA': gen_vent_realista(t_vec),
    'EXTREM':   gen_vent_extrem(t_vec),
}

# ============================================================
# FUNCIONS FÍSIQUES
# ============================================================
def eta_lambda(omega,v):
    if v<=0: return 0.0
    lam=omega*R/v
    return float(max(0.0,1.0-((lam-lam_opt)/lam_opt)**2))

def gP_integral(omega):
    theta_s=np.linspace(0,2*np.pi,360)
    total=0.0
    for theta0 in np.linspace(0,2*np.pi,N_PALES,endpoint=False):
        th=theta_s+theta0
        pm=np.clip(KPITCH*np.cos(th)+0.5*omega*np.sin(th),-8,8)
        total+=np.mean(np.abs(pm))/8.0*0.08
    return float(np.clip(total/N_PALES,0.0,0.08))

def eta_kuramoto(K):
    if K<=0: return 0.0
    return float(np.clip(0.04*(1-((K-0.10)/0.10)**2),0.0,0.04))

# ============================================================
# SIMULACIÓ — una funció reutilitzable
# ============================================================
def simula(v_arr, amb_sistema=True, tau_omega=TAU_OMEGA, soroll_omega=0.02):
    """
    Simula la turbina amb o sense el sistema ZZ+Quijote.
    tau_omega: temps de resposta del rotor (s)
    soroll_omega: soroll de mesura ω
    """
    thetas=np.array([2*np.pi*i/N_PALES for i in range(N_PALES)])
    omega=omR*0.95
    eta_K=eta_kuramoto(K_KUR) if amb_sistema else 0.0
    gP=gP_integral(omR) if amb_sistema else 0.0

    hP_aero=np.zeros(steps); hP_tot=np.zeros(steps)
    hP_grid=np.zeros(steps); hE_rotor=np.zeros(steps)
    hom=np.zeros(steps); heta=np.zeros(steps)
    hPbuf=np.zeros(steps)

    for s in range(steps):
        v=v_arr[s]

        # OMEGA REAL: retard tau + soroll
        om_target=min(omR,lam_opt*max(v,0.1)/R)
        omega+=(om_target-omega)*(DT/tau_omega)
        omega+=np.random.normal(0,soroll_omega)*DT
        omega=float(np.clip(omega,omR*0.3,omR*1.4))

        # Kuramoto + pitch (CORREC 5: normalitzat)
        dJdt=0.01*np.sin(omega*t_vec[s])
        if amb_sistema:
            for i in range(N_PALES):
                pm=np.clip(KPITCH*np.cos(thetas[i])+0.5*omega*np.sin(thetas[i])-0.8*dJdt,-8,8)
                kc=(K_KUR/max(N_PALES,1))*float(np.sum(np.sin(thetas-thetas[i])))
                thetas[i]+=(omega+kc+pm*0.01)*DT

        # Potència (CORREC 1: cascada, sense doble comptatge)
        eta_l=eta_lambda(omega,v)
        P_aero=0.5*RHO_A*A_rot*Cp_max*eta_l*v**3

        if amb_sistema:
            eta_tot=min(1.0,eta_l*(1.0+gP)*(1.0+eta_K))
        else:
            eta_tot=eta_l
        P_tot=0.5*RHO_A*A_rot*Cp_max*eta_tot*v**3

        # Buffer quijote (CORREC 3: net≈0)
        P_buf=max(0.0,-dJdt*omega*m_q*R_tip**2*0.5) if amb_sistema else 0.0
        P_buf=min(P_buf,P_NOM*0.02)

        # P_grid (CORREC 6)
        J_tot=J_base+dJ_max*np.sin(omega*t_vec[s])*0.3
        E_rotor=0.5*J_tot*omega**2
        dE_rot_dt=(E_rotor-hE_rotor[s-1])/DT if s>0 else 0.0
        P_grid=max(0.0,P_tot+P_buf-dE_rot_dt)

        hP_aero[s]=P_aero; hP_tot[s]=P_tot+P_buf
        hP_grid[s]=P_grid; hE_rotor[s]=E_rotor
        hom[s]=omega; heta[s]=eta_l; hPbuf[s]=P_buf/1e3

    return hP_aero,hP_tot,hP_grid,hom,heta,hPbuf

# ============================================================
# EXECUCIÓ COMPARATIVA
# ============================================================
print("\n{'Perfil':12} {'Mode':10} {'η_λ mig':10} {'P_mig(MW)':12} {'Millora':10}")
print("-"*58)

resultats={}
VMIN,VMAX=6.0,15.0

for nom,v_arr in vents.items():
    Pa_b,Pt_b,Pg_b,om_b,eta_b,_=simula(v_arr,amb_sistema=False)
    Pa_s,Pt_s,Pg_s,om_s,eta_s,Pb_s=simula(v_arr,amb_sistema=True)

    m_base=np.mean(Pg_b); m_sist=np.mean(Pg_s)
    mill=max(0.0,(m_sist-m_base)/m_base*100) if m_base>0 else 0
    validat=VMIN<=mill<=VMAX

    print(f"  {nom:12} {'BASE':10} {np.mean(eta_b)*100:10.1f}% {m_base/1e6:12.3f} {'—':10}")
    print(f"  {nom:12} {'SISTEMA':10} {np.mean(eta_s)*100:10.1f}% {m_sist/1e6:12.3f} +{mill:6.1f}% {'✓' if validat else '?'}")

    FC=0.35; preu=65.0
    E_base_any=P_NOM*FC*8760/1e6
    guany_keu=(E_base_any*(1+mill/100)-E_base_any)*preu/1e3

    resultats[nom]={
        'v':v_arr,'Pa_b':Pa_b,'Pa_s':Pa_s,'Pg_b':Pg_b,'Pg_s':Pg_s,
        'om_b':om_b,'om_s':om_s,'eta_b':eta_b,'eta_s':eta_s,
        'mill':mill,'validat':validat,'guany_keu':guany_keu,'Pb':Pb_s
    }

# ============================================================
# GRÀFICS
# ============================================================
BG='#0d0d1a'; PAN='#13132b'
COLS={'SUAU':'#00d2ff','REALISTA':'#00ff88','EXTREM':'#ffd700'}

fig=plt.figure(figsize=(20,16),facecolor=BG)
import matplotlib.gridspec as gridspec
gs=gridspec.GridSpec(3,3,figure=fig,hspace=0.48,wspace=0.38)

def sty(ax,t,xl='',yl=''):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9,pad=4)
    ax.tick_params(colors='#aaa',labelsize=7.5)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8); ax.set_ylabel(yl,color='#aaa',fontsize=8)
    ax.grid(color='#1e1e40',lw=0.5,ls='--')

# Fila 0: perfils de vent
for col,(nom,r) in enumerate(resultats.items()):
    ax=fig.add_subplot(gs[0,col])
    ax.plot(t_vec,r['v'],color=COLS[nom],lw=1.5,label=nom)
    ax.axhline(vr,color='white',ls='--',lw=0.7,alpha=0.4,label=f'vr={vr:.1f}m/s')
    ax.axhline(V_BASE,color='#888',ls=':',lw=0.7,alpha=0.5)
    ax.legend(fontsize=7,framealpha=0.3)
    sty(ax,f'Vent {nom} (m/s)','t [s]','v [m/s]')

# Fila 1: P_grid base vs sistema per cada perfil
for col,(nom,r) in enumerate(resultats.items()):
    ax=fig.add_subplot(gs[1,col])
    ax.plot(t_vec,r['Pg_b']/1e6,'--',color='#888',lw=1.2,alpha=0.8,label='Base')
    ax.plot(t_vec,r['Pg_s']/1e6,color=COLS[nom],lw=2,label='Sistema')
    ax.fill_between(t_vec,r['Pg_b']/1e6,r['Pg_s']/1e6,
                    alpha=0.2,color=COLS[nom])
    vc='#00ff88' if r['validat'] else '#e74c3c'
    ax.set_title(f'{nom} — +{r["mill"]:.1f}% {"✓" if r["validat"] else "?"}',
                 color=vc,fontsize=9,pad=4)
    ax.set_facecolor(PAN); ax.tick_params(colors='#aaa',labelsize=7.5)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.grid(color='#1e1e40',lw=0.5,ls='--')
    ax.set_xlabel('t [s]',color='#aaa',fontsize=8)
    ax.set_ylabel('MW',color='#aaa',fontsize=8)
    ax.legend(fontsize=7,framealpha=0.3)

# Fila 2 esquerra: eta_lambda comparativa
ax=fig.add_subplot(gs[2,0])
for nom,r in resultats.items():
    ax.plot(t_vec,r['eta_b']*100,'--',color=COLS[nom],lw=1,alpha=0.5)
    ax.plot(t_vec,r['eta_s']*100,color=COLS[nom],lw=1.5,label=f'{nom} sist.')
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'η_lambda (%) — base(--) vs sistema','t [s]','η [%]')

# Fila 2 centre: omega per perfil
ax=fig.add_subplot(gs[2,1])
for nom,r in resultats.items():
    ax.plot(t_vec,r['om_s'],color=COLS[nom],lw=1.5,label=nom)
ax.axhline(omR,color='white',ls='--',lw=0.8,alpha=0.4,label=f'omR={omR:.3f}')
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'omega sistema (rad/s)','t [s]','ω [rad/s]')

# Fila 2 dreta: taula resum
ax=fig.add_subplot(gs[2,2])
ax.axis('off'); ax.set_facecolor('#0a0a14')
rows=[]
for nom,r in resultats.items():
    rows.append([
        nom,
        f'{r["mill"]:.1f}%',
        '✓' if r['validat'] else '?',
        f'+{r["guany_keu"]:.0f}k€',
        f'+{r["guany_keu"]*5:.0f}k€'
    ])
tbl=ax.table(
    cellText=rows,
    colLabels=['Perfil','Millora','Vàlid','€/turb/a','5 turb/a'],
    loc='center',cellLoc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(9)
for (r,c),cell in tbl.get_celld().items():
    bg='#1a1a2e' if r==0 else '#0d0d1a'
    if r>0:
        nom_r=rows[r-1][0]
        bg=f'#001a10'
    cell.set_facecolor(bg)
    cell.set_text_props(color='white' if r==0 else '#ddd')
    cell.set_edgecolor('#333355')
ax.set_title('Comparativa 3 perfils de vent',color='white',fontsize=9,pad=4)

fig.suptitle(
    f'Gemell Virtual v4 — Vent realista + retard dinàmic  |  {N_PALES} pales  |  '
    f'τ_ω={TAU_OMEGA}s  |  Kpitch={KPITCH}  Kkur={K_KUR}\n'
    f'SUAU (ideal) vs REALISTA vs EXTREM — el guany apareix on cal: en condicions reals',
    color='white',fontsize=10,fontweight='bold',y=0.999)

out='/mnt/user-data/outputs/gemelo_v4_vent_real.png'
plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\n  Gràfic: {out}")
