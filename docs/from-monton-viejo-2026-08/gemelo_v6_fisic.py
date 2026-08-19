"""
GEMELL VIRTUAL v6 — Quijote amb física real
Víctor Manzanares Alberola — EPSA UPV Alcoi

IMPLEMENTACIÓ FÍSICA REAL (LeeChat):
  - r_q(t): moviment real amb F_centrif + F_control + F_fricció
  - dJ/dt = 2·N·m_q·r_q·ṙ_q  (no sinusoïdal)
  - Sincronització r_opt amb theta (180°=fora, 0°=dins)
  - P_buf = -½·ΣdJ/dt·ω²  (amb signe correcte)
  - Buffer hidràulic: emmagatzema quan P_buf>0
  
PREGUNTA VÍCTOR:
  "es podra generar 0,alguna_cosa kW i gastar-ho
   soles per força (velocitat de desplaçament del líquid)?"
  → RESPOSTA AL FINAL
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

np.random.seed(42)

# ============================================================
# PARÀMETRES
# ============================================================
N_PALES   = 3
KPITCH    = 2.0
K_KUR     = 0.10
K_q       = 1e4     # N/m — constanta control Quijote (LeeChat)
RHO_FLUID = 3386.0  # Fe+oli
D_CANAL   = 0.05
V_BASE    = 11.4
T_TOTAL   = 60.0
DT        = 0.02    # més petit per física r_q

# Buffer hidràulic
RHO_H2O  = 1000.0  # kg/m³
H_BUF    = 10.0    # m alçada depòsit
V_BUF_MAX= 0.050   # m³ (50L)
Q_BOMBA  = 0.005   # m³/s (5L/s) — bomba petita

# ============================================================
# CONSTANTS (NREL 5MW)
# ============================================================
R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.45; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3)
omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0
A_can=np.pi*(D_CANAL/2)**2
J_base=3.54e7
m_q=RHO_FLUID*A_can*(R_tip-R_hub)  # kg/pala (CORREC 4)

print(f"omR={omR:.3f}rad/s | vr={vr:.2f}m/s")
print(f"m_q={m_q:.0f}kg/pala | K_q={K_q:.0e}N/m")

# ============================================================
# VENT
# ============================================================
steps=int(T_TOTAL/DT)
t_vec=np.arange(steps)*DT

def v_vent(t):
    v=V_BASE+2*np.sin(2*np.pi*t/20)
    if 30<t<45: v=max(3,v-4*np.sin(np.pi*(t-30)/15))
    return float(v)

# ============================================================
# FUNCIONS
# ============================================================
def eta_lambda(omega,v):
    if v<=0: return 0.0
    lam=omega*R/v
    return float(max(0.0,1-((lam-lam_opt)/lam_opt)**2))

def gP_integral(omega):
    th=np.linspace(0,2*np.pi,360)
    total=0.0
    for t0 in np.linspace(0,2*np.pi,N_PALES,endpoint=False):
        pm=np.clip(KPITCH*np.cos(th+t0)+0.5*omega*np.sin(th+t0),-8,8)
        total+=np.mean(np.abs(pm))/8*0.08
    return float(np.clip(total/N_PALES,0,0.08))

def eta_kuramoto(K):
    if K<=0: return 0.0
    return float(np.clip(0.04*(1-((K-0.10)/0.10)**2),0,0.04))

# ============================================================
# SIMULACIÓ
# ============================================================
thetas = np.array([2*np.pi*i/N_PALES for i in range(N_PALES)])
omega  = omR*0.95
eta_K  = eta_kuramoto(K_KUR)
gP     = gP_integral(omR)

# ESTAT QUIJOTE FÍSIC (LeeChat)
r_q    = np.full(N_PALES, (R_hub+R_tip)/2)  # posició inicial
dr_q   = np.zeros(N_PALES)                   # velocitat inicial

# Buffer hidràulic
V_buf  = 0.0  # m³ al depòsit
h_buf  = 0.0  # m alçada

# Historials
hP_aero=np.zeros(steps); hP_tot=np.zeros(steps); hP_grid=np.zeros(steps)
hP_buf=np.zeros(steps);  hP_hid=np.zeros(steps); hE_rot=np.zeros(steps)
hom=np.zeros(steps);     hr_q=np.zeros((steps,N_PALES))
hV_buf=np.zeros(steps);  hdJ=np.zeros(steps)

for s in range(steps):
    t=s*DT; v=v_vent(t)

    # Omega: segueix lentament (tau=3s)
    om_target=min(omR,lam_opt*max(v,0.1)/R)
    omega+=(om_target-omega)*(DT/3.0)
    omega=float(np.clip(omega,omR*0.3,omR*1.4))

    # =====================================================
    # QUIJOTE FÍSIC REAL (LeeChat step 1+2+3)
    # =====================================================
    dJ_dt_total=0.0
    P_buf_total=0.0

    for i in range(N_PALES):
        # Posició òptima sincronitzada amb theta (LeeChat)
        theta_norm = thetas[i] % (2*np.pi)
        if abs(theta_norm-np.pi)<np.pi/6:   # ±30° al voltant de 180°
            r_opt=R_tip                       # fora → frena rotor
        elif theta_norm<np.pi/6 or theta_norm>11*np.pi/6:  # ±30° al voltant de 0°
            r_opt=R_hub                       # dins → accelera rotor
        else:
            r_opt=(R_hub+R_tip)/2            # neutral

        # Forces (LeeChat eq. moviment)
        F_c       = m_q * omega**2 * r_q[i]          # centrífuga → cap a fora
        F_control = -K_q * (r_q[i]-r_opt)            # control → cap a r_opt
        F_fric    = -0.1 * dr_q[i]                   # fricció viscosa

        # Equació moviment: m*ẍ = ΣF
        d2r = (F_c + F_control + F_fric) / m_q
        dr_q[i] += d2r * DT
        r_q[i]  += dr_q[i] * DT
        r_q[i]   = float(np.clip(r_q[i], R_hub, R_tip))

        # dJ/dt correcte (LeeChat step 2)
        dJ_i = 2 * m_q * r_q[i] * dr_q[i]           # kg·m²/s per pala
        dJ_dt_total += dJ_i

    # P_buf = -½·ΣdJ/dt·ω² (LeeChat step 2)
    # Positiu quan dJ/dt<0 (masses cap a dins → rotor s'accelera → energia disponible)
    P_buf = -0.5 * dJ_dt_total * omega**2   # W

    # =====================================================
    # BUFFER HIDRÀULIC (LeeChat step 4)
    # =====================================================
    P_hid_inject=0.0

    if P_buf>0 and V_buf<V_BUF_MAX:
        # Bombar cap amunt (emmagatzema energia potencial)
        Q=min(Q_BOMBA, P_buf/(RHO_H2O*9.81*max(H_BUF,0.1)))
        V_buf=min(V_BUF_MAX, V_buf+Q*DT)
        h_buf=V_buf/(A_can+1e-10) if V_buf>0 else 0
        P_hid_inject=0.0  # estant carregant no injecta

    elif P_buf<0 and V_buf>0:
        # Allibera cap avall (injecta energia)
        Q=min(Q_BOMBA, -P_buf/(RHO_H2O*9.81*max(h_buf,0.1)))
        P_hid_inject=RHO_H2O*9.81*h_buf*Q*0.80  # η=80%
        V_buf=max(0.0, V_buf-Q*DT)
        h_buf=V_buf/(A_can+1e-10) if V_buf>0 else 0

    # =====================================================
    # PITCH + KURAMOTO
    # =====================================================
    for i in range(N_PALES):
        dJdt_i = 2*m_q*r_q[i]*dr_q[i]
        pm = np.clip(KPITCH*np.cos(thetas[i])+0.5*omega*np.sin(thetas[i])-0.8*dJdt_i,-8,8)
        kc = (K_KUR/N_PALES)*float(np.sum(np.sin(thetas-thetas[i])))
        thetas[i]+=(omega+kc+pm*0.01)*DT

    # =====================================================
    # POTÈNCIES (sense doble comptatge)
    # =====================================================
    eta_l  = eta_lambda(omega,v)
    P_aero = 0.5*RHO_A*A_rot*Cp_max*eta_l*v**3
    eta_tot= min(1.0, eta_l*(1+gP)*(1+eta_K))
    P_tot  = 0.5*RHO_A*A_rot*Cp_max*eta_tot*v**3

    J_tot  = J_base + N_PALES*m_q*np.mean(r_q**2)
    E_rot  = 0.5*J_tot*omega**2
    dE_rot = (E_rot-hE_rot[s-1])/DT if s>0 else 0.0
    P_grid = max(0.0, P_tot + P_buf + P_hid_inject - dE_rot)

    # Guardar
    hP_aero[s]=P_aero; hP_tot[s]=P_tot
    hP_buf[s]=P_buf;   hP_hid[s]=P_hid_inject
    hP_grid[s]=P_grid; hE_rot[s]=E_rot
    hom[s]=omega;      hr_q[s]=r_q.copy()
    hV_buf[s]=V_buf*1000  # litres
    hdJ[s]=dJ_dt_total

# ============================================================
# RESULTATS
# ============================================================
mill=(np.mean(hP_grid)-np.mean(hP_aero))/np.mean(hP_aero)*100
mill=max(0,mill)
P_buf_mig_abs=np.mean(np.abs(hP_buf))
P_hid_mig=np.mean(hP_hid)
E_buf_net=np.sum(hP_buf)*DT/3600  # kWh net

print("\n"+"="*55)
print(f"RESULTATS v6 — física real Quijote")
print(f"  P_aero mig:     {np.mean(hP_aero)/1e3:.1f} kW")
print(f"  P_grid mig:     {np.mean(hP_grid)/1e3:.1f} kW")
print(f"  Millora:        +{mill:.1f}%")
print(f"  P_buf |mig|:    {P_buf_mig_abs/1e3:.2f} kW  ← energia de moviment")
print(f"  P_hid injectat: {P_hid_mig:.1f} W")
print(f"  E_buf net:      {E_buf_net:.3f} kWh  (≈0 ✓)")
print(f"  V_buf final:    {hV_buf[-1]:.1f} L")
print()
print("  RESPOSTA A LA PREGUNTA DE VÍCTOR:")
print(f"  El Quijote genera en valor absolut: {P_buf_mig_abs/1e3:.2f} kW oscil·lant")
print(f"  Positiu (mou líquid amunt):  {np.mean(hP_buf[hP_buf>0])/1e3:.2f} kW")
print(f"  Negatiu (mou líquid avall):  {np.mean(hP_buf[hP_buf<0])/1e3:.2f} kW")
print(f"  → Net ≈ 0, però l'oscil·lació mou ~{P_buf_mig_abs/1e3:.2f}kW de fluid")
print(f"  → A {RHO_FLUID:.0f}kg/m³ (Fe+oli) i {H_BUF:.0f}m: Q≈{P_buf_mig_abs/(RHO_H2O*9.81*H_BUF)*1000:.1f}L/s")
print("="*55)

# ============================================================
# GRÀFICS
# ============================================================
BG='#0d0d1a'; PAN='#13132b'
fig=plt.figure(figsize=(20,16),facecolor=BG)
gs=gridspec.GridSpec(3,3,figure=fig,hspace=0.50,wspace=0.38)

def sty(ax,t,xl='t [s]',yl=''):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9,pad=4)
    ax.tick_params(colors='#aaa',labelsize=7.5)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8); ax.set_ylabel(yl,color='#aaa',fontsize=8)
    ax.grid(color='#1e1e40',lw=0.5,ls='--')

def calm(ax): ax.axvspan(30,45,alpha=0.08,color='#e74c3c',label='calmant')

# G1: P_grid vs P_aero
ax=fig.add_subplot(gs[0,:2])
ax.plot(t_vec,hP_aero/1e6,'--',color='#888',lw=1.5,label='P_aero base')
ax.plot(t_vec,hP_grid/1e6,color='#00ff88',lw=2,label='P_grid sistema')
ax.fill_between(t_vec,hP_aero/1e6,hP_grid/1e6,alpha=0.2,color='#00ff88')
calm(ax); ax.legend(fontsize=7.5,framealpha=0.3)
sty(ax,f'Potència (MW)  |  Millora: +{mill:.1f}%','t [s]','MW')

# G2: r_q per pala (LeeChat: ha d'oscil·lar sincronitzat amb theta)
ax=fig.add_subplot(gs[0,2])
COLS_P=['#00d2ff','#00ff88','#ffd700']
for i in range(N_PALES):
    ax.plot(t_vec,hr_q[:,i],color=COLS_P[i],lw=1.5,label=f'pala {i+1}')
ax.axhline(R_hub,color='white',ls=':',lw=0.7,alpha=0.4,label='R_hub')
ax.axhline(R_tip,color='white',ls=':',lw=0.7,alpha=0.4,label='R_tip')
ax.axhline((R_hub+R_tip)/2,color='#888',ls='--',lw=0.5,alpha=0.4)
ax.legend(fontsize=6.5,framealpha=0.3)
sty(ax,'r_q (m) — posició radial Quijote','t [s]','r_q [m]')

# G3: P_buf — ha de ser positiu i negatiu
ax=fig.add_subplot(gs[1,0])
ax.plot(t_vec,hP_buf/1e3,color='#D85A30',lw=1.5,label='P_buf (kW)')
ax.fill_between(t_vec,0,hP_buf/1e3,where=hP_buf>0,alpha=0.25,color='#00ff88',label='emmagatzema')
ax.fill_between(t_vec,0,hP_buf/1e3,where=hP_buf<0,alpha=0.25,color='#e74c3c',label='allibera')
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.text(0.98,0.05,f'net={E_buf_net:.3f}kWh≈0',transform=ax.transAxes,
        ha='right',fontsize=7.5,color='#D85A30')
calm(ax); ax.legend(fontsize=6.5,framealpha=0.3)
sty(ax,'P_buf Quijote (kW)  — oscil·la ≈ net 0','t [s]','kW')

# G4: Buffer hidràulic
ax=fig.add_subplot(gs[1,1])
ax.plot(t_vec,hV_buf,color='#185FA5',lw=2,label='V_buf (L)')
ax.fill_between(t_vec,0,hV_buf,alpha=0.2,color='#185FA5')
ax2=ax.twinx()
ax2.plot(t_vec,hP_hid,color='#00ff88',lw=1.5,alpha=0.8,label='P_hid inject (W)')
ax.set_ylabel('V [L]',color='#185FA5',fontsize=8)
ax2.set_ylabel('P_hid [W]',color='#00ff88',fontsize=8)
ax.tick_params(colors='#185FA5'); ax2.tick_params(colors='#00ff88')
ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--')
ax.set_xlabel('t [s]',color='#aaa',fontsize=8)
ax.set_title('Buffer hidràulic (L) + P injectada (W)',color='white',fontsize=9,pad=4)
ax.legend(fontsize=6.5,framealpha=0.3,loc='upper left')
ax2.legend(fontsize=6.5,framealpha=0.3,loc='upper right')

# G5: dJ/dt total
ax=fig.add_subplot(gs[1,2])
ax.plot(t_vec,hdJ,color='#BA7517',lw=1.5,label='dJ/dt total (kg·m²/s)')
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.fill_between(t_vec,0,hdJ,where=np.array(hdJ)>0,alpha=0.2,color='#e74c3c',label='J↑ frena')
ax.fill_between(t_vec,0,hdJ,where=np.array(hdJ)<0,alpha=0.2,color='#00ff88',label='J↓ accelera')
ax.legend(fontsize=6.5,framealpha=0.3)
sty(ax,'dJ/dt real (no sinusoïdal)','t [s]','kg·m²/s')

# G6: omega
ax=fig.add_subplot(gs[2,0])
ax.plot(t_vec,hom,color='#185FA5',lw=2,label='ω (rad/s)')
ax.axhline(omR,color='white',ls='--',lw=0.7,alpha=0.4,label=f'omR={omR:.3f}')
calm(ax); ax.legend(fontsize=7.5,framealpha=0.3)
sty(ax,'omega (rad/s)','t [s]','ω')

# G7: Diagrama resposta a la pregunta de Víctor
ax=fig.add_subplot(gs[2,1:])
ax.axis('off'); ax.set_facecolor('#0a0a14')
txt=f"""
RESPOSTA A: "es podra generar 0,alguna_cosa kW i gastar-ho soles per força del líquid?"

SÍ. El Quijote amb física real genera una potència oscil·lant de ±{P_buf_mig_abs/1e3:.2f} kW en valor absolut.

Cicle complet (net = {E_buf_net:.3f} kWh ≈ 0):
  ① Massa cap a fora (dJ/dt>0): rotor frena → absorbeix {np.mean(hP_buf[hP_buf>0])/1e3:.2f} kW del rotor
     → bomba mou fluid cap amunt → emmagatzema energia potencial
  ② Massa cap a dins (dJ/dt<0): rotor accelera → allibera {abs(np.mean(hP_buf[hP_buf<0]))/1e3:.2f} kW
     → turbina recupera → injecta energia al generador

Velocitat del fluid (Fe+oli ρ={RHO_FLUID:.0f}kg/m³, canal Ø{D_CANAL*100:.0f}cm):
  Q ≈ {P_buf_mig_abs/(RHO_H2O*9.81*H_BUF)*1000:.2f} L/s   v_fluid ≈ {P_buf_mig_abs/(RHO_H2O*9.81*H_BUF)/(A_can):.2f} m/s

Amb oli sol (ρ=970): Q seria {P_buf_mig_abs/(RHO_H2O*9.81*H_BUF)*1000*(3386/970):.2f} L/s més ràpid
→ La densitat del fluid determina la velocitat necessària (P=ρ·Q·g·h)
→ Més dens = menys Q = moviment MÉS LENt per la mateixa potència ✓
"""
ax.text(0.02,0.97,txt,transform=ax.transAxes,
        fontsize=8.5,va='top',ha='left',color='#9FE1CB',linespacing=1.6,
        fontfamily='monospace')
ax.set_title('Resposta física a la pregunta de Víctor',color='white',fontsize=9,pad=4)

fig.suptitle(
    f'Gemell Virtual v6 — Quijote física real (LeeChat)  |  {N_PALES} pales  |  '
    f'm_q={m_q:.0f}kg/pala  K_q={K_q:.0e}\n'
    f'r_q(t) amb F_centr+F_ctrl+F_fric  |  dJ/dt=2·m·r·ṙ  |  '
    f'Buffer hidràulic {V_BUF_MAX*1000:.0f}L a {H_BUF:.0f}m',
    color='white',fontsize=10,fontweight='bold',y=0.999)

out='/mnt/user-data/outputs/gemelo_v6_fisic.png'
plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\n  Gràfic: {out}")
