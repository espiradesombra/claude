"""
GEMELL VIRTUAL v5 — Quijote com a control actiu d'inèrcia
Víctor Manzanares Alberola — EPSA UPV Alcoi

CANVIS vs v4 (patch LeeChat):
  1. J(t) controlat per error_lambda i dω/dt (no sinusoïdal)
  2. Dinàmica omega real: J·ω̇ = T_aero - T_gen - ω·J_dot
  3. Terme -ω·J_dot implementat (ací viu Quijote)
  4. Control anticipatiu: ràfega→J↑, calmant→J↓
  5. Mesura variància P_grid (no sols mitja)
  6. Pitch asimètric: BETA suavitza 0° sense igualar 180°
  7. 3 tests obligatoris: ràfega / calmant / turbulència OU
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
BETA      = 0.5    # asimetria pitch: 0=simètric, 1=màx asimètric
K_KUR     = 0.10
RHO_FLUID = 3386.0
D_CANAL   = 0.05
V_BASE    = 11.4
DT        = 0.05

# Control Quijote
Kq1       = 2.0e5  # guany error_lambda → J_dot
Kq2       = 1.0e5  # guany amortiment  → J_dot
J_DOT_MAX = 5e4    # kg·m²/s  (límit físic actuador)
TAU_ACT   = 1.0    # s (retard actuador fluid)

# Turbina
R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.45; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3)
omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0
A_can=np.pi*(D_CANAL/2)**2
J_base=3.54e7

# Límits J
m_q=RHO_FLUID*A_can*(R_tip-R_hub)
dJ_max=N_PALES*m_q*(R_tip**2-R_hub**2)
J_MIN=J_base; J_MAX=J_base+dJ_max

print(f"omR={omR:.3f}rad/s | vr={vr:.2f}m/s")
print(f"J_base={J_base:.2e} | J_max={J_MAX:.2e} | ΔJ={dJ_max:.2e} ({dJ_max/J_base*100:.0f}%)")

# ============================================================
# FUNCIONS
# ============================================================
def Cp_dinamica(lam, pitch_bias=0.0):
    """
    Cp(λ, pitch) dinàmic — el pitch mou λ_opt efectiu.
    pitch_bias>0 → λ_opt efectiu es desplaça lleugerament.
    Gaussiana al voltant de λ_opt.
    """
    lam_opt_eff = lam_opt + 0.3*pitch_bias  # pitch mou λ_opt
    sigma_lam   = 2.5                        # amplada corba Cp
    Cp = Cp_max * np.exp(-((lam-lam_opt_eff)/sigma_lam)**2)
    return float(np.clip(Cp, 0.0, Cp_max))

def pitch_asimetric(theta_i, omega, dJdt):
    """
    Pitch asimètric (LeeChat + Víctor):
    - 180° (avall): pitch+ màxim
    - 0°  (amunt):  pitch- parcial (BETA controla quant)
    """
    pm_base = KPITCH*np.cos(theta_i) + 0.5*omega*np.sin(theta_i) - 0.8*dJdt
    # Asimetria: redueix la component negativa a 0°
    if pm_base < 0:
        pm_base *= (1.0 - BETA * max(0.0, -np.cos(theta_i)))
    return float(np.clip(pm_base, -8.0, 8.0))

def eta_kuramoto(K):
    if K<=0: return 0.0
    return float(np.clip(0.04*(1-((K-0.10)/0.10)**2), 0.0, 0.04))

def low_pass(x_new, x_old, tau, dt):
    alpha = dt/(tau+dt)
    return alpha*x_new + (1-alpha)*x_old

# ============================================================
# PERFILS DE VENT — 3 tests obligatoris
# ============================================================
T=120.0; steps=int(T/DT); t_vec=np.arange(steps)*DT

def gen_rafega(t_vec):
    """Test 1: ràfega brusca a t=30s"""
    v=np.full(len(t_vec), V_BASE)
    for i,t in enumerate(t_vec):
        v[i]=V_BASE
        if 30<t<35: v[i]+=6.0*np.sin(np.pi*(t-30)/5)   # ràfega brusca 6m/s
        v[i]=max(3.0,v[i])
    return v

def gen_calmant(t_vec):
    """Test 2: calmant sobtat a t=40s"""
    v=np.full(len(t_vec), V_BASE)
    for i,t in enumerate(t_vec):
        v[i]=V_BASE
        if 40<t<55: v[i]=max(3.0,V_BASE-7.0)            # caiguda -7m/s brusca
        if t>55:    v[i]=V_BASE+2.0*np.sin(2*np.pi*t/30) # recuperació
        v[i]=max(3.0,v[i])
    return v

def gen_turbulencia(t_vec):
    """Test 3: turbulència OU forta (IEC classe A)"""
    ou=np.zeros(len(t_vec))
    for i in range(1,len(t_vec)):
        ou[i]=ou[i-1]-0.8*ou[i-1]*DT+1.8*np.sqrt(DT)*np.random.randn()
    v=np.zeros(len(t_vec))
    for i,t in enumerate(t_vec):
        v[i]=V_BASE+2.0*np.sin(2*np.pi*t/40)+ou[i]
        if 60<t<80: v[i]-=4.0
        v[i]=max(2.0,v[i])
    return v

tests={
    'Ràfega brusca':   gen_rafega(t_vec),
    'Calmant sobtat':  gen_calmant(t_vec),
    'Turbulència OU':  gen_turbulencia(t_vec),
}

# ============================================================
# SIMULACIÓ PRINCIPAL
# ============================================================
def simula(v_arr, amb_sistema=True):
    # Estat inicial
    omega    = omR*0.95
    omega_old= omega
    J        = J_base
    J_dot    = 0.0
    J_dot_f  = 0.0  # filtrat (retard actuador)
    thetas   = np.array([2*np.pi*i/N_PALES for i in range(N_PALES)])
    eta_K    = eta_kuramoto(K_KUR) if amb_sistema else 0.0
    pitch_bias=0.0

    # Historials
    hPaero=np.zeros(steps); hPtot=np.zeros(steps); hPgrid=np.zeros(steps)
    hom=np.zeros(steps);    hJ=np.zeros(steps);    hJdot=np.zeros(steps)
    hPbuf=np.zeros(steps);  hlam=np.zeros(steps)

    for s in range(steps):
        t=s*DT; v=v_arr[s]

        # λ actual i error
        lam_act = omega*R/max(v,0.1)
        error_lam = lam_opt - lam_act  # >0 → rotor lent → J↓

        # dω/dt (amortiment)
        domega_dt = (omega-omega_old)/DT if s>0 else 0.0
        omega_old = omega

        # =========================================
        # PATCH LEEHAT: J controlat (no sinusoïdal)
        # =========================================
        if amb_sistema:
            # Llei de control: error_lambda + amortiment
            J_dot_target = (Kq1*(-error_lam) + Kq2*(-domega_dt))
            J_dot_target = float(np.clip(J_dot_target, -J_DOT_MAX, J_DOT_MAX))
            # Retard actuador fluid (low-pass tau=1s)
            J_dot_f = low_pass(J_dot_target, J_dot_f, TAU_ACT, DT)
            J = float(np.clip(J + J_dot_f*DT, J_MIN, J_MAX))
        else:
            J=J_base; J_dot_f=0.0

        # =========================================
        # DINÀMICA OMEGA REAL (patch LeeChat):
        # J·ω̇ = T_aero - T_gen - ω·J_dot
        # =========================================
        pitch_bias = float(np.mean([KPITCH*np.cos(th) for th in thetas])) if amb_sistema else 0.0
        Cp = Cp_dinamica(lam_act, pitch_bias)
        P_aero = 0.5*RHO_A*A_rot*Cp*v**3

        # Generador: segueix P_nom fins al límit
        P_gen  = min(P_aero, P_NOM*1.05)
        T_aero = P_aero/max(omega,0.1)
        T_gen  = P_gen /max(omega,0.1)

        # EQ. MOVIMENT AMB TERME -ω·J_dot
        omega_dot = (T_aero - T_gen - omega*J_dot_f) / max(J, J_base*0.5)
        omega = float(np.clip(omega + omega_dot*DT, omR*0.3, omR*1.5))

        # Kuramoto + pitch asimètric
        dJdt_sig = J_dot_f/max(dJ_max,1.0)  # normalitzat [-1,1]
        if amb_sistema:
            for i in range(N_PALES):
                pm = pitch_asimetric(thetas[i], omega, dJdt_sig)
                kc = (K_KUR/N_PALES)*float(np.sum(np.sin(thetas-thetas[i])))
                thetas[i]+=(omega+kc+pm*0.01)*DT

        # P_grid real
        dE_rot_dt = J_dot_f*0.5*omega**2 + J*omega*omega_dot  # dE/dt = d(½Jω²)/dt
        P_grid = max(0.0, P_gen - dE_rot_dt)

        hPaero[s]=P_aero; hPtot[s]=P_gen; hPgrid[s]=P_grid
        hom[s]=omega; hJ[s]=J; hJdot[s]=J_dot_f
        hPbuf[s]=(J-J_base)/dJ_max*100  # % J usat
        hlam[s]=lam_act

    return hPaero,hPtot,hPgrid,hom,hJ,hJdot,hPbuf,hlam

# ============================================================
# EXECUTA ELS 3 TESTS
# ============================================================
print()
print(f"{'Test':20} {'η_P_grid':10} {'Millora':10} {'Δstd':10} {'Valid?'}")
print("-"*58)

resultats={}
COLS={'Ràfega brusca':'#ffd700','Calmant sobtat':'#00d2ff','Turbulència OU':'#ff6b35'}

for nom,v_arr in tests.items():
    Pb,Pt_b,Pg_b,om_b,J_b,Jd_b,buf_b,lam_b = simula(v_arr, False)
    Ps,Pt_s,Pg_s,om_s,J_s,Jd_s,buf_s,lam_s = simula(v_arr, True)

    mill  = max(0.0,(np.mean(Pg_s)-np.mean(Pg_b))/np.mean(Pg_b)*100)
    std_b = np.std(Pg_b); std_s = np.std(Pg_s)
    dstd  = (std_b-std_s)/std_b*100  # reducció variància (+ = millor)
    valid = mill>=2.0 and dstd>=5.0

    print(f"  {nom:20} {np.mean(Pg_s)/1e6:.3f}MW {mill:+9.1f}% {dstd:+9.1f}% {'✓' if valid else '?'}")

    resultats[nom]={
        'v':v_arr,'Pg_b':Pg_b,'Pg_s':Pg_s,
        'om_b':om_b,'om_s':om_s,
        'J_s':J_s,'Jd_s':Jd_s,
        'lam_b':lam_b,'lam_s':lam_s,
        'mill':mill,'dstd':dstd,'valid':valid,
        'std_b':std_b,'std_s':std_s,
        'col':COLS[nom],'buf':buf_s,
    }

# ============================================================
# GRÀFICS
# ============================================================
BG='#0d0d1a'; PAN='#13132b'
fig=plt.figure(figsize=(22,18),facecolor=BG)
gs=gridspec.GridSpec(4,3,figure=fig,hspace=0.52,wspace=0.38)

def sty(ax,t,xl='t [s]',yl=''):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=8.5,pad=4)
    ax.tick_params(colors='#aaa',labelsize=7.5)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8); ax.set_ylabel(yl,color='#aaa',fontsize=8)
    ax.grid(color='#1e1e40',lw=0.5,ls='--')

# Fila 0: Vent dels 3 tests
for col,(nom,r) in enumerate(resultats.items()):
    ax=fig.add_subplot(gs[0,col])
    ax.plot(t_vec,r['v'],color=r['col'],lw=1.5)
    ax.axhline(vr,color='white',ls='--',lw=0.7,alpha=0.4,label=f'vr={vr:.1f}')
    ax.axhline(V_BASE,color='#888',ls=':',lw=0.5,alpha=0.5)
    ax.legend(fontsize=7,framealpha=0.3)
    sty(ax,f'Vent — {nom}','t [s]','v [m/s]')

# Fila 1: P_grid base vs sistema
for col,(nom,r) in enumerate(resultats.items()):
    ax=fig.add_subplot(gs[1,col])
    ax.plot(t_vec,r['Pg_b']/1e6,'--',color='#888',lw=1.2,alpha=0.8,label='Base')
    ax.plot(t_vec,r['Pg_s']/1e6,color=r['col'],lw=2,label='Quijote ctrl')
    ax.fill_between(t_vec,r['Pg_b']/1e6,r['Pg_s']/1e6,alpha=0.2,color=r['col'])
    vc='#00ff88' if r['valid'] else '#e74c3c'
    ax.set_title(
        f'{nom}\n+{r["mill"]:.1f}% E  |  Δstd={r["dstd"]:+.1f}% {"✓" if r["valid"] else "?"}',
        color=vc,fontsize=8.5,pad=4)
    ax.set_facecolor(PAN); ax.tick_params(colors='#aaa',labelsize=7.5)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.grid(color='#1e1e40',lw=0.5,ls='--')
    ax.set_xlabel('t [s]',color='#aaa',fontsize=8)
    ax.set_ylabel('MW',color='#aaa',fontsize=8)
    ax.legend(fontsize=7,framealpha=0.3)

# Fila 2: J(t) controlat + lambda
for col,(nom,r) in enumerate(resultats.items()):
    ax=fig.add_subplot(gs[2,col]); ax2=ax.twinx()
    ax.plot(t_vec,r['J_s']/1e6,color=r['col'],lw=2,label='J(t) (MJ·m²)')
    ax.axhline(J_base/1e6,color='white',ls='--',lw=0.7,alpha=0.4,label='J_base')
    ax.axhline(J_MAX/1e6, color='#888',ls=':',lw=0.7,alpha=0.4,label='J_max')
    ax2.plot(t_vec,r['lam_s'],color='#ff6b35',lw=1.2,alpha=0.7,label='λ(t)')
    ax2.axhline(lam_opt,color='#ffd700',ls='--',lw=0.7,alpha=0.5,label=f'λ_opt={lam_opt}')
    ax.set_ylabel('J [MJ·m²]',color=r['col'],fontsize=8)
    ax2.set_ylabel('λ',color='#ff6b35',fontsize=8)
    ax.tick_params(colors=r['col']); ax2.tick_params(colors='#ff6b35')
    ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--')
    ax.set_xlabel('t [s]',color='#aaa',fontsize=8)
    ax.set_title(f'J(t) controlat + λ — {nom}',color='white',fontsize=8.5,pad=4)
    ax.legend(fontsize=6.5,framealpha=0.3,loc='upper left')
    ax2.legend(fontsize=6.5,framealpha=0.3,loc='upper right')

# Fila 3: omega + taula resum
for col,(nom,r) in enumerate(resultats.items()):
    if col < 2:
        ax=fig.add_subplot(gs[3,col])
        ax.plot(t_vec,r['om_b'],'--',color='#888',lw=1.2,alpha=0.7,label='Base ω')
        ax.plot(t_vec,r['om_s'],color=r['col'],lw=2,label='Quijote ω')
        ax.axhline(omR,color='white',ls='--',lw=0.7,alpha=0.4)
        sty(ax,f'omega (rad/s) — {nom}','t [s]','ω [rad/s]')
        ax.legend(fontsize=7,framealpha=0.3)

# Taula resum
ax=fig.add_subplot(gs[3,2]); ax.axis('off'); ax.set_facecolor('#0a0a14')
rows=[]
for nom,r in resultats.items():
    rows.append([
        nom[:14],
        f'{r["mill"]:+.1f}%',
        f'{r["dstd"]:+.1f}%',
        f'{r["std_b"]/1e6:.2f}→{r["std_s"]/1e6:.2f}MW',
        '✓' if r['valid'] else '?'
    ])
tbl=ax.table(
    cellText=rows,
    colLabels=['Test','ΔE','Δstd','std_base→sist','OK'],
    loc='center',cellLoc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5)
for (r,c),cell in tbl.get_celld().items():
    bg='#1a1a2e' if r==0 else '#0d0d1a'
    if r>0 and c==4:
        bg='#003322' if rows[r-1][4]=='✓' else '#1a1000'
    cell.set_facecolor(bg)
    cell.set_text_props(color='white' if r==0 else '#ddd')
    cell.set_edgecolor('#333355')
ax.set_title('Resum 3 tests obligatoris',color='white',fontsize=9,pad=4)

fig.suptitle(
    f'Gemell Virtual v5 — Quijote Control Actiu  |  {N_PALES} pales  |  '
    f'Kpitch={KPITCH} β={BETA} K_kur={K_KUR}  Kq1={Kq1:.0e} Kq2={Kq2:.0e}\n'
    f'J·ω̇ = T_aero − T_gen − ω·J_dot  |  '
    f'Criteri: ΔE>2% AND Δstd>5% → Quijote funciona',
    color='white',fontsize=10,fontweight='bold',y=0.999)

out='/mnt/user-data/outputs/gemelo_v5_quijote_control.png'
plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\n  Gràfic: {out}")
