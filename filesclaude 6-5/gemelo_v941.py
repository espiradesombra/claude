"""
GEMELL VIRTUAL v9.4.1 — Model consolidat per al paper
Victor Manzanares Alberola — EPSA UPV Alcoi

MILLORES vs v9.4 (LeeChat):
  1. K_q = 1e5 (de 5e4) — retorn cap a r_opt mes rapid
  2. r_opt amb compensacio F_c: r_opt += K_comp * F_c/(K_q*(Rtip-Rhub))
  3. KQ = 0.6 (de 0.8) — menys competencia pitch vs quijote
  4. K_KUR adaptatiu: K0 + K1*(omR-om)/omR + K2*|valley|
  5. Valley filtrada: low-pass tau=0.1s
  6. Maquina d'estats: STABLE / VALLEY / RECOVERY

OBJECTIU: E_net tendeix a 0, P_gen > 100W, comparativa 3p vs 7p
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

np.random.seed(42)

# ============================================================
# PARAMETRES v9.4.1
# ============================================================
R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.482; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3); omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0; J_base=3.54e7

# Control Quijote v9.4.1
K_q    = 1e5    # N/m  (MILLORA 1: de 5e4)
K_comp = 0.1    # (MILLORA 2: compensacio F_c)
c_fric = 80.0   # kg/s fricció viscosa
DR_MAX = 0.5    # m/s

# Control ZYPYZAPE
KPITCH = 2.0; KOMEGA = 0.5; KQ = 0.6  # (MILLORA 3: de 0.8)

# Kuramoto adaptatiu (MILLORA 4)
K0_kur=0.10; K1_kur=0.05; K2_kur=0.01

# Valley filtre (MILLORA 5)
TAU_VAL=0.1

# Maquina d'estats (MILLORA 6)
THRESH=0.1
STABLE=0; VALLEY=1; RECOVERY=2

# Fluid Fe+oli
RHO_FL=3386.0; D_CANAL=0.05
A_can=np.pi*(D_CANAL/2)**2

# Buffer hidraulic
H_BUF=20.0; V_BUF_MAX=0.10; Q_BOMBA=0.010
P_MAX_ACC=20e6; LLINDAR=20.0; ETA_GEN=0.85

# Simulacio
V_BASE=11.4; T_TOT=80.0; DT=0.05
steps=int(T_TOT/DT); t_vec=np.arange(steps)*DT

print(f"v9.4.1: omR={omR:.4f} | K_q={K_q:.0e} | KQ={KQ} | K_comp={K_comp}")
print(f"omega*R_hub = {omR*R_hub:.2f} m/s >> DR_MAX={DR_MAX}")

# Vent OU
ou=np.zeros(steps)
for i in range(1,steps):
    ou[i]=ou[i-1]-0.5*ou[i-1]*DT+0.8*np.sqrt(DT)*np.random.randn()
v_arr=np.array([max(3, V_BASE+2*np.sin(2*np.pi*s*DT/20)+ou[s]*0.6
                    -(5*np.sin(np.pi*(s*DT-35)/15) if 35<s*DT<50 else 0))
                for s in range(steps)])

# ============================================================
# FUNCIONS
# ============================================================
def eta_lam(om,v):
    if v<=0: return 0.0
    return float(max(0,1-((om*R/v-lam_opt)/lam_opt)**2))

def m_q_f(rho,D,r1,r2):
    return rho*np.pi*(D/2)**2*(r2-r1)

# ============================================================
# SIMULACIO
# ============================================================
def simula(N, label):
    m_q=m_q_f(RHO_FL,D_CANAL,R_hub,R_tip)

    thetas=np.array([2*np.pi*i/N for i in range(N)])
    omega=omR*0.95
    r_q=np.full(N,(R_hub+R_tip)/2)
    dr_q=np.zeros(N)
    V_buf=0.0; P_acc=0.0
    valley_f=0.0  # valley filtrada
    state=STABLE

    hPa=np.zeros(steps); hPg=np.zeros(steps)
    hPb=np.zeros(steps); hPgen=np.zeros(steps)
    hom=np.zeros(steps); hrq=np.zeros((steps,N))
    hVb=np.zeros(steps); hdJ=np.zeros(steps)
    hPacc=np.zeros(steps); hlim=np.zeros(steps)
    hstate=np.zeros(steps); hval=np.zeros(steps)
    hKkur=np.zeros(steps)

    for s in range(steps):
        v=v_arr[s]
        omega+=(min(omR,lam_opt*max(v,0.1)/R)-omega)*(DT/3.0)
        omega=float(np.clip(omega,omR*0.3,omR*1.4))

        dJ=0.0; lim_min=999.0
        for i in range(N):
            th=thetas[i]%(2*np.pi)

            # r_opt suau + MILLORA 2: compensacio F_c
            r_opt0=R_tip-(R_tip-R_hub)*(1+np.cos(th))/2
            F_c=m_q*omega**2*r_q[i]
            r_opt=r_opt0+K_comp*(F_c/(K_q*(R_tip-R_hub)+1e-9))
            r_opt=float(np.clip(r_opt,R_hub,R_tip))

            # MILLORA 6: maquina d'estats modifica r_opt
            if state==VALLEY:
                r_opt=R_tip    # forcar massa cap a fora (carregar buffer)
            elif state==RECOVERY:
                r_opt=R_hub    # forcar massa cap a dins (accelerar rotor)

            F_ct=-K_q*(r_q[i]-r_opt)
            F_fr=-c_fric*dr_q[i]
            dr_q[i]+=(F_c+F_ct+F_fr)/m_q*DT

            # Limit simetric + limit posicional
            dr_q[i]=float(np.clip(dr_q[i],-DR_MAX,DR_MAX))
            lp=omega*r_q[i]; lim_min=min(lim_min,lp)
            dr_q[i]=float(np.clip(dr_q[i],-lp,lp))

            r_q[i]=float(np.clip(r_q[i]+dr_q[i]*DT,R_hub,R_tip))
            dJ+=2*m_q*r_q[i]*dr_q[i]

        P_buf=-0.5*dJ*omega**2

        # MILLORA 5: valley filtrada
        valley_raw=-omega*dJ
        alpha_v=DT/(TAU_VAL+DT)
        valley_f=alpha_v*valley_raw+(1-alpha_v)*valley_f

        # MILLORA 6: transicio d'estats
        if valley_f>THRESH or omega<omR*0.85:
            state=VALLEY
        elif valley_f<-THRESH and omega>omR*0.85:
            state=RECOVERY
        elif abs(valley_f)<=THRESH and omega>=omR*0.92:
            state=STABLE

        # MILLORA 4: K_KUR adaptatiu
        K_kur=K0_kur+K1_kur*max(0,(omR-omega)/omR)+K2_kur*abs(valley_f)
        K_kur=float(np.clip(K_kur,0,0.3))

        # Pitch + Kuramoto
        for i in range(N):
            dJi=2*m_q*r_q[i]*dr_q[i]
            # MILLORA 3: KQ=0.6
            pm=np.clip(KPITCH*np.cos(thetas[i])+KOMEGA*omega*np.sin(thetas[i])-KQ*dJi,-8,8)
            kc=(K_kur/N)*float(np.sum(np.sin(thetas-thetas[i])))
            thetas[i]+=(omega+kc+pm*0.01)*DT

        # Buffer hidraulic + valvules
        Pgen=0.0; Phid=0.0
        if P_buf>LLINDAR and V_buf<V_BUF_MAX:
            P_acc=min(P_acc+P_buf*DT/V_BUF_MAX,P_MAX_ACC)
            Q=min(Q_BOMBA,P_buf/max(P_acc,1e3))
            V_buf=min(V_BUF_MAX,V_buf+Q*DT)
        elif P_buf<-LLINDAR and V_buf>0:
            P_acc=max(P_acc+P_buf*DT/V_BUF_MAX,0)
            Q=min(Q_BOMBA,-P_buf/max(P_acc,1e3))
            Pgen=P_acc*Q*ETA_GEN
            Phid=P_acc*Q*0.82
            V_buf=max(0,V_buf-Q*DT)

        etL=eta_lam(omega,v)
        Pa=0.5*RHO_A*A_rot*Cp_max*etL*v**3
        etaK=float(np.clip(0.04*min(1.2,N/3),0,0.048))
        Pg=max(0,Pa*(1+0.04)*(1+etaK)+max(0,P_buf)+Phid+Pgen)

        hPa[s]=Pa; hPg[s]=Pg; hPb[s]=P_buf
        hPgen[s]=Pgen; hom[s]=omega; hrq[s]=r_q.copy()
        hVb[s]=V_buf*1000; hdJ[s]=dJ; hPacc[s]=P_acc/1e5
        hlim[s]=lim_min; hstate[s]=state
        hval[s]=valley_f; hKkur[s]=K_kur

    return dict(Pa=hPa,Pg=hPg,Pb=hPb,Pgen=hPgen,om=hom,
                rq=hrq,Vb=hVb,dJ=hdJ,Pacc=hPacc,lim=hlim,
                state=hstate,val=hval,Kkur=hKkur,
                N=N,m_q=m_q,label=label)

print("Simulant v9.4.1...")
res={}
for N,lbl in [(3,'3p'),(7,'7p')]:
    r=simula(N,lbl)
    Pa=np.mean(r['Pa']); Pg=np.mean(r['Pg'])
    mill=(Pg-Pa)/Pa*100
    Pb=np.mean(np.abs(r['Pb']))
    Pg_=np.mean(r['Pgen'])
    E_net=np.sum(r['Pb'])*DT/3600
    # distribucio estats
    cnt={STABLE:0,VALLEY:0,RECOVERY:0}
    for s in r['state']: cnt[int(s)]+=1
    print(f"\n  {lbl} ({N} pales):")
    print(f"    Millora P_grid:  +{mill:.2f}%")
    print(f"    P_buf |mig|:     {Pb/1e3:.2f} kW")
    print(f"    P_gen:           {Pg_:.0f} W")
    print(f"    E_net:           {E_net:.1f} kWh")
    print(f"    Estats: STABLE={cnt[STABLE]/steps*100:.0f}%  VALLEY={cnt[VALLEY]/steps*100:.0f}%  RECOVERY={cnt[RECOVERY]/steps*100:.0f}%")
    print(f"    K_kur mig:       {np.mean(r['Kkur']):.4f}")
    print(f"    omega*R_hub min: {np.min(r['lim']):.2f} m/s >> DR_MAX={DR_MAX}")
    res[lbl]=dict(r=r,mill=mill,Pb=Pb,Pgen=Pg_,E_net=E_net)

# ============================================================
# GRAFICS
# ============================================================
BG='#0d0d1a'; PAN='#13132b'
C3='#00d2ff'; C7='#00ff88'; CV='#ffd700'; CX='#e74c3c'
STATE_COLS={STABLE:'#185FA5',VALLEY:'#e74c3c',RECOVERY:'#00ff88'}
STATE_NOMS={STABLE:'STABLE',VALLEY:'VALLEY',RECOVERY:'RECOVERY'}

fig=plt.figure(figsize=(24,20),facecolor=BG)
gs=gridspec.GridSpec(4,3,figure=fig,hspace=0.50,wspace=0.38)

def sty(ax,t,xl='t[s]',yl=''):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9.5,pad=5,fontweight='bold')
    ax.tick_params(colors='#aaa',labelsize=8)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8.5); ax.set_ylabel(yl,color='#aaa',fontsize=8.5)
    ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)

def sh(ax): ax.axvspan(35,50,alpha=0.08,color=CX)

# G1: Vent
ax=fig.add_subplot(gs[0,0])
ax.plot(t_vec,v_arr,color=CV,lw=1.8)
ax.fill_between(t_vec,0,v_arr,alpha=0.1,color=CV)
ax.axhline(vr,color='white',ls='--',lw=1,alpha=0.5,label=f'vr={vr:.1f}m/s')
sh(ax); ax.legend(fontsize=8,framealpha=0.3)
sty(ax,'G1 — Vent realista (OU + calmant)','t[s]','v[m/s]')

# G2: P_grid
ax=fig.add_subplot(gs[0,1:])
r3=res['3p']['r']; r7=res['7p']['r']
ax.plot(t_vec,r3['Pa']/1e6,'--',color='#888',lw=1.2,alpha=0.6,label='Base aero')
ax.plot(t_vec,r3['Pg']/1e6,color=C3,lw=2,label=f"3p +{res['3p']['mill']:.2f}%  Pgen={res['3p']['Pgen']:.0f}W")
ax.plot(t_vec,r7['Pg']/1e6,color=C7,lw=2,label=f"7p +{res['7p']['mill']:.2f}%  Pgen={res['7p']['Pgen']:.0f}W")
sh(ax); ax.legend(fontsize=8,framealpha=0.3)
sty(ax,'G2 — P_grid (MW) v9.4.1 — 3p vs 7p','t[s]','MW')

# G3: Maquina d'estats
ax=fig.add_subplot(gs[1,0])
state_arr=r7['state'].astype(int)
for st,col in STATE_COLS.items():
    mask=state_arr==st
    ax.fill_between(t_vec,st-0.4,st+0.4,where=mask,color=col,alpha=0.7,label=STATE_NOMS[st])
ax.set_yticks([STABLE,VALLEY,RECOVERY])
ax.set_yticklabels(['STABLE','VALLEY','RECOVERY'],fontsize=8)
sh(ax); ax.legend(fontsize=7.5,framealpha=0.3,loc='upper right')
sty(ax,'G3 — Maquina d\'estats (7p)','t[s]','Estat')

# G4: Valley filtrada + K_kur adaptatiu
ax=fig.add_subplot(gs[1,1])
ax2=ax.twinx()
ax.plot(t_vec,r7['val'],color='#D85A30',lw=1.5,label='valley (filtr.)')
ax.axhline(THRESH,color=CX,ls=':',lw=0.8,alpha=0.7,label=f'+{THRESH}')
ax.axhline(-THRESH,color=C7,ls=':',lw=0.8,alpha=0.7,label=f'-{THRESH}')
ax2.plot(t_vec,r7['Kkur'],color='#9b59b6',lw=1.5,label='K_kur adapt.')
ax.set_ylabel('valley',color='#D85A30',fontsize=8); ax.tick_params(colors='#D85A30')
ax2.set_ylabel('K_kur',color='#9b59b6',fontsize=8); ax2.tick_params(colors='#9b59b6')
ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)
ax.set_xlabel('t[s]',color='#aaa',fontsize=8.5)
ax.set_title('G4 — Valley filtrada + K_kur adaptatiu (7p)',color='white',fontsize=9.5,pad=5,fontweight='bold')
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2,l1+l2,fontsize=7,framealpha=0.3)
sh(ax)

# G5: r_q per pala
ax=fig.add_subplot(gs[1,2])
COLS_P=[C3,C7,CV,'#ff6b35','#D85A30','#9b59b6','#e91e63']
for i in range(min(3,r7['N'])):
    ax.plot(t_vec,r7['rq'][:,i],color=COLS_P[i],lw=1.5,alpha=0.85,label=f'pala {i+1}')
ax.axhline(R_hub,color='#555',ls=':',lw=0.8,alpha=0.5)
ax.axhline(R_tip,color='#555',ls=':',lw=0.8,alpha=0.5)
ax.legend(fontsize=7.5,framealpha=0.3); sh(ax)
sty(ax,'G5 — r_q(t) per pala (7p)','t[s]','m')

# G6: dJ/dt i P_buf
ax=fig.add_subplot(gs[2,0])
ax2=ax.twinx()
ax.plot(t_vec,r3['dJ'],color='#D85A30',lw=1.2,alpha=0.8,label='dJ/dt 3p')
ax.plot(t_vec,r7['dJ'],color='#BA7517',lw=1.2,alpha=0.7,label='dJ/dt 7p')
ax.fill_between(t_vec,0,r7['dJ'],where=np.array(r7['dJ'])>0,alpha=0.1,color=CX)
ax.fill_between(t_vec,0,r7['dJ'],where=np.array(r7['dJ'])<0,alpha=0.1,color=C7)
ax2.plot(t_vec,r7['Pb']/1e3,color=C7,lw=1.5,alpha=0.9,label='P_buf 7p (kW)')
ax.set_ylabel('dJ/dt [kg m2/s]',color='#D85A30',fontsize=8); ax.tick_params(colors='#D85A30')
ax2.set_ylabel('P_buf [kW]',color=C7,fontsize=8); ax2.tick_params(colors=C7)
ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)
ax.set_xlabel('t[s]',color='#aaa',fontsize=8.5)
ax.set_title('G6 — dJ/dt i P_buf Quijote',color='white',fontsize=9.5,pad=5,fontweight='bold')
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2,l1+l2,fontsize=7,framealpha=0.3)
sh(ax)

# G7: P_gen
ax=fig.add_subplot(gs[2,1])
ax.plot(t_vec,r3['Pgen'],color=C3,lw=1.5,label=f"3p P_gen mig={res['3p']['Pgen']:.0f}W")
ax.plot(t_vec,r7['Pgen'],color=C7,lw=1.5,label=f"7p P_gen mig={res['7p']['Pgen']:.0f}W")
ax.fill_between(t_vec,0,r7['Pgen'],alpha=0.15,color=C7)
ax.axhline(res['3p']['Pgen'],color=C3,ls='--',lw=0.8,alpha=0.6)
ax.axhline(res['7p']['Pgen'],color=C7,ls='--',lw=0.8,alpha=0.6)
ax.legend(fontsize=8,framealpha=0.3); sh(ax)
sty(ax,'G7 — P_gen hidraulic recuperat (eta=85%)','t[s]','W')

# G8: omega
ax=fig.add_subplot(gs[2,2])
ax.plot(t_vec,r3['om'],color=C3,lw=1.8,label='3p omega')
ax.plot(t_vec,r7['om'],color=C7,lw=1.8,label='7p omega')
ax.axhline(omR,color='white',ls='--',lw=1,alpha=0.4,label=f'omR={omR:.3f}')
ax.axhline(omR*0.85,color=CX,ls=':',lw=0.8,alpha=0.6,label='llindar VALLEY')
ax.legend(fontsize=7.5,framealpha=0.3); sh(ax)
sty(ax,'G8 — omega(t) — estabilitzacio rotor','t[s]','rad/s')

# G9: Taula comparativa v9.4 vs v9.4.1
ax=fig.add_subplot(gs[3,:])
ax.axis('off'); ax.set_facecolor('#0a0a14')
rows=[
    ['Parametre','3p v9.4','3p v9.4.1','7p v9.4','7p v9.4.1','Millora v9.4→v9.4.1'],
    ['Millora P_grid','+1.4%',f"+{res['3p']['mill']:.2f}%",
                       '+1.5%',f"+{res['7p']['mill']:.2f}%",'K_q + estats'],
    ['P_buf |mig|','12 kW',f"{res['3p']['Pb']/1e3:.1f}kW",
                   '16 kW',f"{res['7p']['Pb']/1e3:.1f}kW",'K_q=1e5'],
    ['P_gen','128 W',f"{res['3p']['Pgen']:.0f}W",
             '226 W',f"{res['7p']['Pgen']:.0f}W",'acumulador'],
    ['E_net','-30kWh',f"{res['3p']['E_net']:.0f}kWh",
              '-66kWh',f"{res['7p']['E_net']:.0f}kWh",'K_comp+estats'],
    ['K_kur','0.10 fix','adaptatiu','0.10 fix','adaptatiu','K0+K1*dom+K2*val'],
    ['omega*R_hub','6.6m/s >> 0.5','idem','idem','idem','sempre OK'],
]
tbl=ax.table(cellText=rows[1:],colLabels=rows[0],loc='center',cellLoc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(9)
for (rr,c),cell in tbl.get_celld().items():
    bg='#1F5C9E' if rr==0 else ('#0d1a0d' if rr%2==1 else '#0d0d1a')
    cell.set_facecolor(bg)
    cell.set_text_props(color='white')
    cell.set_edgecolor('#333355')

fig.suptitle(
    'QUIJOTE + ZYPYZAPE v9.4.1 — Model consolidat per al paper\n'
    'Victor Manzanares Alberola — EPSA UPV Alcoi — 2026\n'
    f'K_q=1e5 | K_comp={K_comp} | KQ={KQ} | K_kur adaptatiu | Valley filtrada | '
    'Maquina d\'estats STABLE/VALLEY/RECOVERY',
    color='white',fontsize=10,fontweight='bold',y=0.999)

out='/mnt/user-data/outputs/gemelo_v941_complet.png'
plt.savefig(out,dpi=155,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\nGrafic: {out}")
import shutil
shutil.copy('/home/claude/gemelo_v941.py','/mnt/user-data/outputs/gemelo_v941.py')
