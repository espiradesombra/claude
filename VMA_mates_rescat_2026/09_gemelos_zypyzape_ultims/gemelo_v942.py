"""
QUIJOTE + ZYPYZAPE v9.4.2
Victor Manzanares Alberola — EPSA UPV Alcoi

MILLORES vs v9.4.1 (LeeChat):
  1. K_q = 8e4 (de 1e5)
  2. THRESHOLD = 0.2 (de 0.1)
  3. TAU_VAL = 0.2 (de 0.1)
  4. K_KUR_MAX = 0.20
  5. Estat TRANSITION
  6. KQ = 0.5 (de 0.6)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

np.random.seed(42)

R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.482; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3); omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0

K_q=8e4; K_comp=0.1; c_fric=80.0; DR_MAX=0.5
KPITCH=2.0; KOMEGA=0.5; KQ=0.5
K0_kur=0.10; K1_kur=0.05; K2_kur=0.01; K_KUR_MAX=0.20
TAU_VAL=0.2; THRESH=0.2
STABLE=0; VALLEY=1; RECOVERY=2; TRANSITION=3
RHO_FL=3386.0; D_CANAL=0.05; A_can=np.pi*(D_CANAL/2)**2
V_BUF_MAX=0.10; Q_BOMBA=0.010; P_MAX_ACC=20e6; LLINDAR=20.0; ETA_GEN=0.85
V_BASE=11.4; T_TOT=80.0; DT=0.05
steps=int(T_TOT/DT); t_vec=np.arange(steps)*DT

print(f"v9.4.2: K_q={K_q:.0e} THRESH={THRESH} TAU={TAU_VAL} K_KUR_MAX={K_KUR_MAX}")

ou=np.zeros(steps)
for i in range(1,steps):
    ou[i]=ou[i-1]-0.5*ou[i-1]*DT+0.8*np.sqrt(DT)*np.random.randn()
v_arr=np.array([max(3,V_BASE+2*np.sin(2*np.pi*s*DT/20)+ou[s]*0.6
                    -(5*np.sin(np.pi*(s*DT-35)/15) if 35<s*DT<50 else 0))
                for s in range(steps)])

def eta_lam(om,v):
    if v<=0: return 0.0
    return float(max(0,1-((om*R/v-lam_opt)/lam_opt)**2))

def simula(N):
    m_q=RHO_FL*A_can*(R_tip-R_hub)
    thetas=np.array([2*np.pi*i/N for i in range(N)])
    omega=omR*0.95; r_q=np.full(N,(R_hub+R_tip)/2); dr_q=np.zeros(N)
    V_buf=0.0; P_acc=0.0; valley_f=0.0; state=STABLE

    hPa=np.zeros(steps); hPg=np.zeros(steps); hPgen=np.zeros(steps)
    hPb=np.zeros(steps); hom=np.zeros(steps); hstate=np.zeros(steps)
    hKkur=np.zeros(steps); hval=np.zeros(steps); hdJ=np.zeros(steps)

    for s in range(steps):
        v=v_arr[s]
        omega+=(min(omR,lam_opt*max(v,0.1)/R)-omega)*(DT/3.0)
        omega=float(np.clip(omega,omR*0.3,omR*1.4))
        dJ=0.0
        for i in range(N):
            th=thetas[i]%(2*np.pi)
            r_opt0=R_tip-(R_tip-R_hub)*(1+np.cos(th))/2
            F_c=m_q*omega**2*r_q[i]
            r_opt=r_opt0+K_comp*(F_c/(K_q*(R_tip-R_hub)+1e-9))
            r_opt=float(np.clip(r_opt,R_hub,R_tip))
            if state==VALLEY: r_opt=R_tip
            elif state==RECOVERY: r_opt=R_hub
            elif state==TRANSITION: r_opt=r_opt0
            F_ct=-K_q*(r_q[i]-r_opt); F_fr=-c_fric*dr_q[i]
            dr_q[i]+=(F_c+F_ct+F_fr)/m_q*DT
            dr_q[i]=float(np.clip(dr_q[i],-DR_MAX,DR_MAX))
            dr_q[i]=float(np.clip(dr_q[i],-omega*r_q[i],omega*r_q[i]))
            r_q[i]=float(np.clip(r_q[i]+dr_q[i]*DT,R_hub,R_tip))
            dJ+=2*m_q*r_q[i]*dr_q[i]

        P_buf=-0.5*dJ*omega**2
        valley_raw=-omega*dJ
        alpha_v=DT/(TAU_VAL+DT)
        valley_f=alpha_v*valley_raw+(1-alpha_v)*valley_f

        if valley_f>THRESH and omega<omR*0.90: state=VALLEY
        elif valley_f<-THRESH and omega>omR*0.95: state=RECOVERY
        elif abs(valley_f)<=THRESH*0.5 and 0.95*omR<=omega<=1.05*omR: state=STABLE
        else: state=TRANSITION

        K_kur=float(np.clip(K0_kur+K1_kur*max(0,(omR-omega)/omR)+K2_kur*abs(valley_f),0,K_KUR_MAX))

        for i in range(N):
            dJi=2*m_q*r_q[i]*dr_q[i]
            pm=float(np.clip(KPITCH*np.cos(thetas[i])+KOMEGA*omega*np.sin(thetas[i])-KQ*dJi,-8,8))
            kc=(K_kur/N)*float(np.sum(np.sin(thetas-thetas[i])))
            thetas[i]+=(omega+kc+pm*0.01)*DT

        Pgen=0.0
        if P_buf>LLINDAR and V_buf<V_BUF_MAX:
            P_acc=min(P_acc+P_buf*DT/V_BUF_MAX,P_MAX_ACC)
            Q=min(Q_BOMBA,P_buf/max(P_acc,1e3))
            V_buf=min(V_BUF_MAX,V_buf+Q*DT)
        elif P_buf<-LLINDAR and V_buf>0:
            P_acc=max(P_acc+P_buf*DT/V_BUF_MAX,0)
            Q=min(Q_BOMBA,-P_buf/max(P_acc,1e3))
            Pgen=P_acc*Q*ETA_GEN
            V_buf=max(0,V_buf-Q*DT)

        etL=eta_lam(omega,v)
        Pa=0.5*RHO_A*A_rot*Cp_max*etL*v**3
        Pg=max(0,Pa*(1+0.04)+max(0,P_buf)+Pgen)
        hPa[s]=Pa; hPg[s]=Pg; hPgen[s]=Pgen; hPb[s]=P_buf
        hom[s]=omega; hstate[s]=state; hKkur[s]=K_kur; hval[s]=valley_f; hdJ[s]=dJ

    return dict(Pa=hPa,Pg=hPg,Pgen=hPgen,Pb=hPb,om=hom,
                state=hstate,Kkur=hKkur,val=hval,dJ=hdJ,m_q=m_q)

print("Simulant...")
res={'3p':simula(3),'7p':simula(7)}

STATE_N={STABLE:'STABLE',VALLEY:'VALLEY',RECOVERY:'RECOVERY',TRANSITION:'TRANSITION'}
print("\n"+"="*60)
for lbl,r in res.items():
    mill=(np.mean(r['Pg'])-np.mean(r['Pa']))/np.mean(r['Pa'])*100
    Pb=np.mean(np.abs(r['Pb']))
    Pgen=np.mean(r['Pgen'])
    E_net=np.sum(r['Pb'])*DT/3600
    cnt={s:0 for s in [STABLE,VALLEY,RECOVERY,TRANSITION]}
    for s in r['state']: cnt[int(s)]+=1
    print(f"\n{lbl}:")
    print(f"  Millora P_grid:  +{mill:.2f}%")
    print(f"  P_buf |mig|:     {Pb/1e3:.2f} kW")
    print(f"  P_gen:           {Pgen:.0f} W")
    print(f"  E_net:           {E_net:.1f} kWh")
    estats=" | ".join(f"{STATE_N[s]}={cnt[s]/steps*100:.0f}%" for s in [STABLE,VALLEY,RECOVERY,TRANSITION])
    print(f"  Estats: {estats}")
    print(f"  K_kur mig:       {np.mean(r['Kkur']):.4f}")
    res[lbl].update({'mill':mill,'Pb_abs':Pb,'Pgen_mig':Pgen,'E_net':E_net})

# ============================================================
# GRAFICS: comparativa v9.4 / v9.4.1 / v9.4.2
# ============================================================
BG='#0d0d1a'; PAN='#13132b'
C3='#00d2ff'; C7='#00ff88'; CV='#ffd700'; CX='#e74c3c'
STATE_COLS={STABLE:'#185FA5',VALLEY:'#e74c3c',RECOVERY:'#00ff88',TRANSITION:'#ffd700'}

fig=plt.figure(figsize=(22,18),facecolor=BG)
gs=gridspec.GridSpec(3,3,figure=fig,hspace=0.50,wspace=0.38)

def sty(ax,t,xl='t[s]',yl=''):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9.5,pad=5,fontweight='bold')
    ax.tick_params(colors='#aaa',labelsize=8)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8.5); ax.set_ylabel(yl,color='#aaa',fontsize=8.5)
    ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)
def sh(ax): ax.axvspan(35,50,alpha=0.08,color=CX)

r3=res['3p']; r7=res['7p']

ax=fig.add_subplot(gs[0,0])
ax.plot(t_vec,v_arr,color=CV,lw=1.8)
ax.fill_between(t_vec,0,v_arr,alpha=0.1,color=CV)
ax.axhline(vr,color='white',ls='--',lw=1,alpha=0.5,label=f'vr={vr:.1f}m/s')
sh(ax); ax.legend(fontsize=8,framealpha=0.3)
sty(ax,'G1 — Vent (OU + calmant)','t[s]','v[m/s]')

ax=fig.add_subplot(gs[0,1:])
ax.plot(t_vec,r3['Pa']/1e6,'--',color='#888',lw=1.2,alpha=0.6,label='Base')
ax.plot(t_vec,r3['Pg']/1e6,color=C3,lw=2,label=f"3p +{r3['mill']:.2f}%  Pgen={r3['Pgen_mig']:.0f}W")
ax.plot(t_vec,r7['Pg']/1e6,color=C7,lw=2,label=f"7p +{r7['mill']:.2f}%  Pgen={r7['Pgen_mig']:.0f}W")
sh(ax); ax.legend(fontsize=8,framealpha=0.3)
sty(ax,'G2 — P_grid (MW) v9.4.2','t[s]','MW')

ax=fig.add_subplot(gs[1,0])
state_arr7=r7['state'].astype(int)
for st,col in STATE_COLS.items():
    ax.fill_between(t_vec,st-0.4,st+0.4,where=state_arr7==st,color=col,alpha=0.75)
ax.set_yticks([0,1,2,3]); ax.set_yticklabels(['STABLE','VALLEY','RECOVERY','TRANSITION'],fontsize=7.5)
sh(ax)
sty(ax,'G3 — Maquina estats (7p)','t[s]','Estat')

ax=fig.add_subplot(gs[1,1])
ax2=ax.twinx()
ax.plot(t_vec,r7['val'],color='#D85A30',lw=1.5,label='valley')
ax.axhline(THRESH,color=CX,ls=':',lw=0.8,alpha=0.7)
ax.axhline(-THRESH,color=C7,ls=':',lw=0.8,alpha=0.7)
ax2.plot(t_vec,r7['Kkur'],color='#9b59b6',lw=1.5,label='K_kur')
ax2.axhline(K_KUR_MAX,color='#ffd700',ls='--',lw=0.8,alpha=0.6,label=f'K_KUR_MAX={K_KUR_MAX}')
ax.set_ylabel('valley',color='#D85A30',fontsize=8); ax.tick_params(colors='#D85A30')
ax2.set_ylabel('K_kur',color='#9b59b6',fontsize=8); ax2.tick_params(colors='#9b59b6')
ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)
ax.set_xlabel('t[s]',color='#aaa',fontsize=8.5)
ax.set_title('G4 — Valley + K_kur adapt. (7p)',color='white',fontsize=9.5,pad=5,fontweight='bold')
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2,l1+l2,fontsize=7,framealpha=0.3)
sh(ax)

ax=fig.add_subplot(gs[1,2])
ax.plot(t_vec,r3['Pgen'],color=C3,lw=1.5,label=f"3p mig={r3['Pgen_mig']:.0f}W")
ax.plot(t_vec,r7['Pgen'],color=C7,lw=1.5,label=f"7p mig={r7['Pgen_mig']:.0f}W")
ax.fill_between(t_vec,0,r7['Pgen'],alpha=0.15,color=C7)
ax.legend(fontsize=8,framealpha=0.3); sh(ax)
sty(ax,'G5 — P_gen hidraulic (W)','t[s]','W')

ax=fig.add_subplot(gs[2,0])
ax2=ax.twinx()
ax.plot(t_vec,r7['dJ'],color='#D85A30',lw=1.2,alpha=0.8,label='dJ/dt')
ax2.plot(t_vec,r7['Pb']/1e3,color=C7,lw=1.5,label='P_buf (kW)')
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.set_ylabel('dJ/dt',color='#D85A30',fontsize=8); ax.tick_params(colors='#D85A30')
ax2.set_ylabel('P_buf [kW]',color=C7,fontsize=8); ax2.tick_params(colors=C7)
ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)
ax.set_xlabel('t[s]',color='#aaa',fontsize=8.5)
ax.set_title('G6 — dJ/dt + P_buf (7p)',color='white',fontsize=9.5,pad=5,fontweight='bold')
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2,l1+l2,fontsize=7,framealpha=0.3)
sh(ax)

ax=fig.add_subplot(gs[2,1])
ax.plot(t_vec,r3['om'],color=C3,lw=1.8,label='3p')
ax.plot(t_vec,r7['om'],color=C7,lw=1.8,label='7p')
ax.axhline(omR,color='white',ls='--',lw=1,alpha=0.4,label=f'omR={omR:.3f}')
ax.axhline(omR*0.90,color=CX,ls=':',lw=0.8,alpha=0.6,label='llindar VALLEY (90%)')
ax.legend(fontsize=7.5,framealpha=0.3); sh(ax)
sty(ax,'G7 — omega(t)','t[s]','rad/s')

ax=fig.add_subplot(gs[2,2])
ax.axis('off'); ax.set_facecolor('#0a0a14')
rows=[
    ['','v9.4','v9.4.1','v9.4.2'],
    ['K_q','5e4','1e5','8e4'],
    ['THRESH','—','0.1','0.2'],
    ['TAU_VAL','0.1','0.1','0.2'],
    ['K_KUR_MAX','—','0.30','0.20'],
    ['KQ','0.8','0.6','0.5'],
    ['Millora 3p','+1.4%','+8.4%',f"+{r3['mill']:.1f}%"],
    ['Millora 7p','+1.5%','+9.5%',f"+{r7['mill']:.1f}%"],
    ['P_gen 3p','128W','17W',f"{r3['Pgen_mig']:.0f}W"],
    ['P_gen 7p','226W','40W',f"{r7['Pgen_mig']:.0f}W"],
    ['E_net 3p','-30kWh','-44kWh',f"{r3['E_net']:.0f}kWh"],
    ['STABLE 7p','0%','0%','?%'],
]
tbl=ax.table(cellText=rows[1:],colLabels=rows[0],loc='center',cellLoc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5)
for (rr,c),cell in tbl.get_celld().items():
    bg='#1F5C9E' if rr==0 else ('#0d1a0d' if c==3 and rr>0 else '#0d0d1a')
    cell.set_facecolor(bg)
    cell.set_text_props(color='white')
    cell.set_edgecolor('#333355')
ax.set_title('G8 — Comparativa v9.4 / v9.4.1 / v9.4.2',color='white',fontsize=9.5,pad=5,fontweight='bold')

fig.suptitle(
    f'QUIJOTE v9.4.2 — K_q=8e4 | THRESH=0.2 | TAU=0.2 | K_KUR_MAX=0.20 | TRANSITION\n'
    f'Victor Manzanares Alberola — EPSA UPV Alcoi — 2026',
    color='white',fontsize=10,fontweight='bold',y=0.999)

out='/mnt/user-data/outputs/gemelo_v942_complet.png'
plt.savefig(out,dpi=155,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\nGrafic: {out}")
import shutil; shutil.copy('/home/claude/gemelo_v942.py','/mnt/user-data/outputs/gemelo_v942.py')
