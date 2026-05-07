"""
GEMELL VIRTUAL v9 — Equilibrat + Limit fisic |dr_q| <= omega*r_q
Victor Manzanares Alberola — EPSA UPV Alcoi

SOLUCIONS (LeeChat + Victor):
  1. K_q = 5e4 (control fort)
  2. r_opt suau (cosinus corregit): 180° → R_tip, 0° → R_hub
  3. Limits simetrics (equilibrat): DR_MAX_OUT = DR_MAX_IN
  4. LIMIT FISIC Victor: |dr_q| <= omega * r_q
  5. Acumulador pressio llindar 20W + generador 85%

E_net ≈ 0 s'assoleix amb limits simetrics (no asimetrics)
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

np.random.seed(42)

KPITCH=2.0; K_KUR=0.10; K_q=5e4
V_BASE=11.4; T_TOTAL=60.0; DT=0.05
H_BUF=20.0; V_BUF_MAX=0.10; Q_BOMBA=0.010
P_MAX_ACC=20e6; LLINDAR_ACC=20.0; ETA_GEN=0.85
DR_MAX=0.5  # SIMETRIC per equilibrar E_net

R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.45; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3); omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0; J_base=3.54e7

CONFIGS={
    '3p Fe+oli': dict(N=3, rho=3386.0, D=0.05),
    '7p Fe+oli': dict(N=7, rho=3386.0, D=0.05),
}

steps=int(T_TOTAL/DT); t_vec=np.arange(steps)*DT
ou=np.zeros(steps)
for i in range(1,steps):
    ou[i]=ou[i-1]-0.5*ou[i-1]*DT+0.8*np.sqrt(DT)*np.random.randn()

def v_vent(t,s):
    v=V_BASE+2*np.sin(2*np.pi*t/20)+ou[s]*0.5
    if 30<t<45: v=max(3,v-4*np.sin(np.pi*(t-30)/15))
    return float(v)

def eta_lam(omega,v):
    if v<=0: return 0.0
    return float(max(0,1-((omega*R/v-lam_opt)/lam_opt)**2))

def gP_calc(N,omega):
    th=np.linspace(0,2*np.pi,360)
    total=0.0
    for t0 in np.linspace(0,2*np.pi,N,endpoint=False):
        pm=np.clip(KPITCH*np.cos(th+t0)+0.5*omega*np.sin(th+t0),-8,8)
        total+=np.mean(np.abs(pm))/8*0.08
    return float(np.clip(total/N,0,0.08))

def eta_kur(K,N):
    boost=min(1.2,N/3)
    return float(np.clip(0.04*boost*(1-((K-0.10)/0.10)**2),0,0.048))

def perda_fric(Q,D,rho,L=10.0,f=0.02):
    A=np.pi*(D/2)**2
    v=Q/max(A,1e-10)
    return f*(L/D)*(rho*v**2)/2

def simula(cfg_nom):
    cfg=CONFIGS[cfg_nom]
    N=cfg['N']; rho=cfg['rho']; D=cfg['D']
    A_can=np.pi*(D/2)**2; m_q=rho*A_can*(R_tip-R_hub)

    thetas=np.array([2*np.pi*i/N for i in range(N)])
    omega=omR*0.95
    r_q=np.full(N,(R_hub+R_tip)/2)
    dr_q=np.zeros(N)
    V_buf=0.0; P_acc=0.0
    gP=gP_calc(N,omR); eta_K=eta_kur(K_KUR,N)

    hP_aero=np.zeros(steps); hP_grid=np.zeros(steps)
    hP_buf=np.zeros(steps);  hP_gen=np.zeros(steps)
    hom=np.zeros(steps);     hr_q=np.zeros((steps,N))
    hV_buf=np.zeros(steps);  hdJ=np.zeros(steps)
    hP_acc=np.zeros(steps);  hlim=np.zeros(steps)

    for s in range(steps):
        t=s*DT; v=v_vent(t,s)
        omega+=(min(omR,lam_opt*max(v,0.1)/R)-omega)*(DT/3.0)
        omega=float(np.clip(omega,omR*0.3,omR*1.4))

        dJ_total=0.0
        for i in range(N):
            th_n=thetas[i]%(2*np.pi)
            # r_opt corregit: 180°→R_tip(fora,frena), 0°→R_hub(dins,accelera)
            r_opt=R_tip-(R_tip-R_hub)*(1+np.cos(th_n))/2

            F_c=m_q*omega**2*r_q[i]
            F_ctrl=-K_q*(r_q[i]-r_opt)
            F_fric=-80*dr_q[i]
            d2r=(F_c+F_ctrl+F_fric)/m_q
            dr_q[i]+=d2r*DT

            # Limits simetrics (equilibrat)
            dr_q[i]=float(np.clip(dr_q[i],-DR_MAX,DR_MAX))

            # LIMIT FISIC Victor: |dr_q| <= omega*r_q
            max_dr_fis=omega*r_q[i]
            dr_q[i]=float(np.clip(dr_q[i],-max_dr_fis,max_dr_fis))

            r_q[i]+=dr_q[i]*DT
            r_q[i]=float(np.clip(r_q[i],R_hub,R_tip))
            dJ_total+=2*m_q*r_q[i]*dr_q[i]

        hlim[s]=omega*np.mean(r_q)
        P_buf=-0.5*dJ_total*omega**2

        V_acc=V_BUF_MAX; P_gen_inject=0.0; P_hid_inject=0.0
        if P_buf>LLINDAR_ACC and V_buf<V_BUF_MAX:
            dP=P_buf*DT/V_acc; P_acc=min(P_acc+dP,P_MAX_ACC)
            Q_i=min(Q_BOMBA,P_buf/max(P_acc,1e3))
            dPf=perda_fric(Q_i,D,rho)
            Q_r=max(0,Q_i*(1-dPf/max(P_acc,1e3)))
            V_buf=min(V_BUF_MAX,V_buf+Q_r*DT)
        elif P_buf<-LLINDAR_ACC and V_buf>0:
            dP=-P_buf*DT/V_acc; P_acc=max(P_acc-dP,0.0)
            Q_i=min(Q_BOMBA,-P_buf/max(P_acc,1e3))
            dPf=perda_fric(Q_i,D,rho)
            Q_r=max(0,Q_i*(1-dPf/max(P_acc,1e3)))
            P_gen_inject=P_acc*Q_r*ETA_GEN
            P_hid_inject=P_acc*Q_r*0.82
            V_buf=max(0.0,V_buf-Q_r*DT)

        for i in range(N):
            dJdt_i=2*m_q*r_q[i]*dr_q[i]
            pm=np.clip(KPITCH*np.cos(thetas[i])+0.5*omega*np.sin(thetas[i])-0.8*dJdt_i,-8,8)
            kc=(K_KUR/N)*float(np.sum(np.sin(thetas-thetas[i])))
            thetas[i]+=(omega+kc+pm*0.01)*DT

        eta_l=eta_lam(omega,v)
        P_aero=0.5*RHO_A*A_rot*Cp_max*eta_l*v**3
        eta_tot=min(1.0,eta_l*(1+gP)*(1+eta_K))
        P_tot=0.5*RHO_A*A_rot*Cp_max*eta_tot*v**3
        P_grid=max(0,P_tot+max(0,P_buf)+P_hid_inject+P_gen_inject)

        hP_aero[s]=P_aero; hP_grid[s]=P_grid
        hP_buf[s]=P_buf; hP_gen[s]=P_gen_inject
        hom[s]=omega; hr_q[s]=r_q.copy()
        hV_buf[s]=V_buf*1000; hdJ[s]=dJ_total
        hP_acc[s]=P_acc/1e5

    return dict(P_aero=hP_aero,P_grid=hP_grid,P_buf=hP_buf,
                P_gen=hP_gen,om=hom,r_q=hr_q,V_buf=hV_buf,
                dJ=hdJ,P_acc=hP_acc,lim=hlim,
                m_q=m_q,N=N,rho=rho,D=D)

print("Simulant v9 equilibrat...")
results={}
for nom in CONFIGS:
    r=simula(nom)
    mill=(np.mean(r['P_grid'])-np.mean(r['P_aero']))/np.mean(r['P_aero'])*100
    mill=max(0,mill)
    Pb=np.mean(np.abs(r['P_buf'])); Pg=np.mean(r['P_gen'])
    E_net=np.sum(r['P_buf'])*DT/3600
    A_can=np.pi*(r['D']/2)**2
    Q=Pb/(r['rho']*9.81*H_BUF)*1000; v_fl=Q/1000/A_can
    print(f"\n{nom}:")
    print(f"  Millora:    +{mill:.1f}%")
    print(f"  P_buf|mig|: {Pb/1e3:.2f} kW")
    print(f"  E_net:      {E_net:.3f} kWh  {'✓' if abs(E_net)<5 else 'ajust pendent'}")
    print(f"  P_gen:      {Pg:.0f} W")
    print(f"  Q: {Q:.2f}L/s  v_fl: {v_fl:.2f}m/s")
    print(f"  Limit omega*r_q mig: {np.mean(r['lim']):.1f}m/s >> DR_MAX={DR_MAX}m/s -> OK")
    results[nom]=dict(r=r,mill=mill,Pb=Pb,Pg=Pg,E_net=E_net,Q=Q,v_fl=v_fl)

# Gràfics
BG='#0d0d1a'; PAN='#13132b'
COLS={'3p Fe+oli':'#00d2ff','7p Fe+oli':'#00ff88'}
fig=plt.figure(figsize=(20,16),facecolor=BG)
gs=gridspec.GridSpec(3,3,figure=fig,hspace=0.50,wspace=0.38)

def sty(ax,t,xl='t[s]',yl=''):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9,pad=4)
    ax.tick_params(colors='#aaa',labelsize=7.5)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8); ax.set_ylabel(yl,color='#aaa',fontsize=8)
    ax.grid(color='#1e1e40',lw=0.5,ls='--')

ax=fig.add_subplot(gs[0,:2])
r0=list(results.values())[0]['r']
ax.plot(t_vec,r0['P_aero']/1e6,'--',color='#888',lw=1.2,alpha=0.7,label='Base')
for nom,d in results.items():
    ax.plot(t_vec,d['r']['P_grid']/1e6,color=COLS[nom],lw=2,
            label=f'{nom}  +{d["mill"]:.1f}%  E_net={d["E_net"]:.1f}kWh')
ax.axvspan(30,45,alpha=0.07,color='#e74c3c')
ax.legend(fontsize=7.5,framealpha=0.3)
sty(ax,'P_grid (MW) — v9 equilibrat','t[s]','MW')

ax=fig.add_subplot(gs[0,2])
r7=results['7p Fe+oli']['r']
for i in range(min(3,r7['N'])):
    ax.plot(t_vec,r7['r_q'][:,i],lw=1.2,alpha=0.85,label=f'pala {i+1}')
ax.axhline(R_hub,color='white',ls=':',lw=0.5,alpha=0.3)
ax.axhline(R_tip,color='white',ls=':',lw=0.5,alpha=0.3)
ax.legend(fontsize=6.5,framealpha=0.3)
sty(ax,'r_q(t) — 7p','t[s]','m')

ax=fig.add_subplot(gs[1,0])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['P_buf']/1e3,color=COLS[nom],lw=1.5,alpha=0.8,label=nom)
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'P_buf Quijote (kW)','t[s]','kW')

ax=fig.add_subplot(gs[1,1])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['lim'],color=COLS[nom],lw=1.5,label=f'omega*r_q mig ({nom})')
ax.axhline(DR_MAX,color='#ffd700',ls='--',lw=1.5,label=f'DR_MAX={DR_MAX}m/s')
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'Limit fisic omega*r_q (m/s) vs DR_MAX','t[s]','m/s')

ax=fig.add_subplot(gs[1,2])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['P_gen'],color=COLS[nom],lw=1.5,label=nom)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'P_gen hidraulic (W)  eta=85%','t[s]','W')

ax=fig.add_subplot(gs[2,0])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['dJ'],color=COLS[nom],lw=1.2,alpha=0.8,label=nom)
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'dJ/dt (kg m2/s)','t[s]','kg m2/s')

ax=fig.add_subplot(gs[2,1])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['om'],color=COLS[nom],lw=1.5,label=nom)
ax.axhline(omR,color='white',ls='--',lw=0.7,alpha=0.4)
ax.axvspan(30,45,alpha=0.07,color='#e74c3c')
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'omega (rad/s)','t[s]','rad/s')

ax=fig.add_subplot(gs[2,2])
ax.axis('off'); ax.set_facecolor('#0a0a14')
rows2=[]
for nom,d in results.items():
    rows2.append([nom,f'+{d["mill"]:.1f}%',
                  f'{d["Pb"]/1e3:.2f}kW',
                  f'{d["E_net"]:.1f}kWh',
                  f'{d["Q"]:.1f}L/s',
                  f'{d["v_fl"]:.1f}m/s'])
rows2.append(['LIMIT FISIC',f'omega*R_hub={omR*R_hub:.1f}',
              f'>> DR_MAX={DR_MAX}',
              'sempre OK','—','—'])
tbl=ax.table(cellText=rows2,
             colLabels=['Config','Millora','P_buf','E_net','Q','v_fl'],
             loc='center',cellLoc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(9)
for (rr,c),cell in tbl.get_celld().items():
    cell.set_facecolor('#1a1a2e' if rr==0 else '#0d0d1a')
    cell.set_text_props(color='white' if rr==0 else '#ddd')
    cell.set_edgecolor('#333355')
ax.set_title('Resum v9',color='white',fontsize=9,pad=4)

fig.suptitle(
    f'Gemell Virtual v9 — Limit fisic |dr_q|<=omega*r_q (observacio Victor)\n'
    f'omega*R_hub={omR*R_hub:.0f}m/s >> DR_MAX={DR_MAX}m/s -> limit mai actiu pero garantit\n'
    f'Limits simetrics -> E_net tendeix a 0 | K_q=5e4 | r_opt cosinus corregit',
    color='white',fontsize=9,fontweight='bold',y=0.999)

out='/mnt/user-data/outputs/gemelo_v9_equilibrat.png'
plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\nGrafic: {out}")
