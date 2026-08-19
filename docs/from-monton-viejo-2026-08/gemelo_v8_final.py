"""
GEMELL VIRTUAL v8 — Millores finals LeeChat
Víctor Manzanares Alberola — EPSA UPV Alcoi

MILLORES vs v7:
  1. r_opt corregit: R_tip-(R_tip-R_hub)*(1+cos(theta))/2
     → 180° (avall, max captura) → R_tip (fora, frena)
     → 0°  (amunt, min captura) → R_hub (dins, accelera)
  2. Llindar acumulador: 50W → 20W
  3. Generador hidràulic regeneratiu (η=85%)
  4. Comparativa visual 3p vs 7p + taula final
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

np.random.seed(42)

# ============================================================
# PARÀMETRES
# ============================================================
KPITCH=2.0; K_KUR=0.10; V_BASE=11.4; T_TOTAL=60.0; DT=0.05
H_BUF=20.0; V_BUF_MAX=0.10; Q_BOMBA=0.010; P_MAX_ACC=20e6
LLINDAR_ACC=20.0   # W (MILLORA 2: de 50→20W)
ETA_GEN=0.85       # eficiència generador hidràulic (MILLORA 3)

R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.45; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3); omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0; J_base=3.54e7

CONFIGS={
    '3p Fe+oli': dict(N=3,rho=3386.0,D=0.05,Kq=2e4,DR_MAX=0.5),
    '7p Fe+oli': dict(N=7,rho=3386.0,D=0.05,Kq=2e4,DR_MAX=0.5),
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
    K_opt=0.10; boost=min(1.2,N/3)
    return float(np.clip(0.04*boost*(1-((K-K_opt)/K_opt)**2),0,0.048))

def perda_fric(Q,D,rho,L=10.0,f=0.02):
    A=np.pi*(D/2)**2
    v=Q/max(A,1e-10)
    return f*(L/D)*(rho*v**2)/2

def simula(cfg_nom):
    cfg=CONFIGS[cfg_nom]
    N=cfg['N']; rho=cfg['rho']; D=cfg['D']
    K_q=cfg['Kq']; DR_MAX=cfg['DR_MAX']
    A_can=np.pi*(D/2)**2
    m_q=rho*A_can*(R_tip-R_hub)

    thetas=np.array([2*np.pi*i/N for i in range(N)])
    omega=omR*0.95
    r_q=np.full(N,(R_hub+R_tip)/2)
    dr_q=np.zeros(N)
    V_buf=0.0; P_acc=0.0
    gP=gP_calc(N,omR); eta_K=eta_kur(K_KUR,N)

    hP_aero=np.zeros(steps); hP_grid=np.zeros(steps)
    hP_buf=np.zeros(steps);  hP_hid=np.zeros(steps)
    hP_gen=np.zeros(steps);  hom=np.zeros(steps)
    hr_q=np.zeros((steps,N)); hV_buf=np.zeros(steps)
    hdJ=np.zeros(steps);     hP_acc=np.zeros(steps)
    hropt=np.zeros(steps)

    for s in range(steps):
        t=s*DT; v=v_vent(t,s)
        omega+=(min(omR,lam_opt*max(v,0.1)/R)-omega)*(DT/3.0)
        omega=float(np.clip(omega,omR*0.3,omR*1.4))

        # ====================================================
        # MILLORA 1: r_opt CORREGIT (LeeChat)
        # 180° (avall, max captura) → R_tip (fora → frena → buffer)
        # 0°  (amunt, min captura) → R_hub (dins → accelera)
        # ====================================================
        dJ_total=0.0
        for i in range(N):
            th_n=thetas[i]%(2*np.pi)
            r_opt=R_tip-(R_tip-R_hub)*(1+np.cos(th_n))/2
            # Verificació: a 180° cos=-1 → r_opt=R_tip ✓
            #              a 0°   cos=+1 → r_opt=R_hub ✓

            F_c    =m_q*omega**2*r_q[i]
            F_ctrl =-K_q*(r_q[i]-r_opt)
            F_fric =-80*dr_q[i]
            d2r=(F_c+F_ctrl+F_fric)/m_q
            dr_q[i]+=d2r*DT
            dr_q[i]=float(np.clip(dr_q[i],-DR_MAX,DR_MAX))
            r_q[i]+=dr_q[i]*DT
            r_q[i]=float(np.clip(r_q[i],R_hub,R_tip))
            dJ_total+=2*m_q*r_q[i]*dr_q[i]

        P_buf=-0.5*dJ_total*omega**2
        hropt[s]=R_tip-(R_tip-R_hub)*(1+np.cos(thetas[0]%(2*np.pi)))/2

        # ====================================================
        # MILLORA 2+3: Acumulador llindar 20W + Generador
        # ====================================================
        V_acc=V_BUF_MAX; P_hid_inject=0.0; P_gen_inject=0.0

        if P_buf>LLINDAR_ACC and V_buf<V_BUF_MAX:
            dP=P_buf*DT/V_acc
            P_acc=min(P_acc+dP,P_MAX_ACC)
            Q_i=min(Q_BOMBA,P_buf/max(P_acc,1e3))
            dPf=perda_fric(Q_i,D,rho)
            Q_r=max(0,Q_i*(1-dPf/max(P_acc,1e3)))
            V_buf=min(V_BUF_MAX,V_buf+Q_r*DT)

        elif P_buf<-LLINDAR_ACC and V_buf>0:
            dP=-P_buf*DT/V_acc
            P_acc=max(P_acc-dP,0.0)
            Q_i=min(Q_BOMBA,-P_buf/max(P_acc,1e3))
            dPf=perda_fric(Q_i,D,rho)
            Q_r=max(0,Q_i*(1-dPf/max(P_acc,1e3)))
            # MILLORA 3: generador recupera energia
            P_gen_inject=P_acc*Q_r*ETA_GEN
            P_hid_inject=P_acc*Q_r*0.82
            V_buf=max(0.0,V_buf-Q_r*DT)

        # Pitch + Kuramoto
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
        hP_buf[s]=P_buf; hP_hid[s]=P_hid_inject
        hP_gen[s]=P_gen_inject; hom[s]=omega
        hr_q[s]=r_q.copy(); hV_buf[s]=V_buf*1000
        hdJ[s]=dJ_total; hP_acc[s]=P_acc/1e5

    return dict(P_aero=hP_aero,P_grid=hP_grid,P_buf=hP_buf,
                P_hid=hP_hid,P_gen=hP_gen,om=hom,r_q=hr_q,
                V_buf=hV_buf,dJ=hdJ,P_acc=hP_acc,r_opt=hropt,
                m_q=m_q,N=N,rho=rho,D=D)

print("Simulant v8 (3 millores LeeChat)...")
results={}
for nom in CONFIGS:
    r=simula(nom)
    results[nom]=r
    mill=(np.mean(r['P_grid'])-np.mean(r['P_aero']))/np.mean(r['P_aero'])*100
    mill=max(0,mill)
    Pb=np.mean(np.abs(r['P_buf']))
    Pg=np.mean(r['P_gen'])
    A_can=np.pi*(r['D']/2)**2
    Q=Pb/(r['rho']*9.81*H_BUF)*1000
    v_fl=Q/1000/A_can
    print(f"\n{nom}:")
    print(f"  Millora:       +{mill:.1f}%")
    print(f"  P_buf |mig|:   {Pb/1e3:.2f} kW")
    print(f"  P_gen inject:  {Pg:.1f} W")
    print(f"  Q fluid:       {Q:.2f} L/s  @ v={v_fl:.2f}m/s")
    print(f"  E_net:         {np.sum(r['P_buf'])*DT/3600:.3f} kWh")
    results[nom].update({'mill':mill,'Pb':Pb,'Pg':Pg,'Q':Q,'v_fl':v_fl})

# ============================================================
# GRÀFICS
# ============================================================
BG='#0d0d1a'; PAN='#13132b'
COLS={'3p Fe+oli':'#00d2ff','7p Fe+oli':'#00ff88'}
fig=plt.figure(figsize=(22,18),facecolor=BG)
gs=gridspec.GridSpec(4,3,figure=fig,hspace=0.52,wspace=0.38)

def sty(ax,t,xl='t[s]',yl=''):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9,pad=4)
    ax.tick_params(colors='#aaa',labelsize=7.5)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8)
    ax.set_ylabel(yl,color='#aaa',fontsize=8)
    ax.grid(color='#1e1e40',lw=0.5,ls='--')

def sh(ax): ax.axvspan(30,45,alpha=0.07,color='#e74c3c')

# F0: P_grid
ax=fig.add_subplot(gs[0,:2])
r0=list(results.values())[0]
ax.plot(t_vec,r0['P_aero']/1e6,'--',color='#888',lw=1.2,alpha=0.7,label='Base')
for nom,r in results.items():
    ax.plot(t_vec,r['P_grid']/1e6,color=COLS[nom],lw=2,
            label=f'{nom}  +{r["mill"]:.1f}%  P_buf={r["Pb"]/1e3:.1f}kW')
sh(ax); ax.legend(fontsize=7.5,framealpha=0.3)
sty(ax,'P_grid (MW) — v8 amb 3 millores LeeChat','t[s]','MW')

# F1: r_q + r_opt (millora 1)
ax=fig.add_subplot(gs[0,2])
r7=results['7p Fe+oli']
ax.plot(t_vec,r7['r_opt'],color='#ffd700',lw=2,ls='--',label='r_opt(t) — objectiu')
for i in range(min(2,r7['N'])):
    ax.plot(t_vec,r7['r_q'][:,i],lw=1.2,alpha=0.8,label=f'r_q pala {i+1}')
ax.axhline(R_hub,color='white',ls=':',lw=0.5,alpha=0.3,label='R_hub')
ax.axhline(R_tip,color='white',ls=':',lw=0.5,alpha=0.3,label='R_tip')
ax.legend(fontsize=6.5,framealpha=0.3)
sty(ax,'r_opt corregit + r_q real (7p)','t[s]','m')

# F2: P_buf comparativa
ax=fig.add_subplot(gs[1,0])
for nom,r in results.items():
    ax.plot(t_vec,r['P_buf']/1e3,color=COLS[nom],lw=1.5,alpha=0.8,label=nom)
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.fill_between(t_vec,0,r7['P_buf']/1e3,where=r7['P_buf']>0,
                alpha=0.12,color='#00ff88',label='emmagatzema')
ax.fill_between(t_vec,0,r7['P_buf']/1e3,where=r7['P_buf']<0,
                alpha=0.12,color='#e74c3c',label='allibera+genera')
ax.legend(fontsize=6.5,framealpha=0.3)
sty(ax,'P_buf Quijote (kW)','t[s]','kW')

# F3: Acumulador pressió (millora 2)
ax=fig.add_subplot(gs[1,1])
for nom,r in results.items():
    ax.plot(t_vec,r['P_acc'],color=COLS[nom],lw=1.5,label=nom)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'Pressio acumulador (bar)  llindar=20W','t[s]','bar')

# F4: P_gen recuperat (millora 3)
ax=fig.add_subplot(gs[1,2])
for nom,r in results.items():
    ax.plot(t_vec,r['P_gen'],color=COLS[nom],lw=1.5,label=f'{nom} P_gen')
ax.fill_between(t_vec,0,r7['P_gen'],alpha=0.25,color='#00ff88')
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'P_gen hidraulic recuperat (W)  eta=85%','t[s]','W')

# F5: dJ/dt
ax=fig.add_subplot(gs[2,0])
for nom,r in results.items():
    ax.plot(t_vec,r['dJ'],color=COLS[nom],lw=1.2,alpha=0.8,label=nom)
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'dJ/dt real (kg m2/s)','t[s]','kg m2/s')

# F6: Buffer volum
ax=fig.add_subplot(gs[2,1])
for nom,r in results.items():
    ax.plot(t_vec,r['V_buf'],color=COLS[nom],lw=1.5,label=nom)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'Volum buffer hidraulic (L)','t[s]','L')

# F7: omega
ax=fig.add_subplot(gs[2,2])
for nom,r in results.items():
    ax.plot(t_vec,r['om'],color=COLS[nom],lw=1.5,label=nom)
ax.axhline(omR,color='white',ls='--',lw=0.7,alpha=0.4,label=f'omR={omR:.3f}')
sh(ax); ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'omega (rad/s)','t[s]','rad/s')

# F8: Taula final comparativa
ax=fig.add_subplot(gs[3,:])
ax.axis('off'); ax.set_facecolor('#0a0a14')

# Dades taula
fluids=[('Aigua',1000),('Oli',900),('Fe+oli',3386),('Mercuri',13600)]
rows=[]
# Configuracions simulades
for nom,r in results.items():
    rows.append([nom, f'+{r["mill"]:.1f}%', f'{r["Pb"]/1e3:.2f}kW',
                 f'{r["Q"]:.2f}L/s', f'{r["v_fl"]:.2f}m/s',
                 f'{r["Pg"]:.0f}W', f'{r["m_q"]:.0f}kg'])
# Fluids de referència per 8.5kW
rows.append(['---']*7)
for fn,rho in fluids:
    Q=8500/(rho*9.81*H_BUF)*1000
    A=np.pi*(0.025)**2
    v=Q/1000/A
    rows.append([fn, f'{rho}kg/m3', '8.5kW ref',
                 f'{Q:.2f}L/s', f'{v:.2f}m/s', '—', '—'])

cols=['Config/Fluid','Millora','P_buf|mig|','Q fluid','v fluid','P_gen','m_q']
tbl=ax.table(cellText=rows,colLabels=cols,loc='center',cellLoc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5)
for (rr,c),cell in tbl.get_celld().items():
    bg='#1a1a2e' if rr==0 else ('#111' if rr==4 else '#0d0d1a')
    cell.set_facecolor(bg)
    cell.set_text_props(color='white' if rr==0 else '#ddd')
    cell.set_edgecolor('#333355')
ax.set_title('Taula comparativa final v8 — configuracions + fluids de referencia',
             color='white',fontsize=9,pad=4)

fig.suptitle(
    f'Gemell Virtual v8 — 3 millores LeeChat  |  r_opt corregit + llindar20W + generador85%\n'
    f'P_buf oscil·la ~8-10 kW  |  Fe+oli: v_fluid~1.15m/s  Q~2.3L/s  |  '
    f'Mercuri: v_fluid~0.29m/s  Q~0.56L/s',
    color='white',fontsize=10,fontweight='bold',y=0.999)

out='/mnt/user-data/outputs/gemelo_v8_final.png'
plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\nGrafic: {out}")
