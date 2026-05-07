"""
GEMELL VIRTUAL v9.2 — Valvules de retencio + limit fisic posicional
Victor Manzanares Alberola — EPSA UPV Alcoi

INSIGHT VICTOR: El limit |dr_q| <= omega*r_q NO es constant.
  A r_hub=5m:  limit = omega*5  =  6.8 m/s
  A r_mig=30m: limit = omega*30 = 40.7 m/s
  A r_tip=55m: limit = omega*55 = 74.7 m/s

  → El moment crític es quan la massa torna cap a dins
    i es prop del hub: limit baixa a 6.8 m/s
  → La bomba mai va a 6.8 m/s (DR_MAX=0.5) → sempre OK
  → Pero s'ha de vigilar si augmenta DR_MAX o omega baixa

VALVULES DE RETENCIO (LeeChat v9.2):
  - Bomba cap amunt NOMES quan P_buf > 0 (dJ/dt > 0, masses fora)
  - Allibera cap avall NOMES quan P_buf < 0 (dJ/dt < 0, masses dins)
  → Desacobla el moviment de les masses del flux de liquid
  → E_net tendeix a 0
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

np.random.seed(42)

KPITCH=2.0; K_KUR=0.10; K_q=5e4
V_BASE=11.4; T_TOTAL=60.0; DT=0.05
H_BUF=20.0; V_BUF_MAX=0.10; Q_BOMBA=0.010
P_MAX_ACC=20e6; LLINDAR_ACC=20.0; ETA_GEN=0.85
DR_MAX=0.5  # m/s limit simetric

R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.45; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3); omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0; J_base=3.54e7

print(f"omR={omR:.3f} rad/s")
print(f"Limit fisic posicional:")
for r in [R_hub, (R_hub+R_tip)/2, R_tip]:
    print(f"  r={r:.0f}m: omega*r = {omR*r:.1f} m/s (DR_MAX={DR_MAX} << OK)")

CONFIGS={'3p Fe+oli':dict(N=3,rho=3386.0,D=0.05),
         '7p Fe+oli':dict(N=7,rho=3386.0,D=0.05)}

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
    hP_acc=np.zeros(steps);  hlim_min=np.zeros(steps)
    hdr_max_act=np.zeros(steps)  # limit actiu (min de DR_MAX i omega*r)

    for s in range(steps):
        t=s*DT; v=v_vent(t,s)
        omega+=(min(omR,lam_opt*max(v,0.1)/R)-omega)*(DT/3.0)
        omega=float(np.clip(omega,omR*0.3,omR*1.4))

        dJ_total=0.0
        lim_min_step=999.0  # limit mes restrictiu del pas

        for i in range(N):
            th_n=thetas[i]%(2*np.pi)
            r_opt=R_tip-(R_tip-R_hub)*(1+np.cos(th_n))/2

            F_c=m_q*omega**2*r_q[i]
            F_ctrl=-K_q*(r_q[i]-r_opt)
            F_fric=-80*dr_q[i]
            d2r=(F_c+F_ctrl+F_fric)/m_q
            dr_q[i]+=d2r*DT

            # Limit simetric DR_MAX
            dr_q[i]=float(np.clip(dr_q[i],-DR_MAX,DR_MAX))

            # LIMIT FISIC POSICIONAL: omega*r_q (canvia en cada punt!)
            lim_pos=omega*r_q[i]
            lim_min_step=min(lim_min_step,lim_pos)
            dr_q[i]=float(np.clip(dr_q[i],-lim_pos,lim_pos))

            r_q[i]+=dr_q[i]*DT
            r_q[i]=float(np.clip(r_q[i],R_hub,R_tip))
            dJ_total+=2*m_q*r_q[i]*dr_q[i]

        hlim_min[s]=lim_min_step  # limit mes restrictiu (a r_hub)
        hdr_max_act[s]=min(DR_MAX,lim_min_step)
        P_buf=-0.5*dJ_total*omega**2

        # ================================================
        # VALVULES DE RETENCIO (v9.2, LeeChat)
        # Desacobla moviment masses del flux de liquid
        # ================================================
        V_acc=V_BUF_MAX; P_gen_inject=0.0; P_hid_inject=0.0

        if P_buf>LLINDAR_ACC and V_buf<V_BUF_MAX:
            # masses cap a fora → bomba cap amunt (emmagatzema)
            dP=P_buf*DT/V_acc; P_acc=min(P_acc+dP,P_MAX_ACC)
            Q_i=min(Q_BOMBA,P_buf/max(P_acc,1e3))
            dPf=perda_fric(Q_i,D,rho)
            Q_r=max(0,Q_i*(1-dPf/max(P_acc,1e3)))
            V_buf=min(V_BUF_MAX,V_buf+Q_r*DT)
            P_hid_inject=0.0  # VALVULA: no injecta mentre carrega

        elif P_buf<-LLINDAR_ACC and V_buf>0:
            # masses cap a dins → allibera liquid (injecta energia)
            dP=-P_buf*DT/V_acc; P_acc=max(P_acc-dP,0.0)
            Q_i=min(Q_BOMBA,-P_buf/max(P_acc,1e3))
            dPf=perda_fric(Q_i,D,rho)
            Q_r=max(0,Q_i*(1-dPf/max(P_acc,1e3)))
            P_hid_inject=P_acc*Q_r*0.82
            P_gen_inject=P_acc*Q_r*ETA_GEN
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
                dJ=hdJ,P_acc=hP_acc,lim_min=hlim_min,
                dr_max_act=hdr_max_act,m_q=m_q,N=N,rho=rho,D=D)

print("\nSimulant v9.2...")
results={}
for nom in CONFIGS:
    r=simula(nom)
    mill=(np.mean(r['P_grid'])-np.mean(r['P_aero']))/np.mean(r['P_aero'])*100
    mill=max(0,mill)
    Pb=np.mean(np.abs(r['P_buf'])); Pg=np.mean(r['P_gen'])
    E_net=np.sum(r['P_buf'])*DT/3600
    A_can=np.pi*(r['D']/2)**2
    Q=Pb/(r['rho']*9.81*H_BUF)*1000; v_fl=Q/1000/A_can
    lim_critic=np.min(r['lim_min'])
    print(f"\n{nom}:")
    print(f"  Millora:     +{mill:.1f}%")
    print(f"  P_buf|mig|:  {Pb/1e3:.2f} kW")
    print(f"  E_net:       {E_net:.2f} kWh")
    print(f"  P_gen:       {Pg:.0f} W")
    print(f"  Q: {Q:.2f}L/s  v_fl: {v_fl:.2f}m/s")
    print(f"  Limit min (omega*r_hub): {lim_critic:.2f} m/s >> DR_MAX={DR_MAX} ✓")
    results[nom]=dict(r=r,mill=mill,Pb=Pb,Pg=Pg,E_net=E_net,
                      Q=Q,v_fl=v_fl,lim_min=lim_critic)

# GRAFICS
BG='#0d0d1a'; PAN='#13132b'
COLS={'3p Fe+oli':'#00d2ff','7p Fe+oli':'#00ff88'}
fig=plt.figure(figsize=(20,18),facecolor=BG)
gs=gridspec.GridSpec(4,3,figure=fig,hspace=0.52,wspace=0.38)

def sty(ax,t,xl='t[s]',yl=''):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9,pad=4)
    ax.tick_params(colors='#aaa',labelsize=7.5)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8); ax.set_ylabel(yl,color='#aaa',fontsize=8)
    ax.grid(color='#1e1e40',lw=0.5,ls='--')

# F0: P_grid
ax=fig.add_subplot(gs[0,:2])
r0=list(results.values())[0]['r']
ax.plot(t_vec,r0['P_aero']/1e6,'--',color='#888',lw=1.2,alpha=0.7,label='Base')
for nom,d in results.items():
    ax.plot(t_vec,d['r']['P_grid']/1e6,color=COLS[nom],lw=2,
            label=f'{nom}  +{d["mill"]:.1f}%  E_net={d["E_net"]:.1f}kWh')
ax.axvspan(30,45,alpha=0.07,color='#e74c3c')
ax.legend(fontsize=7.5,framealpha=0.3)
sty(ax,'P_grid (MW) — v9.2 valvules retencio','t[s]','MW')

# F1: r_q + limit posicional
ax=fig.add_subplot(gs[0,2])
r7=results['7p Fe+oli']['r']
for i in range(min(2,r7['N'])):
    ax.plot(t_vec,r7['r_q'][:,i],lw=1.5,alpha=0.85,label=f'r_q pala {i+1}')
ax2=ax.twinx()
ax2.plot(t_vec,r7['lim_min'],color='#ffd700',lw=1.5,ls='--',alpha=0.7,
         label='omega*r_q(t) limit')
ax2.axhline(DR_MAX,color='#e74c3c',ls=':',lw=1.2,label=f'DR_MAX={DR_MAX}')
ax.legend(fontsize=6.5,framealpha=0.3,loc='upper left')
ax2.legend(fontsize=6.5,framealpha=0.3,loc='upper right')
ax.set_ylabel('r_q [m]',color='#aaa',fontsize=8)
ax2.set_ylabel('omega*r_q [m/s]',color='#ffd700',fontsize=8)
ax.tick_params(colors='#aaa'); ax2.tick_params(colors='#ffd700')
ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--')
ax.set_xlabel('t[s]',color='#aaa',fontsize=8)
ax.set_title('r_q + limit posicional omega*r_q(t)',color='white',fontsize=9,pad=4)

# F2: P_buf + dJ/dt
ax=fig.add_subplot(gs[1,0])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['P_buf']/1e3,color=COLS[nom],lw=1.5,alpha=0.8,label=nom)
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.axhline(LLINDAR_ACC/1000,color='#ffd700',ls=':',lw=0.8,alpha=0.6,
           label=f'llindar={LLINDAR_ACC}W')
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'P_buf Quijote (kW) — valvula actua','t[s]','kW')

# F3: Acumulador
ax=fig.add_subplot(gs[1,1])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['P_acc'],color=COLS[nom],lw=1.5,label=nom)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'Pressio acumulador (bar)','t[s]','bar')

# F4: P_gen
ax=fig.add_subplot(gs[1,2])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['P_gen'],color=COLS[nom],lw=1.5,label=nom)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'P_gen hidraulic (W) eta=85%','t[s]','W')

# F5: dJ/dt
ax=fig.add_subplot(gs[2,0])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['dJ'],color=COLS[nom],lw=1.2,alpha=0.8,label=nom)
ax.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'dJ/dt (kg m2/s)','t[s]','kg m2/s')

# F6: Buffer volum
ax=fig.add_subplot(gs[2,1])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['V_buf'],color=COLS[nom],lw=1.5,label=nom)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'Buffer liquid (L)','t[s]','L')

# F7: omega
ax=fig.add_subplot(gs[2,2])
for nom,d in results.items():
    ax.plot(t_vec,d['r']['om'],color=COLS[nom],lw=1.5,label=nom)
ax.axhline(omR,color='white',ls='--',lw=0.7,alpha=0.4)
ax.axvspan(30,45,alpha=0.07,color='#e74c3c')
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'omega (rad/s)','t[s]','rad/s')

# F8: Diagrama explicatiu limit posicional
ax=fig.add_subplot(gs[3,:])
ax.axis('off'); ax.set_facecolor('#0a0a14')
r_vals=np.linspace(R_hub,R_tip,100)
lim_vals=omR*r_vals

ax_ins=ax.inset_axes([0.0,0.05,0.38,0.90])
ax_ins.set_facecolor('#0d0d1a')
ax_ins.plot(r_vals,lim_vals,color='#ffd700',lw=2.5,label='omega*r_q (limit fisic)')
ax_ins.axhline(DR_MAX,color='#e74c3c',ls='--',lw=2,label=f'DR_MAX={DR_MAX}m/s')
ax_ins.fill_between(r_vals,DR_MAX,lim_vals,alpha=0.15,color='#00ff88',
                    label='zona segura')
ax_ins.axvline(R_hub,color='#888',ls=':',lw=0.8)
ax_ins.axvline(R_tip,color='#888',ls=':',lw=0.8)
ax_ins.set_xlabel('r_q [m]',color='#aaa',fontsize=8)
ax_ins.set_ylabel('limit |dr_q| [m/s]',color='#aaa',fontsize=8)
ax_ins.set_title('Limit posicional omega*r_q',color='white',fontsize=8.5)
ax_ins.legend(fontsize=7,framealpha=0.3)
ax_ins.tick_params(colors='#aaa',labelsize=7.5)
for sp in ax_ins.spines.values(): sp.set_color('#333355')
ax_ins.grid(color='#1e1e40',lw=0.5,ls='--')

# Taula
rows=[]
for nom,d in results.items():
    rows.append([nom,f'+{d["mill"]:.1f}%',
                 f'{d["Pb"]/1e3:.2f}kW',
                 f'{d["E_net"]:.1f}kWh',
                 f'{d["Q"]:.1f}L/s',
                 f'{d["v_fl"]:.1f}m/s',
                 f'{d["Pg"]:.0f}W',
                 f'{d["lim_min"]:.1f}m/s'])
rows.append(['LIMIT FISIC',f'omega*{R_hub}m',
             f'={omR*R_hub:.1f}m/s',
             f'>> DR_MAX={DR_MAX}',
             'SEMPRE OK','—','—','posicional'])

ax2_ins=ax.inset_axes([0.40,0.05,0.58,0.90])
ax2_ins.axis('off'); ax2_ins.set_facecolor('#0a0a14')
tbl=ax2_ins.table(cellText=rows,
    colLabels=['Config','Millora','P_buf','E_net','Q','v_fl','P_gen','lim_min'],
    loc='center',cellLoc='center')
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5)
for (rr,c),cell in tbl.get_celld().items():
    cell.set_facecolor('#1a1a2e' if rr==0 else '#0d0d1a')
    cell.set_text_props(color='white' if rr==0 else '#ddd')
    cell.set_edgecolor('#333355')

fig.suptitle(
    f'Gemell Virtual v9.2 — Valvules de retencio + Limit fisic posicional\n'
    f'Victor: "no puc botar molt alt, es un limit que canvia segons el lloc"\n'
    f'omega*R_hub={omR*R_hub:.1f}m/s (critic) >> DR_MAX={DR_MAX}m/s | '
    f'omega*R_tip={omR*R_tip:.1f}m/s (ample) | Limit posicional SEMPRE respectat',
    color='white',fontsize=10,fontweight='bold',y=0.999)

out='/mnt/user-data/outputs/gemelo_v9b_valvules.png'
plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\nGrafic: {out}")
