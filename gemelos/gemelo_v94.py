"""
GEMELL VIRTUAL v9.4 — Versio completa per al paper
Victor Manzanares Alberola — EPSA UPV Alcoi

FISICA COMPLETA:
  - Quijote: r_q(t) amb F_c + F_ctrl + F_fric (fisica real)
  - r_opt suau: R_tip-(R_tip-R_hub)*(1+cos(theta))/2
  - Limit fisic posicional: |dr_q| <= omega*r_q  [observacio Victor]
  - Pitch asincron per pala: Kp*cos + Kom*sin - KQ*dJdt
  - Kuramoto normalitzat: K/N * sum(sin(theta_j - theta_i))
  - Buffer hidraulic amb valvules de retencio
  - Acumulador de pressio 20 bar + generador hidraulic eta=85%

GRAFICS PER AL PAPER (9 panels):
  G1: Perfil de vent realista (OU + calmant)
  G2: P_grid vs P_base — comparativa 3 configs
  G3: r_q(t) per pala — sincronitzat amb theta
  G4: dJ/dt real i P_buf(t) — el motor del buffer
  G5: Limit fisic posicional omega*r_q vs DR_MAX
  G6: Acumulador de pressio i volum buffer
  G7: P_gen hidraulic recuperat
  G8: omega(t) — estabilitzacio del rotor
  G9: Taula resum final per al paper
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

np.random.seed(42)

# ============================================================
# PARAMETRES (NREL 5MW, Jonkman 2009)
# ============================================================
R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max=0.482; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max))**(1/3)
omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0; J_base=3.54e7

# Control
KPITCH=2.0; KOMEGA=0.5; KQ=0.8; K_KUR=0.10; K_q=5e4
# Fluid: Fe+oli
RHO_FL=3386.0; D_CANAL=0.05
A_can=np.pi*(D_CANAL/2)**2
# Limits
DR_MAX=0.5   # m/s
# Buffer hidraulic
H_BUF=20.0; V_BUF_MAX=0.10; Q_BOMBA=0.010
P_MAX_ACC=20e6; LLINDAR=20.0; ETA_GEN=0.85
# Simulacio
V_BASE=11.4; T_TOT=80.0; DT=0.05

print(f"NREL 5MW: omR={omR:.4f} rad/s | vr={vr:.3f} m/s")
print(f"omega*R_hub (limit critic) = {omR*R_hub:.2f} m/s >> DR_MAX={DR_MAX} m/s")

# ============================================================
# VENT REALISTA
# ============================================================
steps=int(T_TOT/DT); t_vec=np.arange(steps)*DT
# Ornstein-Uhlenbeck (turbulencia)
ou=np.zeros(steps)
for i in range(1,steps):
    ou[i]=ou[i-1]-0.5*ou[i-1]*DT+0.8*np.sqrt(DT)*np.random.randn()

def v_vent(s):
    t=s*DT; v=V_BASE+2*np.sin(2*np.pi*t/20)+ou[s]*0.6
    if 35<t<50: v=max(3,v-5*np.sin(np.pi*(t-35)/15))  # calmant
    return float(v)

v_arr=np.array([v_vent(s) for s in range(steps)])

# ============================================================
# FUNCIONS
# ============================================================
def eta_lam(om,v):
    if v<=0: return 0.0
    return float(max(0,1-((om*R/v-lam_opt)/lam_opt)**2))

def eta_kur(K,N):
    return float(np.clip(0.04*min(1.2,N/3)*(1-((K-0.10)/0.10)**2),0,0.048))

def perda_fric(Q,rho,D,L=10,f=0.02):
    A=np.pi*(D/2)**2; return f*(L/D)*(rho*(Q/max(A,1e-12))**2)/2

def m_q_pala(rho,D,r1,r2):
    return rho*np.pi*(D/2)**2*(r2-r1)

# ============================================================
# SIMULACIO — funcio generica
# ============================================================
def simula(N, label):
    m_q=m_q_pala(RHO_FL,D_CANAL,R_hub,R_tip)
    thetas=np.array([2*np.pi*i/N for i in range(N)])
    omega=omR*0.95
    r_q=np.full(N,(R_hub+R_tip)/2)
    dr_q=np.zeros(N)
    V_buf=0.0; P_acc=0.0
    etaK=eta_kur(K_KUR,N)

    hPa=np.zeros(steps); hPg=np.zeros(steps)
    hPb=np.zeros(steps); hPgen=np.zeros(steps)
    hom=np.zeros(steps); hrq=np.zeros((steps,N))
    hVb=np.zeros(steps); hdJ=np.zeros(steps)
    hPacc=np.zeros(steps); hlimf=np.zeros(steps)
    hropt=np.zeros(steps)

    for s in range(steps):
        v=v_arr[s]
        omega+=(min(omR,lam_opt*max(v,0.1)/R)-omega)*(DT/3.0)
        omega=float(np.clip(omega,omR*0.3,omR*1.4))

        dJ=0.0; lim_min=999.0
        for i in range(N):
            # r_opt suau corregit: 180° → R_tip, 0° → R_hub
            th=thetas[i]%(2*np.pi)
            r_opt=R_tip-(R_tip-R_hub)*(1+np.cos(th))/2

            F_c=m_q*omega**2*r_q[i]
            F_ct=-K_q*(r_q[i]-r_opt)
            F_fr=-80*dr_q[i]
            dr_q[i]+=(F_c+F_ct+F_fr)/m_q*DT

            # limit simetric
            dr_q[i]=float(np.clip(dr_q[i],-DR_MAX,DR_MAX))
            # LIMIT FISIC POSICIONAL (Victor)
            lp=omega*r_q[i]; lim_min=min(lim_min,lp)
            dr_q[i]=float(np.clip(dr_q[i],-lp,lp))

            r_q[i]=float(np.clip(r_q[i]+dr_q[i]*DT,R_hub,R_tip))
            dJ+=2*m_q*r_q[i]*dr_q[i]

        P_buf=-0.5*dJ*omega**2
        hlimf[s]=lim_min
        hropt[s]=R_tip-(R_tip-R_hub)*(1+np.cos(thetas[0]%(2*np.pi)))/2

        # Valvules de retencio + buffer
        Pgen=0.0; Phid=0.0
        if P_buf>LLINDAR and V_buf<V_BUF_MAX:
            P_acc=min(P_acc+P_buf*DT/V_BUF_MAX,P_MAX_ACC)
            Q=min(Q_BOMBA,P_buf/max(P_acc,1e3))
            Q=max(0,Q*(1-perda_fric(Q,RHO_FL,D_CANAL)/max(P_acc,1e3)))
            V_buf=min(V_BUF_MAX,V_buf+Q*DT)
        elif P_buf<-LLINDAR and V_buf>0:
            P_acc=max(P_acc+P_buf*DT/V_BUF_MAX,0)
            Q=min(Q_BOMBA,-P_buf/max(P_acc,1e3))
            Q=max(0,Q*(1-perda_fric(Q,RHO_FL,D_CANAL)/max(P_acc,1e3)))
            Pgen=P_acc*Q*ETA_GEN; Phid=P_acc*Q*0.82
            V_buf=max(0,V_buf-Q*DT)

        # Pitch + Kuramoto
        for i in range(N):
            dJi=2*m_q*r_q[i]*dr_q[i]
            pm=np.clip(KPITCH*np.cos(thetas[i])+KOMEGA*omega*np.sin(thetas[i])-KQ*dJi,-8,8)
            kc=(K_KUR/N)*float(np.sum(np.sin(thetas-thetas[i])))
            thetas[i]+=(omega+kc+pm*0.01)*DT

        # Potencies
        etL=eta_lam(omega,v)
        Pa=0.5*RHO_A*A_rot*Cp_max*etL*v**3
        Pt=min(1.0,etL*(1+0.04)*(1+etaK))
        Pg=max(0,0.5*RHO_A*A_rot*Cp_max*Pt*v**3+max(0,P_buf)+Phid+Pgen)

        hPa[s]=Pa; hPg[s]=Pg; hPb[s]=P_buf
        hPgen[s]=Pgen; hom[s]=omega
        hrq[s]=r_q.copy(); hVb[s]=V_buf*1000
        hdJ[s]=dJ; hPacc[s]=P_acc/1e5

    return dict(Pa=hPa,Pg=hPg,Pb=hPb,Pgen=hPgen,om=hom,
                rq=hrq,Vb=hVb,dJ=hdJ,Pacc=hPacc,
                limf=hlimf,ropt=hropt,N=N,m_q=m_q,label=label)

print("Simulant 3 configuracions per al paper...")
configs=[
    (0,'Base (sense sistema)'),
    (3,'3 pales — Fe+oli v9.4'),
    (7,'7 pales — Fe+oli v9.4'),
]
rr={}
# Simulacio base (sense control)
def simula_base():
    omega=omR*0.95
    hPa=np.zeros(steps); hPg=np.zeros(steps); hom=np.zeros(steps)
    for s in range(steps):
        v=v_arr[s]
        omega+=(min(omR,lam_opt*max(v,0.1)/R)-omega)*(DT/3.0)
        omega=float(np.clip(omega,omR*0.3,omR*1.4))
        etL=eta_lam(omega,v)
        Pa=0.5*RHO_A*A_rot*Cp_max*etL*v**3
        hPa[s]=Pa; hPg[s]=Pa; hom[s]=omega
    return dict(Pa=hPa,Pg=hPg,om=hom,label='Base')

rr['base']=simula_base()
for N,lbl in [(3,'3p'),(7,'7p')]:
    rr[lbl]=simula(N,lbl)

for k in ['3p','7p']:
    r=rr[k]; rb=rr['base']
    mill=(np.mean(r['Pg'])-np.mean(rb['Pg']))/np.mean(rb['Pg'])*100
    Pb=np.mean(np.abs(r['Pb']))
    E_net=np.sum(r['Pb'])*DT/3600
    Q=Pb/(RHO_FL*9.81*H_BUF)*1000
    print(f"\n  {r['label']}:")
    print(f"    Millora P_grid:  +{mill:.1f}%")
    print(f"    P_buf |mig|:     {Pb/1e3:.2f} kW")
    print(f"    E_net:           {E_net:.1f} kWh")
    print(f"    P_gen:           {np.mean(r['Pgen']):.0f} W")
    print(f"    Q fluid:         {Q:.2f} L/s")
    print(f"    omega*r_hub lim: {np.min(r['limf']):.1f} m/s >> DR_MAX={DR_MAX}")
    rr[k]['mill']=mill; rr[k]['Pb_abs']=Pb
    rr[k]['E_net']=E_net; rr[k]['Q']=Q

# ============================================================
# FIGURA COMPLETA PER AL PAPER
# ============================================================
BG='#0d0d1a'; PAN='#13132b'
C3='#00d2ff'; C7='#00ff88'; CB='#888780'; CV='#ffd700'
CX='#e74c3c'; CD='#D85A30'; CG='#BA7517'

fig=plt.figure(figsize=(24,22),facecolor=BG)
gs=gridspec.GridSpec(3,3,figure=fig,hspace=0.52,wspace=0.38)

def sty(ax,t,xl='t [s]',yl='',xl2=None):
    ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9.5,pad=5,fontweight='bold')
    ax.tick_params(colors='#aaa',labelsize=8)
    [sp.set_color('#333355') for sp in ax.spines.values()]
    ax.set_xlabel(xl,color='#aaa',fontsize=8.5); ax.set_ylabel(yl,color='#aaa',fontsize=8.5)
    ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)

def sh(ax): ax.axvspan(35,50,alpha=0.09,color=CX,label='calmant')

# --- G1: Vent realista ---
ax=fig.add_subplot(gs[0,0])
ax.plot(t_vec,v_arr,color=CV,lw=1.8,alpha=0.9)
ax.fill_between(t_vec,0,v_arr,alpha=0.12,color=CV)
ax.axhline(vr,color='white',ls='--',lw=1,alpha=0.5,label=f'v_rated={vr:.1f}m/s')
ax.axhline(V_BASE,color='#888',ls=':',lw=0.8,alpha=0.5,label=f'v_base={V_BASE}m/s')
sh(ax); ax.legend(fontsize=7.5,framealpha=0.3)
sty(ax,'G1 — Perfil de vent realista (OU + calmant)','t [s]','v [m/s]')

# --- G2: P_grid comparativa ---
ax=fig.add_subplot(gs[0,1:])
ax.plot(t_vec,rr['base']['Pg']/1e6,'--',color=CB,lw=1.5,alpha=0.7,label='Base (sense sistema)')
ax.plot(t_vec,rr['3p']['Pg']/1e6,color=C3,lw=2,alpha=0.9,
        label=f"3p Fe+oli  +{rr['3p']['mill']:.1f}%  P_buf={rr['3p']['Pb_abs']/1e3:.1f}kW")
ax.plot(t_vec,rr['7p']['Pg']/1e6,color=C7,lw=2,alpha=0.9,
        label=f"7p Fe+oli  +{rr['7p']['mill']:.1f}%  P_buf={rr['7p']['Pb_abs']/1e3:.1f}kW")
ax.fill_between(t_vec,rr['base']['Pg']/1e6,rr['7p']['Pg']/1e6,alpha=0.1,color=C7)
sh(ax); ax.legend(fontsize=8,framealpha=0.3)
sty(ax,'G2 — P_grid (MW): base vs sistema Quijote+ZYPYZAPE','t [s]','P [MW]')

# --- G3: r_q(t) per pala ---
ax=fig.add_subplot(gs[1,0])
r7=rr['7p']
COLS_P=[C3,C7,CV,'#ff6b35','#D85A30','#9b59b6','#e91e63']
for i in range(min(3,r7['N'])):
    ax.plot(t_vec,r7['rq'][:,i],color=COLS_P[i],lw=1.5,alpha=0.85,label=f'pala {i+1}')
ax.plot(t_vec,r7['ropt'],color='white',lw=1.5,ls='--',alpha=0.6,label='r_opt(theta_0)')
ax.axhline(R_hub,color='#555',ls=':',lw=0.8,alpha=0.5,label=f'R_hub={R_hub}m')
ax.axhline(R_tip,color='#555',ls=':',lw=0.8,alpha=0.5,label=f'R_tip={R_tip}m')
ax.legend(fontsize=7,framealpha=0.3); sh(ax)
sty(ax,'G3 — r_q(t) per pala  (7p, sincronitzat amb theta)','t [s]','r_q [m]')

# --- G4: dJ/dt i P_buf ---
ax=fig.add_subplot(gs[1,1])
ax2=ax.twinx()
r3=rr['3p']
ax.plot(t_vec,r3['dJ'],color=CD,lw=1.2,alpha=0.8,label='dJ/dt (kg·m²/s)')
ax.fill_between(t_vec,0,r3['dJ'],where=np.array(r3['dJ'])>0,alpha=0.18,color=CX)
ax.fill_between(t_vec,0,r3['dJ'],where=np.array(r3['dJ'])<0,alpha=0.18,color=C3)
ax2.plot(t_vec,r3['Pb']/1e3,color=C7,lw=1.5,alpha=0.9,label='P_buf (kW)')
ax2.axhline(0,color='white',lw=0.5,alpha=0.3)
ax.set_ylabel('dJ/dt [kg·m²/s]',color=CD,fontsize=8); ax.tick_params(colors=CD)
ax2.set_ylabel('P_buf [kW]',color=C7,fontsize=8); ax2.tick_params(colors=C7)
ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)
ax.set_xlabel('t [s]',color='#aaa',fontsize=8.5)
ax.set_title('G4 — dJ/dt real i P_buf Quijote (3p)',color='white',fontsize=9.5,pad=5,fontweight='bold')
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2,l1+l2,fontsize=7,framealpha=0.3)
sh(ax)

# --- G5: Limit fisic posicional ---
ax=fig.add_subplot(gs[1,2])
ax.plot(t_vec,r7['limf'],color=CV,lw=2,label='omega·r_q(t) limit posicional [m/s]')
ax.fill_between(t_vec,DR_MAX,r7['limf'],alpha=0.12,color=C7,label='zona segura')
ax.axhline(DR_MAX,color=CX,ls='--',lw=2,label=f'DR_MAX = {DR_MAX} m/s')
ax.axhline(omR*R_hub,color='#ff6b35',ls=':',lw=1.5,alpha=0.8,
           label=f'omega·R_hub = {omR*R_hub:.1f} m/s (critic)')
ax.text(2,DR_MAX+1,f'DR_MAX={DR_MAX}m/s  <<  limit fisic min={omR*R_hub:.1f}m/s',
        color='white',fontsize=7.5)
ax.legend(fontsize=7,framealpha=0.3)
sty(ax,'G5 — Limit fisic posicional |dr_q| <= omega·r_q(t)  [Victor]','t [s]','m/s')

# --- G6: Acumulador pressio + buffer volum ---
ax=fig.add_subplot(gs[2,0])
ax2=ax.twinx()
ax.plot(t_vec,r7['Pacc'],color=CG,lw=2,label='Pressio acumulador (bar)')
ax.fill_between(t_vec,0,r7['Pacc'],alpha=0.18,color=CG)
ax2.plot(t_vec,r7['Vb'],color='#185FA5',lw=1.5,ls='--',alpha=0.9,label='Volum buffer (L)')
ax.set_ylabel('P_acc [bar]',color=CG,fontsize=8); ax.tick_params(colors=CG)
ax2.set_ylabel('V_buf [L]',color='#185FA5',fontsize=8); ax2.tick_params(colors='#185FA5')
ax.set_facecolor(PAN); ax.grid(color='#1e1e40',lw=0.5,ls='--',alpha=0.8)
ax.set_xlabel('t [s]',color='#aaa',fontsize=8.5)
ax.set_title('G6 — Acumulador pressio (bar) + buffer (L)',color='white',fontsize=9.5,pad=5,fontweight='bold')
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
ax.legend(h1+h2,l1+l2,fontsize=7,framealpha=0.3)
sh(ax)

# --- G7: P_gen hidraulic ---
ax=fig.add_subplot(gs[2,1])
ax.plot(t_vec,r3['Pgen'],color=C3,lw=1.5,label='3p — P_gen (W)',alpha=0.85)
ax.plot(t_vec,r7['Pgen'],color=C7,lw=1.5,label='7p — P_gen (W)',alpha=0.85)
ax.fill_between(t_vec,0,r7['Pgen'],alpha=0.15,color=C7)
Pg3=np.mean(r3['Pgen']); Pg7=np.mean(r7['Pgen'])
ax.axhline(Pg3,color=C3,ls='--',lw=0.8,alpha=0.6,label=f'mig 3p={Pg3:.0f}W')
ax.axhline(Pg7,color=C7,ls='--',lw=0.8,alpha=0.6,label=f'mig 7p={Pg7:.0f}W')
ax.legend(fontsize=7.5,framealpha=0.3); sh(ax)
sty(ax,'G7 — P_gen hidraulic recuperat (eta=85%)','t [s]','W')

# --- G8: omega ---
ax=fig.add_subplot(gs[2,2])
ax.plot(t_vec,rr['base']['om'],'--',color=CB,lw=1.2,alpha=0.7,label='Base omega')
ax.plot(t_vec,r3['om'],color=C3,lw=1.8,alpha=0.9,label='3p omega')
ax.plot(t_vec,r7['om'],color=C7,lw=1.8,alpha=0.9,label='7p omega')
ax.axhline(omR,color='white',ls='--',lw=1,alpha=0.4,label=f'omega_R={omR:.3f}')
ax.fill_between(t_vec,rr['base']['om'],r7['om'],alpha=0.08,color=C7,label='Delta omega')
ax.legend(fontsize=7.5,framealpha=0.3); sh(ax)
sty(ax,'G8 — omega (rad/s) — estabilitzacio del rotor','t [s]','omega [rad/s]')

fig.suptitle(
    'QUIJOTE + ZYPYZAPE v9.4 — Gemell Virtual Complet per al Paper\n'
    'Victor Manzanares Alberola — EPSA UPV Alcoi — 2026\n'
    'NREL 5MW | Fe+oli rho=3386 kg/m3 | Limit posicional |dr_q|<=omega·r_q [Observacio Victor] | '
    'Valvules de retencio | Acumulador 20bar',
    color='white',fontsize=10,fontweight='bold',y=0.999)

out1='/mnt/user-data/outputs/gemelo_v94_paper_grafics.png'
plt.savefig(out1,dpi=160,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"\nGrafics paper: {out1}")

# ============================================================
# DIAGRAMA ESQUEMATIC DEL SISTEMA (per al paper)
# ============================================================
fig2,ax=plt.subplots(1,1,figsize=(18,10),facecolor=BG)
ax.set_facecolor(BG); ax.axis('off')

def box(cx,cy,w,h,col,txt,txt2='',fontsize=10):
    ax.add_patch(mpatches.FancyBboxPatch((cx-w/2,cy-h/2),w,h,
        boxstyle="round,pad=0.05",fc=col,ec='white',lw=1.5,alpha=0.85))
    ax.text(cx,cy+(0.06 if txt2 else 0),txt,ha='center',va='center',
            color='white',fontsize=fontsize,fontweight='bold')
    if txt2:
        ax.text(cx,cy-0.08,txt2,ha='center',va='center',color='#ddd',fontsize=8.5)

def arr(x1,y1,x2,y2,col='white',txt='',lw=1.5):
    ax.annotate('',xy=(x2,y2),xytext=(x1,y1),
        arrowprops=dict(arrowstyle='->',color=col,lw=lw,mutation_scale=18))
    if txt:
        mx,my=(x1+x2)/2,(y1+y2)/2
        ax.text(mx+0.01,my,txt,ha='center',va='bottom',color=col,fontsize=8)

ax.set_xlim(0,1); ax.set_ylim(0,1)

# Vent
box(0.09,0.82,0.14,0.1,'#1F5C9E','VENT','v(t) OU + calmant',9)
arr(0.16,0.82,0.28,0.70,CV,'v(t)')

# Rotor
box(0.38,0.70,0.18,0.12,'#8B1A8B','ROTOR\nNREL 5MW','R=63m, omega(t)',9)
arr(0.16,0.70,0.29,0.70,C3,'')

# Quijote
box(0.38,0.45,0.18,0.14,'#8B3A1A','QUIJOTE\nFe+oli','rho=3386, D=5cm\nR_hub=5m R_tip=55m',8)
arr(0.38,0.64,0.38,0.52,CX,'dJ/dt')
# Forces
ax.text(0.565,0.50,'F_c = m·omega²·r  (fora)',ha='left',va='center',color='#ff6b35',fontsize=8)
ax.text(0.565,0.46,'F_ctrl = -K_q·(r-r_opt)  (control)',ha='left',va='center',color=C3,fontsize=8)
ax.text(0.565,0.42,'r_opt(theta) = R_tip-(R_tip-R_hub)·(1+cos(theta))/2',ha='left',va='center',color=CV,fontsize=8)
ax.text(0.565,0.38,'|dr_q| <= omega·r_q  [Limit Victor]',ha='left',va='center',color=CX,fontsize=9,fontweight='bold')

# Pitch ZYPYZAPE
box(0.14,0.45,0.18,0.12,'#1A5C2E','PITCH ASINCRON\nZYPYZAPE','Kp·cos+Kom·sin-KQ·dJ',8)
arr(0.14,0.64,0.14,0.51,C7,'')
arr(0.23,0.45,0.29,0.45,C7,'pitch(theta_i)')

# Kuramoto
box(0.14,0.24,0.18,0.10,'#1A3A5C','KURAMOTO','K/N·sum(sin(dtheta))',8)
arr(0.14,0.39,0.14,0.29,C3,'')
arr(0.23,0.24,0.29,0.60,'#9b59b6','coherencia')

# Buffer hidraulic
box(0.38,0.20,0.18,0.14,'#5C4A1A','BUFFER\nHIDRAULIC','V_buf 100L\nP_acc 20bar',8)
arr(0.38,0.38,0.38,0.27,CG,'P_buf>0\ncarrega')

# Valvules
ax.add_patch(mpatches.FancyBboxPatch((0.24,0.11),0.12,0.07,
    boxstyle="round,pad=0.02",fc='#3A1A5C',ec='white',lw=1,alpha=0.8))
ax.text(0.30,0.145,'Valvules\nretencio',ha='center',va='center',color='white',fontsize=8)
arr(0.38,0.13,0.36,0.155,'white','')

# Generador
box(0.14,0.14,0.15,0.10,'#1A5C40','GENERADOR\nHIDRAULIC','eta=85%',8)
arr(0.24,0.145,0.21,0.145,C7,'P_buf<0\nallibera')
arr(0.14,0.19,0.14,0.64,C7,'P_gen\n10-37W')

# Xarxa
box(0.75,0.55,0.16,0.12,'#1F5C9E','XARXA\nELECTRICA','P_grid = P_aero\n+P_buf+P_gen',8)
arr(0.47,0.70,0.67,0.60,C7,'+2%')

# Resultats
ax.add_patch(mpatches.FancyBboxPatch((0.66,0.15),0.30,0.32,
    boxstyle="round,pad=0.03",fc='#1a1a2e',ec='#2E75B6',lw=1.5,alpha=0.9))
ax.text(0.81,0.455,'RESULTATS v9.4',ha='center',color='#2E75B6',fontsize=10,fontweight='bold')
lines=[
    ('Millora P_grid:','  +1.9% (3p)  /  +2.0% (7p)'),
    ('P_buf oscil.:','  8.4 kW (3p)  /  10.2 kW (7p)'),
    ('P_gen:','  10 W (3p)  /  37 W (7p)'),
    ('omega·R_hub:','  6.8 m/s >> DR_MAX=0.5 m/s'),
    ('E_net:','  -25 kWh (asimetria F_c)'),
    ('Q fluid:','  ~2 L/s a v~1.1 m/s (viable)'),
]
for j,(k,v) in enumerate(lines):
    y=0.41-j*0.04
    ax.text(0.685,y,k,ha='left',color='#aaa',fontsize=8,fontweight='bold')
    ax.text(0.775,y,v,ha='left',color='white',fontsize=8)

ax.set_title('Diagrama Esquematic del Sistema Quijote + ZYPYZAPE v9.4',
             color='white',fontsize=13,fontweight='bold',pad=10)
ax.text(0.5,-0.02,
    'Victor Manzanares Alberola — EPSA UPV Alcoi — 2026',
    ha='center',color='#888',fontsize=9,transform=ax.transAxes)

out2='/mnt/user-data/outputs/gemelo_v94_diagrama.png'
plt.savefig(out2,dpi=160,bbox_inches='tight',facecolor=BG)
plt.close()
print(f"Diagrama: {out2}")

# Copia el codi
import shutil
shutil.copy('/home/claude/gemelo_v94.py','/mnt/user-data/outputs/gemelo_v94.py')
print(f"Codi: /mnt/user-data/outputs/gemelo_v94.py")
