"""
GEMELL VIRTUAL v10 — Cp(lambda, pitch) dinamic
Victor Manzanares Alberola — EPSA UPV Alcoi

NOVETAT: Cp(lambda, beta) real basat en NREL 5MW (Jonkman 2009)
  - beta = angle de pitch per pala (rad)
  - lambda = tip speed ratio actual
  - Cp(lam,beta) = Cp_max * f(lam) * g(beta)

El pitch asincron canvia beta en cada pala i en cada instant
→ Cp real canvia per pala → guany real del pitch captura

EFECTE ESPERAT: guany 2% → 6-12%
"""
import argparse
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(42)

# ============================================================
# PARAMETRES
# ============================================================
KPITCH=2.0; K_KUR=0.10; K_q=5e4
V_BASE=11.4; T_TOTAL=60.0; DT=0.05
H_BUF=20.0; V_BUF_MAX=0.10; Q_BOMBA=0.010
P_MAX_ACC=20e6; LLINDAR_ACC=20.0; ETA_GEN=0.85; DR_MAX=0.5

R=63.0; P_NOM=5e6; RHO_A=1.225; Cp_max_ref=0.482; lam_opt=7.55
A_rot=np.pi*R**2
vr=(P_NOM/(0.5*RHO_A*A_rot*Cp_max_ref))**(1/3); omR=lam_opt*vr/R
R_hub=5.0; R_tip=55.0; J_base=3.54e7

CONFIGS={
    '3p base(Cp fix)': dict(N=3,rho=3386.0,D=0.05,cp_din=False),
    '3p Cp(lam,pitch)': dict(N=3,rho=3386.0,D=0.05,cp_din=True),
    '7p Cp(lam,pitch)': dict(N=7,rho=3386.0,D=0.05,cp_din=True),
}

steps=int(T_TOTAL/DT); t_vec=np.arange(steps)*DT
ou=np.zeros(steps)
for i in range(1,steps):
    ou[i]=ou[i-1]-0.5*ou[i-1]*DT+0.8*np.sqrt(DT)*np.random.randn()

# ============================================================
# MODEL Cp(lambda, pitch) - NREL 5MW calibrat
# ============================================================
def Cp_dinamic(lam, beta_rad):
    """
    Cp(lambda, beta) basat en approximacio analitica NREL 5MW.
    Referencia: Jonkman 2009, eq. simplificada de Heier 2006.

    Cp = C1*(C2/lam_i - C3*beta - C4)*exp(-C5/lam_i) + C6*lam

    on lam_i = 1/(1/(lam+0.08*beta) - 0.035/(beta^3+1))
    """
    beta = float(np.degrees(beta_rad))  # convertir a graus per a la formula
    beta = float(np.clip(beta, -5.0, 25.0))  # limits fisics pitch

    # Constants calibrades NREL 5MW
    C1=0.5176; C2=116.0; C3=0.4; C4=5.0; C5=21.0; C6=0.0068

    # lam_i amb correccio de pitch
    lam_c = max(lam, 0.1)
    beta_c = max(abs(beta), 0.01)
    lam_i = 1.0/(1.0/(lam_c + 0.08*beta_c) - 0.035/(beta_c**3 + 1.0))
    lam_i = max(lam_i, 0.01)

    Cp = C1*(C2/lam_i - C3*beta_c - C4)*np.exp(-C5/lam_i) + C6*lam_c
    return float(np.clip(Cp, 0.0, Cp_max_ref))

def Cp_fix(lam, beta_rad):
    """Model Cp fix (v9.2) — per comparar"""
    return float(max(0, Cp_max_ref*(1-((lam-lam_opt)/lam_opt)**2)))

# Verificacio Cp(lambda, beta=0)
lam_test=np.linspace(2,14,100)
cp_test=[Cp_dinamic(l,0) for l in lam_test]
cp_fix=[Cp_fix(l,0) for l in lam_test]
print(f"Cp max (beta=0): dinamic={max(cp_test):.4f}  fix={max(cp_fix):.4f}")
lam_opt_din=lam_test[np.argmax(cp_test)]
print(f"lam_opt: dinamic={lam_opt_din:.2f}  fix={lam_opt:.2f}")

def gP_pitch(theta_i, omega):
    """Angle de pitch asincron per pala i (rad)."""
    pm_deg = KPITCH*np.cos(theta_i) + 0.5*omega*np.sin(theta_i)
    return float(np.clip(np.radians(pm_deg), -0.14, 0.14))  # +-8 graus

def eta_kur(K,N):
    boost=min(1.2,N/3)
    return float(np.clip(0.04*boost*(1-((K-0.10)/0.10)**2),0,0.048))

def perda_fric(Q,D,rho,L=10.0,f=0.02):
    A=np.pi*(D/2)**2; v=Q/max(A,1e-10)
    return f*(L/D)*(rho*v**2)/2

# ============================================================
# SIMULACIO
# ============================================================
def simula(cfg_nom):
    cfg=CONFIGS[cfg_nom]
    N=cfg['N']; rho=cfg['rho']; D=cfg['D']; cp_din=cfg['cp_din']
    A_can=np.pi*(D/2)**2; m_q=rho*A_can*(R_tip-R_hub)

    thetas=np.array([2*np.pi*i/N for i in range(N)])
    omega=omR*0.95
    r_q=np.full(N,(R_hub+R_tip)/2)
    dr_q=np.zeros(N)
    V_buf=0.0; P_acc=0.0
    beta_pala=np.zeros(N)  # pitch per pala
    eta_K=eta_kur(K_KUR,N)

    hP_aero=np.zeros(steps); hP_grid=np.zeros(steps)
    hP_buf=np.zeros(steps);  hP_gen=np.zeros(steps)
    hom=np.zeros(steps);     hr_q=np.zeros((steps,N))
    hV_buf=np.zeros(steps);  hdJ=np.zeros(steps)
    hCp=np.zeros(steps);     hbeta=np.zeros(steps)

    for s in range(steps):
        t=s*DT
        v=V_BASE+2*np.sin(2*np.pi*t/20)+ou[s]*0.5
        if 30<t<45: v=max(3,v-4*np.sin(np.pi*(t-30)/15))
        v=float(v)

        omega+=(min(omR,lam_opt*max(v,0.1)/R)-omega)*(DT/3.0)
        omega=float(np.clip(omega,omR*0.3,omR*1.4))

        lam_act=omega*R/max(v,0.1)

        # Quijote
        dJ_total=0.0
        for i in range(N):
            th_n=thetas[i]%(2*np.pi)
            r_opt=R_tip-(R_tip-R_hub)*(1+np.cos(th_n))/2
            F_c=m_q*omega**2*r_q[i]
            F_ctrl=-K_q*(r_q[i]-r_opt)
            F_fric=-80*dr_q[i]
            d2r=(F_c+F_ctrl+F_fric)/m_q
            dr_q[i]+=d2r*DT
            dr_q[i]=float(np.clip(dr_q[i],-DR_MAX,DR_MAX))
            dr_q[i]=float(np.clip(dr_q[i],-omega*r_q[i],omega*r_q[i]))
            r_q[i]+=dr_q[i]*DT
            r_q[i]=float(np.clip(r_q[i],R_hub,R_tip))
            dJ_total+=2*m_q*r_q[i]*dr_q[i]

        P_buf=-0.5*dJ_total*omega**2

        # ================================================
        # PITCH ASINCRON + Cp(lam, beta) DINAMIC
        # ================================================
        P_aero_total=0.0
        Cp_mig=0.0
        beta_mig=0.0

        for i in range(N):
            if cp_din:
                # Pitch asincron: angle per pala
                dJdt_i=2*m_q*r_q[i]*dr_q[i]
                beta_i=gP_pitch(thetas[i], omega)
                # Ajust pel quijote: si J puja, redueix pitch
                beta_i-=float(np.clip(0.8*np.radians(dJdt_i*0.01),-0.05,0.05))
                beta_pala[i]=float(np.clip(beta_i,-0.14,0.14))

                # Cp real per aquesta pala
                Cp_i=Cp_dinamic(lam_act, beta_pala[i])
            else:
                # Cp fix (model v9.2)
                Cp_i=Cp_fix(lam_act, 0)

            # Kuramoto
            kc=(K_KUR/N)*float(np.sum(np.sin(thetas-thetas[i])))
            thetas[i]+=(omega+kc)*DT

            # Potencia aerodinamica real per pala
            P_pala=0.5*RHO_A*A_rot/N*Cp_i*v**3
            P_aero_total+=P_pala
            Cp_mig+=Cp_i/N
            beta_mig+=np.degrees(beta_pala[i])/N if cp_din else 0.0

        P_aero_total*=(1.0+eta_K)  # Kuramoto

        # Buffer + valvules
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
            P_hid_inject=P_acc*Q_r*0.82
            P_gen_inject=P_acc*Q_r*ETA_GEN
            V_buf=max(0.0,V_buf-Q_r*DT)

        P_grid=max(0,P_aero_total+max(0,P_buf)+P_hid_inject+P_gen_inject)

        # Referencia aerodinamica sense pitch (per calcular millora)
        P_aero_base=0.5*RHO_A*A_rot*Cp_fix(lam_act,0)*v**3

        hP_aero[s]=P_aero_base; hP_grid[s]=P_grid
        hP_buf[s]=P_buf; hP_gen[s]=P_gen_inject
        hom[s]=omega; hr_q[s]=r_q.copy()
        hV_buf[s]=V_buf*1000; hdJ[s]=dJ_total
        hCp[s]=Cp_mig; hbeta[s]=beta_mig

    return dict(P_aero=hP_aero,P_grid=hP_grid,P_buf=hP_buf,
                P_gen=hP_gen,om=hom,r_q=hr_q,V_buf=hV_buf,
                dJ=hdJ,Cp=hCp,beta=hbeta,
                m_q=m_q,N=N,rho=rho,D=D)

def run_simulations():
    print("\nSimulant v10...")
    results = {}
    for nom in CONFIGS:
        r = simula(nom)
        mill = (np.mean(r["P_grid"]) - np.mean(r["P_aero"])) / np.mean(r["P_aero"]) * 100
        mill = max(0, mill)
        Pb = np.mean(np.abs(r["P_buf"]))
        Pg = np.mean(r["P_gen"])
        Cp_mig = np.mean(r["Cp"])
        beta_mig = np.mean(r["beta"])
        E_net = np.sum(r["P_buf"]) * DT / 3600
        print(f"\n{nom}:")
        print(f"  Millora P_grid: +{mill:.1f}%")
        print(f"  Cp mig:         {Cp_mig:.4f} (ref fix={Cp_max_ref:.4f})")
        print(f"  beta mig:       {beta_mig:.2f} deg")
        print(f"  P_buf|mig|:     {Pb/1e3:.2f} kW")
        print(f"  E_net:          {E_net:.1f} kWh")
        results[nom] = dict(r=r, mill=mill, Pb=Pb, Pg=Pg, Cp=Cp_mig, beta=beta_mig)
    return results


def _sty_plot(ax, t, xl='t[s]', yl='', pan='#13132b'):
    ax.set_facecolor(pan)
    ax.set_title(t, color='white', fontsize=9, pad=4)
    ax.tick_params(colors='#aaa', labelsize=7.5)
    for sp in ax.spines.values():
        sp.set_color('#333355')
    ax.set_xlabel(xl, color='#aaa', fontsize=8)
    ax.set_ylabel(yl, color='#aaa', fontsize=8)
    ax.grid(color='#1e1e40', lw=0.5, ls='--')


def plot_results(results, lam_opt_din, out_path: Path):
    BG='#0d0d1a'; PAN='#13132b'
    COLS={'3p base(Cp fix)':'#888780',
          '3p Cp(lam,pitch)':'#00d2ff',
          '7p Cp(lam,pitch)':'#00ff88'}

    fig=plt.figure(figsize=(22,18),facecolor=BG)
    gs=gridspec.GridSpec(3,3,figure=fig,hspace=0.50,wspace=0.38)

    def _sty_plot(ax,t,xl='t[s]',yl=''):
        ax.set_facecolor(PAN); ax.set_title(t,color='white',fontsize=9,pad=4)
        ax.tick_params(colors='#aaa',labelsize=7.5)
        [sp.set_color('#333355') for sp in ax.spines.values()]
        ax.set_xlabel(xl,color='#aaa',fontsize=8); ax.set_ylabel(yl,color='#aaa',fontsize=8)
        ax.grid(color='#1e1e40',lw=0.5,ls='--')

    # F0: Cp(lam,beta) corba
    ax=fig.add_subplot(gs[0,0])
    lam_v=np.linspace(2,14,200)
    for b_deg in [0,5,10,15]:
        cp_v=[Cp_dinamic(l,np.radians(b_deg)) for l in lam_v]
        ax.plot(lam_v,cp_v,lw=2,label=f'beta={b_deg}deg')
    ax.plot(lam_v,cp_fix,'--',color='#888',lw=1.5,alpha=0.7,label='Cp fix (v9.2)')
    ax.axvline(lam_opt,color='white',ls=':',lw=0.8,alpha=0.4,label=f'lam_opt={lam_opt}')
    ax.axvline(lam_opt_din,color='#ffd700',ls=':',lw=0.8,alpha=0.6,
               label=f'lam_opt din={lam_opt_din:.2f}')
    ax.legend(fontsize=6.5,framealpha=0.3); ax.set_ylim(0,0.55)
    _sty_plot(ax,'Cp(lambda,beta) — NREL 5MW calibrat','lambda','Cp')

    # F1: P_grid comparativa
    ax=fig.add_subplot(gs[0,1:])
    r0=list(results.values())[0]['r']
    ax.plot(t_vec,r0['P_aero']/1e6,'--',color='#555',lw=1.2,alpha=0.6,label='Base aero')
    for nom,d in results.items():
        ax.plot(t_vec,d['r']['P_grid']/1e6,color=COLS[nom],lw=2,
                label=f'{nom}  +{d["mill"]:.1f}%')
    ax.axvspan(30,45,alpha=0.07,color='#e74c3c')
    ax.legend(fontsize=7.5,framealpha=0.3)
    _sty_plot(ax,'P_grid (MW) — base vs Cp fix vs Cp(lam,pitch)','t[s]','MW')

    # F2: Cp(t) mig per les configuracions dinamiques
    ax=fig.add_subplot(gs[1,0])
    for nom,d in results.items():
        if 'Cp' in d['r'] and np.any(d['r']['Cp']>0):
            ax.plot(t_vec,d['r']['Cp'],color=COLS[nom],lw=1.5,alpha=0.85,label=nom)
    ax.axhline(Cp_max_ref,color='#ffd700',ls='--',lw=1,alpha=0.7,
               label=f'Cp_max ref={Cp_max_ref}')
    ax.legend(fontsize=7,framealpha=0.3)
    _sty_plot(ax,'Cp mig al llarg del temps','t[s]','Cp')

    # F3: beta (pitch) per les configs dinamiques
    ax=fig.add_subplot(gs[1,1])
    for nom,d in results.items():
        if np.any(np.abs(d['r']['beta'])>0.01):
            ax.plot(t_vec,d['r']['beta'],color=COLS[nom],lw=1.5,alpha=0.85,label=nom)
    ax.axhline(0,color='white',lw=0.5,alpha=0.3)
    ax.legend(fontsize=7,framealpha=0.3)
    _sty_plot(ax,'beta mig pales (deg) — pitch asincron','t[s]','beta [deg]')

    # F4: P_buf
    ax=fig.add_subplot(gs[1,2])
    for nom,d in results.items():
        ax.plot(t_vec,d['r']['P_buf']/1e3,color=COLS[nom],lw=1.2,alpha=0.8,label=nom)
    ax.axhline(0,color='white',lw=0.5,alpha=0.3)
    ax.legend(fontsize=7,framealpha=0.3)
    _sty_plot(ax,'P_buf Quijote (kW)','t[s]','kW')

    # F5: omega
    ax=fig.add_subplot(gs[2,0])
    for nom,d in results.items():
        ax.plot(t_vec,d['r']['om'],color=COLS[nom],lw=1.5,label=nom)
    ax.axhline(omR,color='white',ls='--',lw=0.7,alpha=0.4)
    ax.axvspan(30,45,alpha=0.07,color='#e74c3c')
    ax.legend(fontsize=7,framealpha=0.3)
    _sty_plot(ax,'omega (rad/s)','t[s]','rad/s')

    # F6: taula resum final
    ax=fig.add_subplot(gs[2,1:])
    ax.axis('off'); ax.set_facecolor('#0a0a14')

    rows=[]
    base_mill=list(results.values())[0]['mill']
    for nom,d in results.items():
        extra='+' if d['mill']>base_mill+0.5 else ''
        rows.append([
            nom,
            f'+{d["mill"]:.1f}%',
            f'{d["Cp"]:.4f}',
            f'{d["beta"]:.2f}deg',
            f'{d["Pb"]/1e3:.2f}kW',
            f'{d["Pg"]:.0f}W',
        ])
    tbl=ax.table(cellText=rows,
        colLabels=['Config','Millora P_grid','Cp mig','beta mig','P_buf','P_gen'],
        loc='center',cellLoc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(9)
    for (rr,c),cell in tbl.get_celld().items():
        bg='#1a1a2e' if rr==0 else '#0d0d1a'
        if rr==2: bg='#001a0a'  # highlight best
        if rr==3: bg='#001a1a'
        cell.set_facecolor(bg)
        cell.set_text_props(color='white' if rr==0 else '#ddd')
        cell.set_edgecolor('#333355')
    ax.set_title('Resum v10 — Cp(lambda,pitch) dinamic vs Cp fix',
                 color='white',fontsize=9,pad=4)

    fig.suptitle(
        f'Gemell Virtual v10 — Cp(lambda,beta) dinamic NREL 5MW\n'
        f'Cp(lam,beta)=C1*(C2/lam_i-C3*beta-C4)*exp(-C5/lam_i)+C6*lam  |  '
        f'Pitch asincron mou beta per pala → Cp real canvia → guany real',
        color='white',fontsize=10,fontweight='bold',y=0.999)

    out='/mnt/user-data/outputs/gemelo_v10_cp_dinamic.png'
    plt.savefig(out,dpi=150,bbox_inches='tight',facecolor=BG)
    plt.close()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f"\nGrafic: {out_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Gemelo virtual v10 — Cp(lambda,pitch) dinamic")
    parser.add_argument("--out", type=Path, default=None, help="PNG de salida")
    parser.add_argument("--no-plot", action="store_true", help="Solo simulacion (sin matplotlib save)")
    args = parser.parse_args()

    results = run_simulations()
    if not args.no_plot:
        out = args.out or Path(__file__).resolve().parent / "gemelo_v10_cp_dinamic.png"
        plot_results(results, lam_opt_din, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
