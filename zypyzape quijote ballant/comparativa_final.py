#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
comparativa_final.py
═══════════════════════════════════════════════════════════════════════════════
LES 4 COMPARATIVES de Víctor, totes contra MÀXIMA TECNOLOGIA ACTUAL.
Mètriques completes: energia, nadir de freqüència, RoCoF, dispersió de ω.

"Màxima tecnologia actual" (BASE) = turbina moderna amb:
  · control de parell MPPT (τ=Kω²)
  · synthetic inertia electrònica (allibera E_cinètica del rotor via convertidor)
  · resposta en mil·lisegons, sense parts mòbils

Configuracions comparades (44 molins, NREL 5MW):
  BASE   : 44 turbines màxima tecnologia actual
  ZZ     : BASE + ZypyZape (coordinació Kuramoto de fase entre turbines)
  QJ     : BASE + Quijote (massa mòbil a la PUNTA, ball ràpid Δr petit)
  ZZ+QJ  : BASE + ZypyZape + Quijote integrats

Realisme: sensor de vent amb retard+soroll, η actuador, fatiga, pes a punta.

Víctor Manzanares Alberola — EPSA UPV Alcoi
═══════════════════════════════════════════════════════════════════════════════
"""
import numpy as np

G=9.81; RHO=1.225; CP_MAX=0.482; LAM_OPT=7.55; N_PALES=3; R=63.0
J_ROTOR=4.0e7; A_ROT=np.pi*R**2
P_NOM=0.5*RHO*A_ROT*CP_MAX*12.0**3
N_MOL=44
# Quijote: pes a la PUNTA, ball curt ràpid (la idea de Víctor)
R_PUNTA=54.0; DR_BALL=0.5      # mou només 0.5m prop de la punta
M_PUNTA=1000.0                 # massa gran a la punta
DR_MAX=1.0; DDR_MAX=2.0        # a la punta pot ser ràpid (ball curt)
ETA_ACT=0.85; C_FRIC=200.0; K_FATIGA=0.15
RETARD_S=1.0; SOROLL_V=0.4; TAU_FILT=2.0

def Cp(lam):
    if lam<=0: return 0.0
    return max(0.0,CP_MAX*(1.0-((lam-LAM_OPT)/LAM_OPT)**2))

def simula(config, dt=0.01, T=15.0, seed=1):
    """
    config: 'base' | 'zz' | 'qj' | 'zzqj'
    Retorna mètriques: E_total, nadir, rocof_pic, dispersio_omega
    """
    rng=np.random.RandomState(seed)
    n=int(T/dt)
    omega=np.full(N_MOL,LAM_OPT*10/R)+rng.randn(N_MOL)*0.02
    theta=rng.rand(N_MOL)*2*np.pi
    r=np.full(N_MOL,R_PUNTA); rdot=np.zeros(N_MOL)
    v_base=10.0+rng.randn(N_MOL)*0.5
    v_filt=v_base.copy(); buf=[v_base.copy() for _ in range(max(1,int(RETARD_S/dt)))]
    K_mppt=0.5*RHO*A_ROT*R**3*CP_MAX/LAM_OPT**3
    f_grid=50.0; H_sys=3.0; S_sys=N_MOL*P_NOM
    K_ZZ=0.5 if config in ('zz','zzqj') else 0.0
    has_QJ = config in ('qj','zzqj')
    f_hist=np.zeros(n); rocof_hist=np.zeros(n); disp_hist=np.zeros(n)
    E_total=0.0; E_act_total=0.0; W_fat_total=0.0
    rng2=np.random.RandomState(7)
    for k in range(n):
        t=k*dt
        v=v_base+0.5*np.sin(2*np.pi*t/8)
        # sensor realista
        vr=v+SOROLL_V*rng2.randn(N_MOL); buf.append(vr); vd=buf.pop(0)
        v_filt+=(vd-v_filt)*dt/TAU_FILT
        # event de freqüència a t=2s
        df=-0.6*np.exp(-(t-2)/5)*(1-np.exp(-(t-2)/0.3)) if t>=2 else 0.0
        theta_mean=np.angle(np.mean(np.exp(1j*theta)))
        P_inject_total=0.0
        for i in range(N_MOL):
            lam=omega[i]*R/max(v[i],0.1)
            P_vent=min(0.5*RHO*A_ROT*Cp(lam)*v[i]**3,P_NOM*1.3)
            tau_vent=P_vent/max(omega[i],0.1)
            tau_gen_base=min(K_mppt*omega[i]**2,P_NOM/max(omega[i],0.1))
            # synthetic inertia electrònica (TOTES, és l'estat de l'art)
            P_si=np.clip(-0.4e6*df,0,0.15*P_NOM)
            tau_gen=tau_gen_base+P_si/max(omega[i],0.1)
            # Quijote: mou massa a la punta si event
            rr=0.0
            if has_QJ:
                if df<-0.02 and r[i]>R_PUNTA-DR_BALL: rc=-DR_MAX
                elif df>-0.005 and r[i]<R_PUNTA: rc=DR_MAX
                else: rc=0.0
                # ZypyZape modula la fase del ball
                if config=='zzqj':
                    rc=rc*(1+0.3*np.cos(theta[i]-theta_mean))
                dv=np.clip(rc-rdot[i],-DDR_MAX*dt,DDR_MAX*dt); rdot[i]+=dv
                r_new=np.clip(r[i]+rdot[i]*dt,R_PUNTA-DR_BALL,R_PUNTA); rr=(r_new-r[i])/dt
                r[i]=r_new
            massa=M_PUNTA if has_QJ else 0.0
            Jdot=N_PALES*2*massa*r[i]*rr if has_QJ else 0.0
            P_quij=-omega[i]*Jdot if has_QJ else 0.0
            J=J_ROTOR+N_PALES*massa*r[i]**2
            # ZypyZape: acoblament de fase
            kura=K_ZZ*np.sin(theta_mean-theta[i])
            P_inject_total += (tau_gen-tau_gen_base)*omega[i] + P_quij
            E_total += tau_gen*omega[i]*dt
            if has_QJ:
                Fc=massa*omega[i]**2*r[i]; Fg=-massa*G*np.sin(theta[i])
                Fa=-(Fc+Fg)-C_FRIC*rr; P=Fa*rr
                E_act_total+=(P/ETA_ACT if P>0 else P*ETA_ACT)*dt
                W_fat_total+=K_FATIGA*abs(Fc*rr)*dt
            omega[i]=max(omega[i]+(tau_vent-tau_gen-Jdot*omega[i])/J*dt,0.1)
            theta[i]+=(omega[i]+kura)*dt
        P_des=-0.05*S_sys if t>=2 else 0.0
        dfdt=(P_des+P_inject_total)*50.0/(2*H_sys*S_sys)
        f_grid+=dfdt*dt
        f_hist[k]=f_grid; rocof_hist[k]=dfdt; disp_hist[k]=np.std(omega)
    idx=int(2/dt)
    E_net=(E_total-max(E_act_total,0)-W_fat_total)/3.6e6  # kWh
    return {
        'E_net_kWh':E_net,
        'nadir':np.min(f_hist[idx:]),
        'rocof_pic':np.min(rocof_hist[idx:idx+int(1/dt)]),
        'disp_omega':np.mean(disp_hist[idx:]),
    }

print("█"*72)
print("  LES 4 COMPARATIVES — totes contra MÀXIMA TECNOLOGIA ACTUAL")
print("█"*72)
print("\n  44 molins NREL 5MW. Event de freqüència a t=2s. Pes Quijote a la")
print("  PUNTA (ball curt ràpid). Sensor realista + fatiga + η actuador.\n")

# mitjana sobre diversos seeds per robustesa
configs=['base','zz','qj','zzqj']
noms={'base':'Màx. tec. actual (BASE)','zz':'BASE + ZypyZape',
      'qj':'BASE + Quijote (punta)','zzqj':'BASE + ZypyZape + Quijote'}
res={c:{'E':[],'nadir':[],'rocof':[],'disp':[]} for c in configs}
for seed in range(5):
    for c in configs:
        r=simula(c,seed=seed)
        res[c]['E'].append(r['E_net_kWh']); res[c]['nadir'].append(r['nadir'])
        res[c]['rocof'].append(r['rocof_pic']); res[c]['disp'].append(r['disp_omega'])

print(f"  {'configuració':>28}{'E_net kWh':>12}{'nadir Hz':>11}{'RoCoF Hz/s':>12}{'disp ω':>10}")
print("  "+"-"*70)
base=res['base']
for c in configs:
    E=np.mean(res[c]['E']); nad=np.mean(res[c]['nadir'])
    roc=np.mean(res[c]['rocof']); dis=np.mean(res[c]['disp'])
    print(f"  {noms[c]:>28}{E:>12.1f}{nad:>11.4f}{roc:>12.4f}{dis:>10.5f}")
print()

# comparativa relativa vs BASE
print("="*72)
print("  DIFERÈNCIA vs MÀXIMA TECNOLOGIA ACTUAL (BASE)")
print("="*72)
bE=np.mean(base['E']); bN=np.mean(base['nadir']); bR=np.mean(base['rocof']); bD=np.mean(base['disp'])
print(f"\n  {'config':>28}{'ΔE %':>10}{'Δnadir mHz':>13}{'ΔRoCoF %':>11}{'Δdisp %':>10}")
for c in ['zz','qj','zzqj']:
    E=np.mean(res[c]['E']); nad=np.mean(res[c]['nadir'])
    roc=np.mean(res[c]['rocof']); dis=np.mean(res[c]['disp'])
    dE=100*(E-bE)/abs(bE); dN=1000*(nad-bN)
    dR=100*(abs(roc)-abs(bR))/abs(bR); dD=100*(dis-bD)/bD
    print(f"  {noms[c]:>28}{dE:>9.2f}%{dN:>13.1f}{dR:>10.1f}%{dD:>9.1f}%")
print()
print("  Llegenda: Δnadir positiu=millor (cau menys). ΔRoCoF negatiu=millor.")
print("            Δdisp negatiu=millor (parc més sincronitzat).")

print("\n" + "═"*72)
print("  VEREDICTE DE LES 4 COMPARATIVES")
print("═"*72)
# calcular per a conclusió
zz_dN=1000*(np.mean(res['zz']['nadir'])-bN)
qj_dN=1000*(np.mean(res['qj']['nadir'])-bN)
zzqj_dN=1000*(np.mean(res['zzqj']['nadir'])-bN)
zz_dD=100*(np.mean(res['zz']['disp'])-bD)/bD
print(f"""
  1. ZypyZape vs actual:  nadir {zz_dN:+.1f} mHz, dispersió {zz_dD:+.1f}%
     → la coordinació {'ajuda' if zz_dN>1 or zz_dD<-1 else 'aporta poc'} sobre l'estat de l'art.

  2. Quijote vs actual:   nadir {qj_dN:+.1f} mHz
     → la massa a la punta {'aporta' if qj_dN>1 else 'aporta poc'} sobre la inèrcia electrònica
       que la turbina moderna JA té.

  3. ZypyZape+Quijote vs actual: nadir {zzqj_dN:+.1f} mHz
     → la combinació {'supera' if zzqj_dN>max(zz_dN,qj_dN)+0.5 else 'no supera clarament'} la suma de les parts.

  CLAU: la "màxima tecnologia actual" ja inclou inèrcia sintètica
  electrònica (resposta en ms). Quijote (mecànic) competeix contra
  això en desavantatge de velocitat. ZypyZape (coordinació) és la
  part que pot aportar valor real, perquè és software sobre el control.
""")
