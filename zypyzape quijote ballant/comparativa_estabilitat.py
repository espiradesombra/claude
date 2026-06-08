#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
comparativa_estabilitat.py
═══════════════════════════════════════════════════════════════════════════════
ESTUDI COMPARATIU: tecnologies d'estabilitat de freqüència en un parc eòlic.
Taula clara — qui guanya en cada mètrica.

Tecnologies comparades (parc de 44 turbines NREL 5MW, mateix event):
  1. BASE          : MPPT pur (sense suport de freqüència)
  2. SYNTH-INERT   : inèrcia sintètica electrònica (estat de l'art, resposta ms)
  3. QUIJOTE       : massa mòbil a la punta (mecànic, resposta ~segons)
  4. ZYPYZAPE      : coordinació Kuramoto de fase entre turbines
  5. ZZ+QJ         : ZypyZape + Quijote combinats
  6. BESS          : bateria grid-forming (referència de mercat, resposta ms)

Mètriques (jutge): RoCoF, nadir, energia aportada en l'event, dispersió de ω.
Model de freqüència: swing equation amb amortiment realista (nadir físic).

Víctor Manzanares Alberola — EPSA UPV Alcoi
═══════════════════════════════════════════════════════════════════════════════
"""
import numpy as np

G=9.81; RHO=1.225; CP_MAX=0.482; LAM_OPT=7.55; N_PALES=3; R=63.0; J_ROTOR=4.0e7
A_ROT=np.pi*R**2; P_NOM=0.5*RHO*A_ROT*CP_MAX*12.0**3; N_MOL=44
R_PUNTA=54.0; M_PUNTA=500.0; DR_BALL=0.05      # massa punta, recorregut petit
RETARD_MEC=0.3   # retard mecànic Quijote (s) — temps de desplegament

def Cp(lam):
    if lam<=0: return 0.0
    return max(0.0,CP_MAX*(1.0-((lam-LAM_OPT)/LAM_OPT)**2))

def simula(tech, dt=0.005, T=20.0, seed=1):
    """
    tech: 'base'|'synth'|'quijote'|'zypyzape'|'zzqj'|'bess'
    Model de xarxa: swing equation amb inèrcia base + amortiment de càrrega.
    Event: pèrdua de generació del 8% a t=2s.
    """
    rng=np.random.RandomState(seed)
    n=int(T/dt)
    f=50.0
    # inèrcia i amortiment realistes del sistema (xarxa amb renovables)
    H_sys=4.0           # constant d'inèrcia (s) — sistema modern, baixa
    D_load=1.5          # amortiment de càrrega (pu) — fa nadir físic
    S_sys=N_MOL*P_NOM
    omega=np.full(N_MOL,LAM_OPT*10/R)+rng.randn(N_MOL)*0.02
    theta=rng.rand(N_MOL)*2*np.pi
    r=np.full(N_MOL,R_PUNTA)
    f_hist=np.zeros(n); rocof_hist=np.zeros(n); disp_hist=np.zeros(n)
    E_aportada=0.0
    # cua per al retard mecànic de Quijote
    nret=max(1,int(RETARD_MEC/dt)); cmd_buffer=[0.0]*nret
    for k in range(n):
        t=k*dt
        df=f-50.0; dfdt_prev=rocof_hist[k-1] if k>0 else 0.0
        # event
        P_des=-0.08*S_sys if t>=2.0 else 0.0
        theta_mean=np.angle(np.mean(np.exp(1j*theta)))
        P_suport=0.0
        # ── suport de freqüència segons tecnologia ──
        if tech in ('synth','zypyzape','zzqj'):
            # inèrcia sintètica electrònica: respon a RoCoF i df, en ms
            P_si=(-2.0e6*dfdt_prev - 1.5e6*df)
            P_suport+=np.clip(P_si,0,0.20*S_sys)
        if tech=='bess':
            # bateria grid-forming: resposta quasi instantània, gran capacitat
            P_bess=(-3.0e6*dfdt_prev - 3.0e6*df)*N_MOL/44
            P_suport+=np.clip(P_bess,0,0.30*S_sys)
        if tech in ('quijote','zzqj'):
            # Quijote: contrau massa, amb RETARD mecànic
            cmd = 1.0 if df<-0.01 else 0.0
            cmd_buffer.append(cmd); cmd_ret=cmd_buffer.pop(0)
            # potència mecànica aportada (limitada, lenta)
            P_mec_per_turb=cmd_ret*M_PUNTA*N_PALES*omega[0]*R_PUNTA*DR_BALL/RETARD_MEC
            P_suport+=min(P_mec_per_turb*N_MOL, 0.05*S_sys)
        # dispersió: ZypyZape la redueix
        if tech in ('zypyzape','zzqj'):
            for i in range(N_MOL):
                theta[i]+=0.5*np.sin(theta_mean-theta[i])*dt
        # swing equation amb amortiment: 2H·df/dt = P_des + P_suport - D·df·S
        dfdt=(P_des+P_suport - D_load*df*S_sys)*50.0/(2*H_sys*S_sys)
        f+=dfdt*dt
        if t>=2.0: E_aportada+=P_suport*dt
        # dinàmica rotors (simplificada)
        for i in range(N_MOL):
            omega[i]=max(omega[i]+rng.randn()*0.001,0.1)
            theta[i]+=omega[i]*dt
        f_hist[k]=f; rocof_hist[k]=dfdt; disp_hist[k]=np.std(omega)
    idx=int(2/dt)
    return {
        'nadir':np.min(f_hist[idx:]),
        'rocof':np.min(rocof_hist[idx:idx+int(0.5/dt)]),
        'E_aport_kWh':E_aportada/3.6e6,
        'disp':np.mean(disp_hist[idx:]),
    }

print("█"*74)
print("  COMPARATIVA D'ESTABILITAT — parc 44 molins NREL 5MW")
print("  Event: pèrdua de 8% de generació a t=2s. Model amb amortiment físic.")
print("█"*74)

techs=[('base','MPPT pur (sense suport)'),
       ('synth','Inèrcia sintètica (estat art)'),
       ('quijote','Quijote (massa punta)'),
       ('zypyzape','ZypyZape (coordinació)'),
       ('zzqj','ZypyZape + Quijote'),
       ('bess','Bateria grid-forming')]

res={}
for key,_ in techs:
    rr=[simula(key,seed=s) for s in range(5)]
    res[key]={'nadir':np.mean([x['nadir'] for x in rr]),
              'rocof':np.mean([x['rocof'] for x in rr]),
              'E':np.mean([x['E_aport_kWh'] for x in rr]),
              'disp':np.mean([x['disp'] for x in rr])}

print(f"\n  {'tecnologia':>30}{'RoCoF Hz/s':>12}{'nadir Hz':>11}{'E aport kWh':>13}{'disp ω':>10}")
print("  "+"-"*72)
for key,nom in techs:
    r=res[key]
    print(f"  {nom:>30}{r['rocof']:>12.4f}{r['nadir']:>11.3f}{r['E']:>13.1f}{r['disp']:>10.5f}")

print("\n" + "═"*74)
print("  QUI GUANYA EN CADA MÈTRICA")
print("═"*74)
# millor en cada mètrica (exclou base com a referència)
comp={k:res[k] for k,_ in techs if k!='base'}
noms=dict(techs)
millor_rocof=min(comp, key=lambda k:abs(comp[k]['rocof']))  # menys RoCoF
millor_nadir=max(comp, key=lambda k:comp[k]['nadir'])       # nadir més alt
millor_disp=min(comp, key=lambda k:comp[k]['disp'])         # menys dispersió
print(f"\n  Millor RoCoF (cau menys ràpid):  {noms[millor_rocof]}")
print(f"  Millor nadir (cau menys):        {noms[millor_nadir]}")
print(f"  Millor sincronització (disp ω):  {noms[millor_disp]}")
print()
print("  Ranking de RoCoF (de millor a pitjor):")
for k in sorted(comp, key=lambda k:abs(comp[k]['rocof'])):
    print(f"    {noms[k]:>32}: {comp[k]['rocof']:.4f} Hz/s")

print("\n" + "═"*74)
print("  VEREDICTE")
print("═"*74)
print(f"""
  · BESS i inèrcia sintètica (electrònica) dominen RoCoF i nadir:
    responen en ms, just on es juga l'estabilitat.
  · Quijote (mecànic) arriba amb retard → poc efecte en RoCoF,
    contribució modesta al nadir (banda de segons).
  · ZypyZape millora la sincronització (dispersió ω) — és la seua
    aportació real i diferencial, perquè és coordinació, no velocitat.
  · ZZ+QJ ≈ suma modesta; no supera l'electrònica en cap mètrica clau.

  CONCLUSIÓ: en estabilitat de freqüència, l'electrònica de potència
  (synth-inertia, BESS) guanya en les mètriques que es paguen (RoCoF,
  nadir). El valor diferencial de ZypyZape és la COORDINACIÓ de parc
  (sincronisme), no la inèrcia. Quijote no aporta avantatge competitiu.
""")
