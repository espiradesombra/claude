#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
comparativa_rocof.py
═══════════════════════════════════════════════════════════════════════════════
COMPARATIVA DE RoCoF — una sola cosa, ben feta.
Atén les correccions del revisor (GPT):
  · RoCoF és la conclusió MÉS sòlida (separació d'escales temporals).
  · EVITA doble comptabilització: H_sys NO inclou els rotors del parc;
    els rotors són fonts que injecten potència de suport per separat.
  · No barreja nadir (dependent de paràmetres) — només RoCoF inicial.

Compara el RoCoF (taxa de canvi de freqüència en els primers 500 ms) després
d'una pèrdua sobtada de generació, segons quin mecanisme de suport actua:

  BASE   : sense suport ràpid
  SYNTH  : inèrcia sintètica electrònica (τ ≈ 20 ms)
  BESS   : bateria grid-forming (τ ≈ 10 ms)
  QUIJOTE: massa mòbil (τ ≈ 300 ms, retard mecànic de desplegament)

El RoCoF inicial es mesura en la finestra 0-500 ms, on es juga l'estabilitat.

Víctor Manzanares Alberola — EPSA UPV Alcoi
═══════════════════════════════════════════════════════════════════════════════
"""
import numpy as np

# ── Paràmetres de xarxa (SENSE comptar rotors dins de H_sys) ──────────────────
F_NOM = 50.0
S_SYS = 1000e6           # potència base del sistema (1 GVA)
# H_sys = inèrcia de la xarxa EXCLOENT els 44 rotors (síncrons clàssics + càrrega)
# Açò evita la doble comptabilització que advertia el revisor.
H_SYS = 3.0             # s (sistema modern amb poca inèrcia síncrona)
D_LOAD = 1.0            # amortiment de càrrega (pu)
DELTA_P = 0.10          # pèrdua sobtada: 10% de la generació

# ── Mecanismes de suport: cada un té una constant de temps τ ──────────────────
# Modelats com a resposta de primer ordre: P_suport segueix la demanda amb τ
MECANISMES = {
    'BASE'   : {'tau': None,  'gain': 0.0,  'nom':'Sense suport ràpid'},
    'SYNTH'  : {'tau': 0.020, 'gain': 0.8,  'nom':'Inèrcia sintètica (electrònica)'},
    'BESS'   : {'tau': 0.010, 'gain': 1.0,  'nom':'Bateria grid-forming'},
    'QUIJOTE': {'tau': 0.300, 'gain': 0.8,  'nom':'Quijote (massa mòbil)'},
}

def simula_rocof(mec, dt=0.001, T=2.0):
    """
    Swing equation: 2·H_sys·S/F_nom · df/dt = ΔP_event + P_suport − D·Δf·S
    P_suport = mecanisme de suport amb la seua constant de temps τ.
    Mesura el RoCoF (df/dt) en la finestra 0-500 ms.
    """
    n=int(T/dt)
    f=F_NOM; df=0.0
    P_sup=0.0
    cfg=MECANISMES[mec]
    rocof_hist=np.zeros(n); f_hist=np.zeros(n); t_hist=np.zeros(n)
    for k in range(n):
        t=k*dt
        # event a t=0.1s
        P_event = -DELTA_P*S_SYS if t>=0.1 else 0.0
        # demanda de suport: proporcional a la desviació i al RoCoF
        df_grid = f - F_NOM
        if cfg['tau'] is not None:
            P_objectiu = cfg['gain']*(-df_grid)*S_SYS*2.0   # suport proporcional
            P_objectiu = np.clip(P_objectiu, 0, 0.15*S_SYS)
            # resposta de primer ordre amb constant τ (ací entra l'escala temporal!)
            P_sup += (P_objectiu - P_sup)*dt/cfg['tau']
        else:
            P_sup = 0.0
        # swing equation (H_sys NO inclou rotors → sense doble comptabilitat)
        dfdt = (P_event + P_sup - D_LOAD*df_grid*S_SYS)*F_NOM/(2*H_SYS*S_SYS)
        f += dfdt*dt
        rocof_hist[k]=dfdt; f_hist[k]=f; t_hist[k]=t
    # RoCoF: pitjor valor en finestra 0-500ms després de l'event
    idx0=int(0.1/dt); idx1=int(0.6/dt)
    rocof_pitjor = np.min(rocof_hist[idx0:idx1])
    # RoCoF mitjà en els primers 500ms (mètrica estàndard de xarxa)
    rocof_500 = (f_hist[idx1]-f_hist[idx0])/0.5
    nadir = np.min(f_hist[idx0:])
    return rocof_pitjor, rocof_500, nadir

print("█"*72)
print("  COMPARATIVA DE RoCoF — quin mecanisme frena millor la caiguda inicial")
print("█"*72)
print(f"\n  Xarxa: H_sys={H_SYS}s (SENSE rotors, evita doble comptabilitat),")
print(f"  event = pèrdua sobtada del {DELTA_P*100:.0f}% de generació.")
print(f"  Mètrica: RoCoF en la finestra 0-500 ms (on es juga l'estabilitat).\n")

print(f"  {'mecanisme':>34}{'τ resposta':>12}{'RoCoF pic':>12}{'RoCoF 500ms':>14}")
print("  "+"-"*70)
resultats={}
for mec in ['BASE','SYNTH','BESS','QUIJOTE']:
    rp,r500,nad=simula_rocof(mec)
    resultats[mec]=(rp,r500)
    tau=MECANISMES[mec]['tau']
    tau_str=f"{tau*1000:.0f} ms" if tau else "—"
    print(f"  {MECANISMES[mec]['nom']:>34}{tau_str:>12}{rp:>11.3f} {r500:>13.3f}")

print("\n" + "═"*72)
print("  QUI GUANYA EN RoCoF")
print("═"*72)
# ordenar per RoCoF a 500ms (menys negatiu = millor)
ordre=sorted(['SYNTH','BESS','QUIJOTE'], key=lambda m:abs(resultats[m][1]))
print("\n  Ranking (millor = frena més la caiguda en 500ms):")
for i,m in enumerate(ordre,1):
    print(f"    {i}. {MECANISMES[m]['nom']:>34}: RoCoF = {resultats[m][1]:.3f} Hz/s")

base_r=abs(resultats['BASE'][1])
print(f"\n  Millora vs BASE (reducció del RoCoF):")
for m in ['SYNTH','BESS','QUIJOTE']:
    millora=100*(1-abs(resultats[m][1])/base_r)
    print(f"    {MECANISMES[m]['nom']:>34}: {millora:+.1f}%")

print("\n" + "═"*72)
print("  VEREDICTE (conclusió estructural, no numèrica)")
print("═"*72)
print(f"""
  El RoCoF es juga en els primers 500 ms. Cada mecanisme actua amb la
  seua constant de temps τ:
    · BESS i synth-inertia: τ ≈ 10-20 ms → actuen DINS la finestra → frenen molt.
    · Quijote: τ ≈ 300 ms → actua a la VORA de la finestra → frena poc.

  La conclusió és ESTRUCTURAL (separació d'escales temporals), com va
  confirmar el revisor: τ_Quijote ≫ τ_electrònica, per tant Quijote
  arriba tard a la finestra del RoCoF i la seua contribució és petita.

  Açò NO depèn de la calibració fina del model: és ordre de magnitud.
  Encara que Quijote escale en nombre de turbines, el retard NO escala.

  → En RoCoF, l'electrònica de potència guanya per disseny temporal.
    Quijote no pot competir en aquesta finestra (física, no opinió).
""")
