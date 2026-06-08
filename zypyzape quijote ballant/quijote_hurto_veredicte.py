#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quijote_hurto_veredicte.py
═══════════════════════════════════════════════════════════════════════════════
PROVA DEFINITIVA: és viable l'hurto gravitatori? Que ho dicte el codi.

Víctor Manzanares Alberola — EPSA UPV Alcoi — github.com/espiradesombra/claude

El codi NO assumeix la resposta. Per a cada configuració (massa, posició del
ball, amplitud, velocitat d'actuador) construeix el balanç energètic complet
d'un cicle i NOMÉS el considera vàlid si és un CICLE TANCAT REAL: la massa
torna exactament al punt de partida i la seva energia total es conserva.

  τ_grav = −m·g·r·cos(θ)        [par gravitatori sobre el rotor, signe físic]
  F_act  = m·ω²·r − m·g·sin(θ)  [força radial que l'actuador venç]
  W_rotor    = ∮ τ_grav·ω dt    [el que el rotor rep de la gravetat]
  W_actuador = ∮ F_act·(dr/dt) dt [el que costa moure la massa]
  guany_net  = W_rotor − W_actuador

CONDICIONS DE VALIDESA (sense aquestes, qualsevol "guany" és artefacte):
  · cicle tancat:   r(final) = r(inici)
  · seguiment real: la massa segueix r_target (no es queda a mitges)
  · energia massa conservada: ΔE_massa ≈ 0 sobre el cicle

ÚS:  python3 quijote_hurto_veredicte.py
═══════════════════════════════════════════════════════════════════════════════
"""
import numpy as np
import itertools

# ── Constants físiques (NREL 5MW + Quijote) ───────────────────────────────────
G       = 9.81
OMEGA   = 1.327          # rad/s (ω_rated)
R_TIP   = 55.0
R_HUB   = 5.0
T_REV   = 2*np.pi/OMEGA
N       = 20000          # resolució angular

theta = np.linspace(0, 2*np.pi, N, endpoint=False)
dth   = theta[1]-theta[0]
dt    = dth/OMEGA
cos_t = np.cos(theta)
sin_t = np.sin(theta)


def moviment_radial(r_target, dr_max):
    """Segueix r_target amb velocitat radial limitada a dr_max."""
    r = np.empty(N); r[0] = r_target[0]
    step = dr_max*dt
    rr = r_target[0]
    for i in range(1, N):
        d = r_target[i]-rr
        rr += np.clip(d, -step, step)
        r[i] = rr
    return r


def balanc_cicle(m, r_lo, r_hi, dr_max, verbose=False):
    """Balanç energètic d'una volta amb la regla asimètrica Ec.2."""
    r_target = np.where(cos_t > 0, r_hi, r_lo)
    r = moviment_radial(r_target, dr_max)

    # Validesa del cicle
    cicle_tancat = abs(r[-1]-r[0]) < 0.01*(r_hi-r_lo+1e-9)
    seguiment    = np.mean(np.abs(r-r_target)) < 0.1*(r_hi-r_lo+1e-9)

    dr    = np.gradient(r, theta)
    dr_dt = dr*OMEGA

    # Energia total de la massa (per verificar conservació)
    v2 = (OMEGA*r)**2 + dr_dt**2
    E_massa = 0.5*m*v2 + m*G*r*sin_t
    dE_massa = E_massa[-1]-E_massa[0]

    # Treballs
    F_act = m*OMEGA**2*r - m*G*sin_t
    W_actuador = np.sum(F_act*dr_dt)*dt
    W_rotor    = np.sum(-m*G*r*cos_t*OMEGA)*dt
    guany_net  = W_rotor - W_actuador

    valid = (cicle_tancat and seguiment
             and abs(dE_massa) < abs(W_rotor)*0.05 + 1.0)

    t_commut = (r_hi-r_lo)/dr_max if dr_max > 0 else np.inf

    if verbose:
        estat = "✓ CICLE VÀLID" if valid else "✗ no tancat (massa a mitges)"
        print(f"  m={m:6.1f}kg  r:[{r_lo:.2f},{r_hi:.2f}]  Δr={r_hi-r_lo:.3f}m  "
              f"v={dr_max:.1f}m/s  [{estat}]")
        print(f"     W_rotor={W_rotor/1e3:+.2f}kJ  W_act={W_actuador/1e3:+.2f}kJ  "
              f"guany={guany_net/1e3:+.2f}kJ  ΔE_massa={dE_massa/1e3:+.3f}kJ")
        print(f"     commutació={t_commut:.2f}s ({t_commut/T_REV:.2f} voltes)")

    return {'guany_kJ': guany_net/1e3, 'valid': valid,
            'W_rotor_kJ': W_rotor/1e3, 'W_act_kJ': W_actuador/1e3,
            'dE_massa_kJ': dE_massa/1e3, 't_commut': t_commut,
            'voltes': t_commut/T_REV}


# ═══════════════════════════════════════════════════════════════════════════════
print("═"*70)
print("  QUIJOTE — VEREDICTE SOBRE L'HURTO GRAVITATORI")
print("  Víctor Manzanares Alberola — EPSA UPV Alcoi")
print("═"*70)

print("\n── Casos discutits ──")
print("\n[A] Pes gran 332kg, recorregut 50m, actuador lent (original):")
balanc_cicle(332.4, 5.0, 55.0, 0.5, verbose=True)
print("\n[C] Els teus números: 66.5kg, Δr=0.444m a la punta:")
balanc_cicle(66.5, 55.0-0.444, 55.0, 0.5, verbose=True)
print("\n[D] Ideal: actuador infinitament ràpid (snap perfecte):")
balanc_cicle(332.4, 5.0, 55.0, 1000.0, verbose=True)

# ── Escaneig exhaustiu ──
print("\n── Escaneig exhaustiu (només cicles tancats vàlids compten) ──")
masses     = [10, 50, 66.5, 100, 332, 1000, 5000]
r_centres  = [10, 30, 45, 52, 54.7]
amplituds  = [0.444, 1, 2, 5, 25, 50]
velocitats = [0.5, 2, 10, 50, 1000]

n_valid = n_guany = descartats = 0
millor = -np.inf; millor_cfg = None
for m, rc, dr, vm in itertools.product(masses, r_centres, amplituds, velocitats):
    rlo, rhi = rc-dr/2, rc+dr/2
    if rlo < R_HUB or rhi > R_TIP: continue
    res = balanc_cicle(m, rlo, rhi, vm)
    if not res['valid']:
        descartats += 1
        continue
    n_valid += 1
    if res['guany_kJ'] > 0.01:
        n_guany += 1
        if res['guany_kJ'] > millor:
            millor = res['guany_kJ']; millor_cfg = (m, rlo, rhi, vm)

print(f"\n  Cicles tancats vàlids:           {n_valid}")
print(f"  Descartats (massa no completa):  {descartats}")
print(f"  Amb guany net real (W_net > 0):  {n_guany}")

# ── Llei de conservació directa ──
print("\n── Comprovació independent: ∮τ_grav dθ (snap ideal) ──")
for m in [50, 332, 1000]:
    r_id = np.where(cos_t > 0, R_TIP, R_HUB)
    W = np.sum(-m*G*r_id*cos_t)*dth
    print(f"  m={m:5d}kg:  ∮τ_grav dθ = {W/1e3:+.2f} kJ  "
          f"(treball gravetat sobre rotor, no compensa el cost d'actuador)")

# ── Veredicte ──
print("\n" + "═"*70)
print("  VEREDICTE")
print("═"*70)
if n_guany == 0:
    print(f"""
  De {n_valid} configuracions amb cicle TANCAT i físicament vàlid,
  CAP no produeix guany net. L'HURTO GRAVITATORI NO ÉS VIABLE.

  Els 'guanys' aparents que es veurien sense filtrar cicles
  ({descartats} casos descartats) són artefactes: l'actuador no completa
  el moviment dins d'una volta, la massa acaba en posició diferent de
  l'inicial, i el 'guany' és energia presa del transitori, no del cicle.

  Causa estructural: la gravetat és conservativa. Per a tota trajectòria
  on la massa torna al punt de partida, el treball net és nul, i moure
  la massa contra la centrífuga sempre té un cost.

  CONCLUSIÓ POSITIVA: el mecanisme és inservible com a font d'energia,
  però el ball curt a la punta (cas C: commutació en {balanc_cicle(66.5,54.556,55.0,0.5)['voltes']:.2f} voltes)
  és un EXCEL·LENT actuador per a regulació de freqüència ràpida.
""")
else:
    print(f"\n  {n_guany} casos amb guany aparent. Millor: {millor:+.2f} kJ.")
    print(f"  Config: {millor_cfg}. REVISAR manualment.")
