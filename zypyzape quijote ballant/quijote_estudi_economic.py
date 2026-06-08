#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quijote_estudi_economic.py
═══════════════════════════════════════════════════════════════════════════════
ESTUDI ECONÒMIC: drets, costos de fabricació, temps de vida i viabilitat
del sistema Quijote (inèrcia variable mecànica) com a servei de regulació.

Ancorat en dades de mercat reals (2025-2026):
  · Mercat d'inèrcia alemany (gener 2026): remuneració preu fix, paga BESS
    grid-forming. Preu inèrcia sintètica: pic £4,73/MW·s, sovint →0 per excés.
  · FFR (fast frequency response): ~$0,3/MWh en alta penetració renovable.
  · Competència directa: bateries + inversors (sense parts mòbils).

L'estudi NO inventa precisió: usa RANGS i anàlisi de sensibilitat. El número
exacte variarà; el que importa és l'ordre de magnitud i el llindar de viabilitat.

Víctor Manzanares Alberola — EPSA UPV Alcoi
═══════════════════════════════════════════════════════════════════════════════
"""
import numpy as np

print("█"*72)
print("  QUIJOTE — ESTUDI ECONÒMIC (drets, fabricació, vida, viabilitat)")
print("█"*72)

# ─────────────────────────────────────────────────────────────────────────────
# 1. CAPACITAT TÈCNICA: quanta inèrcia/energia pot aportar Quijote
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "═"*72)
print("  1. CAPACITAT TÈCNICA DEL DISPOSITIU")
print("═"*72)
# Per turbina NREL 5MW, massa mòbil m/pala
R_TIP=55.0; R_HUB=5.0; N_PALES=3; OMEGA=1.327
for m in [500, 1000, 2000]:
    # energia cinètica modulable = ½·ΔJ·ω²
    dJ = N_PALES*m*(R_TIP**2 - R_HUB**2)
    E_modulable = 0.5*dJ*OMEGA**2     # J
    E_kWh = E_modulable/3.6e6
    # constant d'inèrcia equivalent H (s) afegida, sobre 5 MVA
    H_afegit = E_modulable/(5e6)
    # MW·s d'inèrcia (energia cinètica disponible)
    MWs = E_modulable/1e6
    print(f"  m={m}kg/pala: ΔJ={dJ:.2e} kg·m²  E_modulable={E_kWh:.2f} kWh  "
          f"= {MWs:.2f} MW·s  (H+={H_afegit:.3f}s)")
print("\n  Referència: una bateria de 5 MW / 1 MWh aporta molt més,")
print("  i respon en mil·lisegons sense parts mòbils.")

# ─────────────────────────────────────────────────────────────────────────────
# 2. INGRESSOS POTENCIALS (mercat de regulació real)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "═"*72)
print("  2. INGRESSOS POTENCIALS (dades mercat 2025-2026)")
print("═"*72)
# Preu inèrcia: pic 4.73 £/MW·s, sovint baix. Usem rang.
# El servei es 'lloga' per disponibilitat, no per activació contínua.
MWs_disp = 0.5*N_PALES*1000*(R_TIP**2-R_HUB**2)*OMEGA**2/1e6  # m=1000kg
print(f"\n  Capacitat oferible (m=1000kg): {MWs_disp:.2f} MW·s per turbina")
print(f"\n  {'escenari preu':>28}{'€/MW·s/any':>14}{'ingrés/turbina/any':>22}")
# remuneració per disponibilitat: assumim pagament anual per MW·s ofert
# rang basat en mercats reals (molt incert)
for nom, preu_anual in [("optimista (mercat escàs)", 2000),
                        ("mitjà", 500),
                        ("realista (excés oferta)", 100),
                        ("pessimista (preu→0)", 20)]:
    ingres = MWs_disp * preu_anual
    print(f"  {nom:>28}{preu_anual:>14}{ingres:>20.0f} €")
print("\n  ⚠ El mercat alemany 2026 paga BESS; un dispositiu mecànic ha de")
print("    competir-hi i certificar-se com a 'grid-forming' equivalent.")

# ─────────────────────────────────────────────────────────────────────────────
# 3. COSTOS DE FABRICACIÓ (CAPEX)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "═"*72)
print("  3. COSTOS DE FABRICACIÓ (CAPEX per turbina)")
print("═"*72)
# Components del sistema Quijote per turbina (3 pales)
components = [
    ("Masses mòbils (3×1000kg, acer/fluid)", 3*1000*3.0),       # €/kg material+forma
    ("Rails/guies radials reforçats (3 pales)", 3*8000),
    ("Actuadors lineals d'alta força (3)", 3*15000),
    ("Sistema hidràulic/elèctric de potència", 25000),
    ("Sensors, encoders, control embarcat", 12000),
    ("Reforç estructural pala (fatiga)", 3*10000),
    ("Electrònica de potència + certificació GFM", 30000),
    ("Integració, instal·lació, proves", 40000),
]
total_capex = sum(c for _,c in components)
print(f"\n  {'component':>42}{'cost (€)':>14}")
for nom, cost in components:
    print(f"  {nom:>42}{cost:>14,.0f}")
print(f"  {'─'*56}")
print(f"  {'CAPEX TOTAL per turbina':>42}{total_capex:>14,.0f} €")
print(f"\n  Rang d'incertesa: ±40% → [{total_capex*0.6:,.0f} – {total_capex*1.4:,.0f}] €")

# ─────────────────────────────────────────────────────────────────────────────
# 4. COSTOS OPERATIUS (OPEX) i VIDA ÚTIL
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "═"*72)
print("  4. OPEX, FATIGA i VIDA ÚTIL")
print("═"*72)
# El moviment cíclic radial genera fatiga → manteniment + vida limitada
cicles_per_any = OMEGA/(2*np.pi)*3600*24*365 * 0.3  # 30% del temps actiu
print(f"\n  Cicles de càrrega radial/any (30% actiu): {cicles_per_any:.2e}")
print(f"  → Fatiga severa: les masses mòbils i rails són peces de desgast.")
vida_anys = 12   # optimista per a mecanisme amb tant cicle (vs 25 turbina)
opex_anual = total_capex*0.06   # 6% CAPEX/any (alt, per parts mòbils)
print(f"  Vida útil estimada del sistema Quijote: {vida_anys} anys")
print(f"    (vs 25 anys de la turbina → cal substituir-lo ~2 vegades)")
print(f"  OPEX anual (manteniment, 6% CAPEX): {opex_anual:,.0f} €/any")
print(f"  Temps de fabricació/integració: 6-12 mesos (enginyeria a mida)")

# ─────────────────────────────────────────────────────────────────────────────
# 5. VAN / PAYBACK (anàlisi de sensibilitat)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "═"*72)
print("  5. VIABILITAT FINANCERA (VAN a 12 anys, descompte 7%)")
print("═"*72)
def van(capex, ingres_anual, opex_anual, vida, r=0.07):
    flux = -capex
    for t in range(1, vida+1):
        flux += (ingres_anual - opex_anual)/(1+r)**t
    return flux

print(f"\n  CAPEX={total_capex:,.0f}€, OPEX={opex_anual:,.0f}€/any, vida={vida_anys}anys")
print(f"\n  {'ingrés regulació/any':>24}{'VAN (€)':>16}{'veredicte':>20}")
for nom, preu_anual in [("optimista", 2000),("mitjà", 500),
                        ("realista", 100),("pessimista", 20)]:
    ingres = MWs_disp*preu_anual
    v = van(total_capex, ingres, opex_anual, vida_anys)
    vd = "VIABLE" if v>0 else "PÈRDUA"
    print(f"  {nom+' ('+str(int(ingres))+'€/any)':>24}{v:>16,.0f}{vd:>20}")

# llindar: quin ingrés anual cal per a VAN=0
# VAN=0 → capex = (ingres-opex)·Σ1/(1+r)^t
factor = sum(1/(1.07)**t for t in range(1,vida_anys+1))
ingres_breakeven = total_capex/factor + opex_anual
preu_breakeven = ingres_breakeven/MWs_disp
print(f"\n  LLINDAR DE RENDIBILITAT (VAN=0):")
print(f"    ingrés mínim = {ingres_breakeven:,.0f} €/any per turbina")
print(f"    = preu de {preu_breakeven:,.0f} €/MW·s/any pel servei")
print(f"    (preu de mercat real: 20-2000 €/MW·s/any, sovint <500)")

# ─────────────────────────────────────────────────────────────────────────────
# 6. VEREDICTE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "═"*72)
print("  VEREDICTE ECONÒMIC")
print("═"*72)
print(f"""
  CAPACITAT:  Quijote aporta ~{MWs_disp:.1f} MW·s d'inèrcia per turbina (m=1000kg).
              Una bateria petita aporta molt més, en ms, sense desgast.

  CAPEX:      ~{total_capex/1000:.0f} k€/turbina (rang {total_capex*0.6/1000:.0f}-{total_capex*1.4/1000:.0f}).
  VIDA:       ~{vida_anys} anys (fatiga per cicle radial) vs 25 de la turbina.
  LLINDAR:    Cal ~{preu_breakeven:,.0f} €/MW·s/any per ser rendible.
              El mercat paga sovint <500 €/MW·s/any → per sota del llindar.

  CONCLUSIÓ HONESTA:
  · Com a NEGOCI de regulació, Quijote mecànic difícilment competeix amb
    bateries grid-forming: menys capacitat, més lent, desgast, vida curta.
  · El mercat d'inèrcia (Alemanya 2026) està dissenyat per a BESS, no per
    a massa mòbil. La certificació grid-forming d'un sistema mecànic és
    un obstacle addicional gran.
  · NOMÉS tindria sentit com a RETROFIT barat si les masses i actuadors
    aprofiten estructura existent i el preu de mercat puja molt (escassetat).

  RECOMANACIÓ: el valor de Quijote NO és comercial com a producte d'inèrcia.
  El seu valor real és ACADÈMIC i de PRIOR ART: el paper honest sobre
  inèrcia variable, l'auditoria energètica rigorosa, i la metodologia de
  refutació. Això és publicable i defensable; el producte físic, no.
""")
