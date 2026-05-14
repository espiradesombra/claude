# ZypyZape + Quijote + Kilómetro
### Víctor Manzanares Alberola — EPSA, Universitat Politècnica de València (Alcoi)

---

## Què és aquest repositori?

Implementació i documentació dels sistemes d'inèrcia sintètica mecànica per a turbines eòliques i xarxes energètiques.

**Idea central:** convertir els aerogeneradors existents en una bateria mecànica distribuïda sense litio, sense electrònica complexa, sense violar cap llei termodinàmica.

---

## Els tres sistemes

### 1. ZypyZape
Control de parc eòlic per acoblament de fase tipus Kuramoto (K = 0.10).
- N turbines s'intercanvien energia cíclicament (f = 0.4 Hz)
- Una accelera mentre l'altra frena → inèrcia sintètica col·lectiva
- Millora nadir de freqüència i ROCOF en pertorbacions de xarxa
- Validat: η_global = 98%, millora nadir = +0.002 Hz/mòdul

### 2. Quijote
Modulació d'inèrcia per massa desplaçable dins de les pales.
- Fluid Fe+oli (ρ = 3386 kg/m³) com a massa lliscant en canal Ø5cm
- r_q ∈ [5m, 55m] — la massa es desplaça radialment
- Equació de moviment: `m_q·r̈ = F_c + F_ctrl + F_fric`
- Límit posicional original: `|ṙ_q| ≤ ω·r_q`
- Validat sobre NREL 5MW: **+1.4–1.5% P_grid**, P_gen = 128–226 W

**3 pales vs 7 pales:**
- ΔE màxim: proporcional a N (7p = 2.33× més que 3p)
- Factor equivalent a √3 elèctric: `Factor_N = 2·sin(π/N)`
- J_total constant amb ball sinusoïdal (propietat matemàtica fonamental)

### 3. Kilómetro
Buffer inercial distribuït per a xarxes energètiques a gran escala.
- N unitats acoblades via ZypyZape → cobertura contínua de fase
- 33 unitats (Plan 33x1): separació 11° → quasi cobertura total del cicle
- No genera energia: redistribueix la que prové de fonts externes (vent, corrents)
- Versió submarina: camps gravetat + flotabilitat com a "rail" del cicle

---

## Fitxers

| Fitxer | Descripció |
|--------|-----------|
| `gemelo_v94.py` | Gemell digital v9.4 — model base validat |
| `gemelo_v941.py` | v9.4.1 — màquina d'estats STABLE/VALLEY/RECOVERY |
| `gemelo_v942.py` | v9.4.2 — TRANSITION, K_q=8e4, THRESH=0.2 |
| `gemelo_virtual_aprendre.py` | Versió didàctica — paràmetres configurables |
| `zypyzape_twin_v4_8_quijote.py` | Gemell complet v4.8 amb Quijote integrat |
| `01_CONTEXT_ZYPYZAPE.txt` | Context complet del model (paràmetres, física) |
| `02_MATH_QUIJOTE_3vs7.txt` | Matemàtica del Quijote: 3 vs 7 pales |
| `QUIJOTE_ZYPYZAPE_paper_llibre_bilingual.docx` | Paper científic (Valencià/Anglès) |
| `Kilometre_ZypyZape_doc_tecnic.docx` | Document tècnic Kilómetro |

---

## Paràmetres clau (v4.8)

```
Turbina:    R=60m, J=5e6 kg·m², S_nom=2.5MW, N=3 pales
Quijote:    M_Q=4kg/pala, r∈[5m,60m], fluid Fe+oli ρ=3386kg/m³
ZypyZape:   K=0.10 (subcrític), f_cicle=0.4Hz, P_ZZ=13% S_nom
Xarxa:      2GW, H_sys=4s, f0=50Hz, dP=-100MW (pertorbació ref.)
```

---

## Resultats validats

| Paràmetre | 3 pales | 7 pales |
|-----------|---------|---------|
| Millora P_grid | +1.4% | +1.5% |
| P_buf oscil·lant | 12 kW | 16 kW |
| P_gen hidràulic | 128 W | 226 W |
| ΔE Quijote/cicle | 28 kJ | 65 kJ |
| Millora nadir | +0.002 Hz/mòdul | — |

---

## Física fonamental

**El sistema NO crea energia.** Redistribueix l'energia del vent (font no conservativa) en el temps, reducint pèrdues i millorant l'estabilitat.

Equació d'estat del conjunt:
```
J_i(t)·ω̇_i = τ_ext,i + τ_coupling,i − τ_loss,i
τ_coupling,i = K·Σⱼ sin(θⱼ − θᵢ)   [Kuramoto mecànic]
∫P_total(t)dt ≤ E_input             [Primera llei TD]
```

---

## Publicació

- Paper Quijote: pendent de DOI (Zenodo / arXiv)
- Document Kilómetro: pendent de DOI

---

## Llicència

© Víctor Manzanares Alberola, 2026. Tots els drets reservats.
Les idees matemàtiques i d'enginyeria són originals de l'autor.
El codi es publica amb finalitats de reproducibilitat científica.

---

*EPSA — Escola Politècnica Superior d'Alcoi | Universitat Politècnica de València*
