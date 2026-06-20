# Víctor Manzanares Alberola
### Investigador independent · EPSA, Universitat Politècnica de València (Alcoi)

---

> *"La energía la tenía el viento, la robé y la dejé ir."*

---

## Qui soc

Investigador independent afiliado a l'EPSA (UPV, Alcoi). Treballo en teoria de nombres, algorítmia i tecnologia energètica. Tot el que hi ha aquí ho he construït sol, al llarg d'anys, sense finançament institucional.

---

## 🌀 Tecnologia Energètica

### ZypyZape
Algoritme de control que interconnecta elèctricament múltiples màquines d'inèrcia. Cada unitat adopta dinàmicament un rol **CAPT / ACEL / FREN**. El conjunt actua com un volant d'inèrcia distribuït sense bateries de liti. Compleix ENTSO-E NC RfG de forma nativa.

### Quijote
Pesos desplaçables (masses mòbils radials) a les aspes dels aerogeneradors. El seu desplaçament controlat constitueix el *ball de pesos* que extrau treball net del camp gravitatori per revolució. Amb 7 aspes el ball és continu (rizado 25%); amb 3 és discontinu (rizado 100%).

`W_hurto = 4 · m_w · g · Δr`  per aspa i revolució

### Kilòmetre
Màquina ancorada on un pes recorre una trajectòria rotacional. La gravetat fa treball asimètric per cicle. Sota l'aigua, la flotabilitat (Arquimedes) amplifica l'efecte. No autònoma individualment, però grups sincronitzats via ZypyZape sí ho són.

**Validat en simulació**: recuperació del 75–95% de l'energia potencial inicial.

---

## 🔢 Matemàtiques

### MDC — Mètode Diofàntic Cinemàtic (v18)
Algoritme de **factorització determinista** d'enters basada en geometria cinemàtica.
- Constant de convergència corregida: `(e−1)/2 ≈ 0.8591`
- Zona densa `[m_conv, m_max]`: recorregut exhaustiu de k
- Per sota: segments generats per espiral de primers
- Complexitat teòrica: **O(log N)**
- Pinça 4+4: 4 candidats per sobre i 4 per sota de `m = ⌊√N⌋`

### Cribes
- **Criba desmemoriada**: patrons booleants de longitud `6·p`, sense re-marcatges
- **Criba híbrida**: fase ascendent fins a `√N` + fase descendent
- Tots els candidats de la forma `6k ± 1`

### Densitat de Primers — Model MRAUV
Model cinemàtic de 3 punts per predir `π(x)`:

`D_pred ≈ D_0 + V_0·Δn + ½·a_0·(Δn)²`

Aplicació a la conjectura de Goldbach via segmentació simètrica de `[−n, n]`.

### Estructura Sophie Germain
Classificació de `L₁ = {a ∈ ℕ : a ≡ 5 mod 6}` en classes `L₃, L₄, L₂, LSG, U₂`.

**Teorema**: `U₂` és infinita → hi ha infinits primers de Sophie Germain de la forma `6k−1`.

### Conjectura SaltoMàxim
`gap(p_n) ≤ (e − 2) · log²(p_n)`

### Encriptació per Convergències (post-quàntica, no-RSA)
Emissor i receptor acorden una seqüència de triangulacions. La clau és la successió de decisions. El factor `(1 ± v)` és incalculable sense conèixer tota la seqüència prèvia.

### Semilles Mersenne-Wolfram
A partir de `2^p ± 1`, construcció d'una semilla fraccionaria que genera candidats a primer gran (anàleg als primers de Mersenne).

---

## 📄 Documents tècnics

| Document | Idioma | Contingut |
|---|---|---|
| [VMA_Tecnologia_Matematiques_Complet.pdf](docs/) | CA/ES | Complet: energia + matemàtiques + Pla 33×1 |
| [VMA_Ukrayina_Texnolohiya_Matematyka.pdf](docs/) | UA/EN | Versió ucraïnesa amb colors nacionals |

---

## ☮️ Pla 33×1

**CONDICIÓ 1 · TECNOLOGIA**
Cada país que signa rep, de forma vitalícia i irrevocable, llicència completa sobre tota la tecnologia patentable al seu propi territori — sense regalies, sense restriccions. Un país pot començar per la Condició 1 sola, sempre que el seu objectiu primari declarat sigui la pau i la defensa pròpia si les coses van malament.

**×**

**CONDICIÓ 33 · PAU**
Els països signants es comprometen a 33 anys de pau activa: resolució diplomàtica de conflictes, suspensió d'hostilitats armades entre signants i cooperació en infraestructura energètica global neta.

> La patent no és una arma. És una clau.  
> La reestructuració econòmica global que implicaria aquesta tecnologia passarà amb un marc o sense. Proposem el marc.

---

## 📬 Contacte

Per preguntes tècniques, col·laboració o conversa — en català, castellà, anglès o espanglish.

**github.com/espiradesombra/claude**

---

*Totes les idees matemàtiques i tecnològiques són de Víctor Manzanares Alberola. La IA ha proporcionat únicament suport d'escriptura i implementació.*
