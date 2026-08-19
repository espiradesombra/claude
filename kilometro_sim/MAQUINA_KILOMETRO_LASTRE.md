# Máquina Kilómetro — Ascensor de lastre con pernos (documentación completa)

Documento de referencia de la máquina tal como se ha modelado y simulado en `kilometro_sim/`, alineada con el diseño VMA (3 flota / 4 se hunde, perneo, enjambre, recarga de distancia en un solo lado).

**No es un motor perpetuo ni un “robo de gravedad”.**  
**Sí es** un sistema de **almacenamiento / conmutación de lastre** en medio fluido, con valor de ingeniería como buffer, embrague de densidad y logística de masas.

---

## 1. Idea en una frase

> Un **objeto** (módulo Kilómetro) recorre una **guía** entre dos cotas.  
> Al **enganchar un peso extra** se hunde y puede entregar energía; al **soltarlo** flota y vuelve casi sin esfuerzo.  
> Los **pernos** solo deciden si cada peso cuelga del **objeto** o del **recorrido**.  
> La energía sale de la **altura de los lastres en el almacén**, no de un truco geométrico impar/par.

---

## 2. Piezas de la máquina

### 2.1 Objeto (módulo / “Kilómetro”)

- Cuerpo que se desplaza a lo largo del **recorrido** (guía).
- Lleva de forma nominal una configuración de pesos tal que:
  - con **3 pesos** (conteo total en el objeto) → **flota** (flotabilidad ≥ neutra),
  - con **4 pesos** → **se hunde**.
- En las **puntas** puede llevar 2 pesos fijos (geometría + posibilidad de giro 180° para presentar pernos al recorrido).
- Puede **recargar distancia / geometría** (extender o rearmar el “recorrido del objeto”) siempre en la **misma estación** (lado S, cota alta en el modelo simulado).

### 2.2 Recorrido (guía / pista)

- Estructura fija o en bucle con al menos dos cotas útiles:
  - **ALTA** (`h₁`),
  - **BAJA** (`h₂`),
  - \(\Delta h = h_1 - h_2 > 0\).
- Actúa como:
  - riel de guiado del objeto,
  - **estantería de lastres** (stock arriba / aparcamiento abajo),
  - soporte cuando un peso está **perneado a la guía** y no al objeto.

### 2.3 Pesos (lastres discretos)

- Masas sueltas, idénticas en el modelo simple (`m_peso`).
- Cada una puede estar en uno de estos sitios:
  1. **Stock ALTA** (listos para enganchar),
  2. **Sobre el objeto** (aumentan densidad del módulo),
  3. **Stock BAJA** (aparcados tras un ciclo de descarga).
- La “batería gravitacional” es literalmente: **cuántos kg hay todavía en cota alta**.

### 2.4 Sistema de pernos (dos por peso)

Cada peso tiene **dos pernos** (embrague de lastre):

| Perno | Une | Estado típico |
|--------|-----|----------------|
| **R** | Peso ↔ **Recorrido** | Peso aparcado en la guía |
| **O** | Peso ↔ **Objeto** | Peso viaja con el módulo |

Reglas de seguridad:

| R | O | Significado |
|---|---|-------------|
| 1 | 0 | Peso solo en recorrido (almacén / aparcamiento) |
| 0 | 1 | Peso solo en objeto (módulo más denso) |
| 1 | 1 | **Make-before-break**: transferencia segura un instante |
| 0 | 0 | Prohibido (peso suelto) |

El perneo **no crea energía**. Solo cambia **quién soporta el kg** a una cota dada:

```text
mismo peso, misma altura h
  colgado del objeto  ≡  colgado del recorrido
  → mismo potencial m g h (misma “factura” gravitatoria)
```

### 2.5 Fluido (mar / agua)

- Proporciona **flotabilidad de Arquímedes**.
- Con 3 pesos el objeto se comporta como casi neutro o flotante → subida / reintegración barata.
- Con 4 pesos la densidad media supera la del fluido → descenso “solo”.
- El fluido es un **riel pasivo de bajo coste de control**, no un generador.
- **No sube solo los lastres** de abajo a arriba: eso lo hace un lift, la guía motorizada u otro proceso.

---

## 3. Umbral 3 flota / 4 se hunde

### 3.1 Por qué es potente

Estar cerca de la flotabilidad neutra hace que **un solo peso discreto** conmute el comportamiento:

```text
n_pesos_en_objeto ≤ 3  →  flota / se reintegra fácil
n_pesos_en_objeto ≥ 4  →  se hunde / se aleja del punto de flotación
```

Eso evita bombas hidráulicas o grandes variaciones de volumen: basta un **click de pernos**.

### 3.2 Interpretación en la simulación

Parámetros por defecto (`sim_enjambre_pesos_impar.py`):

| Símbolo | Valor | Rol |
|---------|-------|-----|
| `n_float` | 3 | máximo para flotar |
| `n_sink` | 4 | mínimo para hundirse |
| `n_puntas` | 2 | pesos en extremos (diseño geométrico) |
| `m_peso` | 10 kg | masa de cada lastre |
| `m_obj_base` | 30 kg | masa estructural del módulo |
| `delta_h` | 15 m | desnivel útil |
| `n_stock_alta` | 5 | lastres iniciales en cota alta |
| `eta_gen` | 0.85 | eficiencia al recuperar en bajada |
| `eta_lift` | 0.90 | eficiencia al subir un lastre |
| `E_perno` | 1.5 J | coste eléctrico por conmutación de perno |
| `drag_frac` | 0.06 | pérdidas hidrodinámicas (fracción) |

Al inicio de cada ciclo el objeto se calibra en **3 pesos** (flota) en **cota alta**.

---

## 4. Ciclo operativo detallado

### 4.1 Modo A — Almacén abierto (sin subir pesos)

Recarga de distancia **siempre en el mismo lado** (estación **S** = cota alta).

```text
[0] Estado reposo
    - Objeto en ALTA, n = 3, flota
    - Stock ALTA = N, stock BAJA = 0 (o lo que haya)
    - Distancia/geometría del objeto lista en estación S

[1] Enganche (perneo)
    - Peso de stock ALTA: R=1,O=0  →  make-before-break  →  R=0,O=1
    - n: 3 → 4
    - Coste: ~4 × E_perno (dos pernos en enganche + dos en soltar al final del ciclo)
    - Objeto se hunde

[2] Bajada
    - Conjunto desciende Δh
    - Se puede recuperar energía (freno regenerativo, generador en el recorrido, etc.)
    - Contabilidad HONESTA: solo el lastre EXTRA que no vuelve con el objeto
      aporta PE “gastable” de inventario:
        W_rec ≈ η_gen · (1 − drag) · m_peso · g · Δh
      (por cada peso extra que se queda abajo; en el modelo simple = 1 peso/ciclo)
    - La PE del cuerpo + 3 pesos que suben después por boyancia
      NO se cuenta como generación neta durable (se anula al flotar)

[3] Soltar / perneo al recorrido en BAJA
    - En cota baja (o tras giro 180° presentando la otra punta):
      O=1,R=0  →  make-before-break  →  O=0,R=1
    - Peso queda en stock BAJA
    - n: 4 → 3
    - Objeto flota otra vez

[4] Subida del objeto
    - Boyancia lleva el módulo a cota ALTA
    - Coste eléctrico del módulo ≈ 0 en el modelo ideal de flotación
    - En estación S: “recarga de distancia” (geometría del Kilómetro otra vez lista)

[5] Inventario
    - Stock ALTA −= 1
    - Stock BAJA += 1
    - El OBJETO parece “como al inicio”
    - El SISTEMA no: hay un kg menos arriba y uno más abajo
```

Repetir hasta `stock_ALTA = 0` → **STOP** (no se puede pasar a 4 pesos).

**Número de ciclos máximos en modo A:**

\[
N_{\text{ciclos}} = N_{\text{stock alta inicial}}
\]

Par o impar del número final = **tamaño del almacén**, no magia física.

### 4.2 Modo B — Ciclo cerrado (lift continuo del lastre)

Igual que A, pero al final de cada ciclo (o de forma continua):

```text
[6] Lift / restitución
    - Subir 1 peso de BAJA → ALTA
    - W_lift = (m_peso · g · Δh) / η_lift
    - Stock se mantiene
    - El sistema puede ciclar indefinidamente en teoría
```

Balance por ciclo (orden de magnitud):

\[
W_{\text{gen}} \approx \eta_{\text{gen}}\, m_{\text{peso}}\, g\, \Delta h
\]
\[
W_{\text{lift}} = \frac{m_{\text{peso}}\, g\, \Delta h}{\eta_{\text{lift}}}
\]

Si \(\eta_{\text{gen}} < 1\) y \(\eta_{\text{lift}} < 1\):

\[
W_{\text{lift}} > W_{\text{gen}}
\quad\Rightarrow\quad
\text{surplus por ciclo} < 0
\]

En la simulación de referencia: ~**−465 J/ciclo** con los parámetros por defecto.

### 4.3 Diagrama de bloques

```text
                    ┌─────────────────────────┐
                    │  ESTACIÓN S (siempre)   │
                    │  recarga distancia      │
                    │  enganche lastre        │
                    └───────────┬─────────────┘
                                │
         ALTA  ═════════════════╪═════════════════
              stock ● ● ● ● ●   │
                                │
                           ┌────▼────┐
                           │ OBJETO  │  3 flota / 4 se hunde
                           │+pernos  │
                           └────┬────┘
                                │ bajada (genera si hay lastre extra)
                                ▼
         BAJA  ═════════════════╪═════════════════
              aparca ● ● ●      │
                                │
                    sube ligero (boyancia) ───► vuelve a S
```

---

## 5. Dos “estados iniciales” (el malentendido clave)

| Estado | Qué incluye | ¿Se restaura cada ciclo en modo A? |
|--------|-------------|-------------------------------------|
| **Objeto** | n=3, flota, geometría/distancia, cota alta | **Sí** (en estación S) |
| **Sistema** | todos los pesos en sus cotas + stocks | **No** (cada ciclo mueve 1 kg ALTA→BAJA) |

El perneo “devuelve el Kilómetro al traje inicial”, **no** devuelve el peso a la cota alta.

Por eso parece magia si solo miras el objeto; se desmonta si miras la **estantería**.

---

## 6. Contabilidad energética (sin trampa)

### 6.1 Campo gravitatorio

En ciclo cerrado de posiciones de **cada masa**:

\[
\oint \mathbf{F}_g \cdot d\mathbf{r} = 0
\]

La gravedad es un “muelle ideal”: lo que da al bajar lo exige al subir **esa misma masa**.

### 6.2 Modo almacén abierto (descarga de batería)

Energía útil aproximadamente:

\[
E_{\text{útil}} \approx
N \cdot \eta_{\text{gen}} \cdot (1 - f_{\text{drag}}) \cdot m_{\text{peso}} \cdot g \cdot \Delta h
- E_{\text{pernos}}
\]

Interpretación: estás **descargando** la PE que ya tenían los \(N\) lastres en la estantería alta  
(alguien los subió antes: red, grúa, otro proceso, oleaje, etc.).

### 6.3 Cierre del almacén

Si al final subes los \(N\) pesos aparcados abajo:

\[
E_{\text{reset}} = N \cdot \frac{m_{\text{peso}} \cdot g \cdot \Delta h}{\eta_{\text{lift}}}
\]

Saldo típico:

\[
E_{\text{útil}} - E_{\text{reset}} < 0
\]

(por \(\eta < 1\) y drag).

### 6.4 Resultados numéricos de la simulación de referencia

Script: `sim_enjambre_pesos_impar.py`  
Parámetros: `m_peso=10 kg`, `Δh=15 m`, `N=5`, `η_gen=0.85`, `η_lift=0.90`.

| Concepto | Valor aprox. |
|----------|----------------|
| \(m g h\) por peso | 1472 J |
| W_rec honesto / ciclo | 1176 J |
| Pernos / ciclo | 6 J |
| Surplus honesto / ciclo (sin lift) | ~1170 J |
| Ciclos hasta STOP (modo A) | **5** (= stock alta) |
| Surplus acum. A sin reset | ~5849 J |
| Coste subir 5 pesos | ~8175 J |
| Saldo A al cerrar almacén | **~−2326 J** |
| Surplus / ciclo modo B (con lift) | **~−465 J** |

---

## 7. “Siempre impar” y paridad

### 7.1 Lo que la sim demostró

- Con stock inicial \(N\), sin lift, el sistema se detiene en el ciclo **\(N\)**.
- Si \(N=5\) paras en 5 (impar); si \(N=6\) paras en 6 (par).
- La paridad **no** es un efecto físico especial.

### 7.2 Lo que sí es “siempre a medias” en el diseño

Cada ciclo de carga mueve lastre en **un solo sentido**:

```text
peso:  ALTA → objeto → BAJA
```

No hay, en ese mismo acto, el retorno del peso.  
Por eso el **inventario** queda siempre “abierto” hasta que:

- se acaba el stock alto, o  
- un lift / la guía / el enjambre **reponen** pesos arriba.

Eso es gestión de almacén, no sobreunidad.

---

## 8. Variantes del diseño (conversación VMA)

### 8.1 Giro 180° del objeto

Durante el movimiento o abajo, el objeto rota para presentar la otra **punta** al recorrido y transferir el peso al “otro extremo” de la guía.

- Útil como **embrague de fase** y presentación de pernos.
- **No** cambia el potencial: solo la altura vertical del kg importa.

### 8.2 Enjambre de módulos

Varios objetos comparten la misma guía y el mismo stock de pesos:

- unos bajan densos,
- otros suben neutros,
- la guía acumula o reparte lastres.

Misma contabilidad global: suma de \(m_i g \Delta h_i\) de cada lastre.

### 8.3 Doble Kilómetro 90° + perneo inercial

Línea paralela de trabajo (otro script: `doble_km_90_perneo.py`):

- dos módulos desfasados,
- perneo/desperneo como **supercondensador mecánico** de picos,
- regeneración **solo en fase gravitatoria favorable** del eje,
- métrica \(P_{\text{pico}}/P_{\text{media}}\), no \(W_{\text{net}}>0\) mágico.

Complementa el lastre hidrostático: uno es **inventario de PE**, el otro es **buffer de potencia/inercia**.

### 8.4 Recorrido reduccionista (radio variable)

Script: `recorrido_reduccionista.py`

- \(r\) pequeño en tramos desfavorables, grande en favorables,
- o rearmar radio **en parado** (sin centrífuga).

Optimiza peajes eje/actuador; **no** crea energía neta en ciclo cerrado.

### 8.5 Fluido + volumen / flotabilidad casi neutra

Si además se varía volumen (tipo glider):

- bajar “solo” cerca del neutro es barato en fuerza,
- el peaje del ciclo de volumen es \(P\cdot\Delta V\) (más caro en profundidad).

Misma lección: control barato ≠ generación neta.

---

## 9. Qué es y qué no es esta máquina

| Afirmación | Veredicto |
|------------|-----------|
| Motor perpetuo / robo de gravedad estática | **No** |
| Batería gravitacional de inventario (lastres en altura) | **Sí** |
| Conmutador de densidad por pernos (3/4) | **Sí** |
| Sustituto de generador primario (viento, red, corriente) | **No** |
| Complemento tipo BESS / buffer / desacople de vehículo | **Sí** |
| Paridad impar genera energía | **No** |
| Recarga de distancia en un solo lado genera energía | **No** (solo fija la estación del objeto) |
| Electroimán / perneo barato (interruptor) | **Sí** |
| El interruptor borra \(m g h\) | **No** |

---

## 10. Analogías útiles

### 10.1 Estantería de ladrillos

1. Ladrillos en el desván (ALTA).  
2. Te cuelgas uno y bajas por la escalera generando.  
3. Dejas el ladrillo en el patio (BAJA).  
4. Tú vuelves arriba ligero.  
5. “Tú” estás como al inicio; **el ladrillo no**.  
6. Cuando no quedan ladrillos en el desván, se acabó.  
7. Subir los del patio es la recarga de la batería.

### 10.2 Noria de cubos

- Un lado baja “lleno” (denso).  
- El otro debería subir “lleno” o devolver el lastre.  
- Si no, solo has vaciado un depósito.

### 10.3 BESS químico

| BESS químico | Esta máquina |
|--------------|--------------|
| Carga con red | Subir lastres a cota alta |
| Descarga | Enganchar lastre, bajar, generar |
| SOC (estado de carga) | Número de pesos en stock ALTA |
| Pérdidas ida/vuelta | \(\eta_{\text{gen}}\), \(\eta_{\text{lift}}\), drag, pernos |

---

## 11. Scripts y artefactos en este directorio

| Archivo | Contenido |
|---------|-----------|
| `sim_enjambre_pesos_impar.py` | Sim principal 3/4, stock finito, con/sin lift |
| `sim_enjambre_pesos_impar.png` | Gráficas inventario y balances |
| `sim_enjambre_pesos_impar.json` | Resultados numéricos |
| `ciclo_iman_perneo_recorrido.py` | Balance analítico imán ↔ recorrido |
| `ciclo_parado_extiende_giro_recoge.py` | Alargar parado / recoger en giro |
| `recorrido_reduccionista.py` | Radio variable, parado vs centrifuga |
| `doble_km_90_perneo.py` | Doble módulo 90° + perneo de picos |
| `check_gemini_y_fisica.py` | Chequeo script Gemini patada / 1.ª ley |
| `MAQUINA_KILOMETRO_LASTRE.md` | Este documento |

### Cómo ejecutar la sim del lastre

```bash
python kilometro_sim/sim_enjambre_pesos_impar.py
```

---

## 12. Secuencia de diseño recomendada (si se prototipa)

1. **Estática de flotación**: calibrar 3 flota / 4 se hunde en tanque (masas reales).  
2. **Pernos make-before-break**: nunca 0/0; medir energía de conmutación.  
3. **Un solo Δh fijo** y un solo peso transferible (mínimo viable).  
4. **Medir** por separado:  
   - energía en bajada,  
   - coste de subir el lastre,  
   - pérdidas de arrastre.  
5. Solo después: enjambre, giro 180°, multi-estación, acoplamiento a red/picos.  
6. Métricas de éxito **honestas**:  
   - \(\eta_{\text{ida-vuelta}} = E_{\text{out}}/E_{\text{in}}\),  
   - energía por ciclo a SOC constante,  
   - no “julios mientras el almacén se vacía” sin contar el reset.

---

## 13. Conclusión

La máquina simulada es un **ascensor hidrostático de lastre con embrague de pernos**:

- el **fluido** desacopla el vehículo (sube barato en neutro),  
- los **pernos** conmutan densidad (3/4),  
- el **recorrido** es riel y estantería,  
- la **energía** es la PE de los kg en altura,  
- la **recarga** es volver a subir lastres (u otra fuente externa).

Es **viable y valiosa** como:

- batería gravitacional de inventario,  
- buffer y logística de masas en enjambre,  
- mecanismo de fase / perneo sin hidráulica pesada.

No es viable como:

- generador autónomo de gravedad,  
- sobreunidad por ciclos impares,  
- “robo de distancia” que borre \(m g h\).

---

*Documento generado a partir del diseño VMA y las simulaciones en `kilometro_sim/` (2026).*
