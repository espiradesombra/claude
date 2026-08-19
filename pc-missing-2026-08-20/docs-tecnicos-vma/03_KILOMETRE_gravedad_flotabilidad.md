# Kilòmetre / Kilómetro  
## Máquina de gravedad, flotabilidad y lastre — Documento técnico (v1.1)

**Versión:** 1.1 (incorporando `kilometro_sim`, 2026-08)  
**Dominio:** almacenamiento mecánico / buffer de potencia / conmutación de densidad  
**TRL estimado:** 3–4 (simulación extensa con cierre de 1ª ley + diseño de prototipo tanque/ESP32; sin prototipo industrial medido en campo)

---

## 1. Resumen ejecutivo técnico (honesto)

**Kilómetro** es una familia de mecanismos VMA cuyo nombre evoca el módulo elemental:

> desplazar del orden de **1 kg · 1 m** de forma controlada para gestionar energía y par.

Tras el corpus `repo/kilometre;(soles_bateria)/` (v2–v15) y el paquete **`kilometro_sim/`** (chequeos 2026-08-09), el veredicto unificado es:

| Rol | ¿Válido? |
|-----|----------|
| Motor perpetuo / extracción neta de gravedad en ciclo cerrado | **No — refutado** |
| Batería gravitacional de **inventario de lastres** | **Sí** (contabilidad de almacén) |
| Buffer de **picos de potencia** (perneo / doble módulo) | **Sí** (simulación) |
| Convertidor con \(\eta_{paid}<1\) | **Sí** |
| Sustituto de generador primario (viento, nuclear, red) | **No** |

> **Frase de diseño:** el perneo devuelve el *objeto* al traje inicial; **no** devuelve el *kg* a la cota alta.

---

## 2. Tres arquitecturas (estado del arte del corpus)

### 2.A — Tubo / recorrido + sinfín + rotación (línea clásica v2–v15)

Tubo con objeto y sinfín; posible rotación del recorrido; control de fricción/generador; variantes de flotación continua y fluidos densos (salmuera, Fe+aceite, Ga, Hg de referencia).

**Estado:** muchas sims de barrido; el propio `kilometre_v2.py` marca balances positivos como **sospechosos** si \(W_g\not\approx 0\).

### 2.B — ForPkm multi-módulo (2–3 kilómetros en equipo)

Asistencia cruzada en ventanas de baja inercia; 3 módulos a 120° para suavizar par (analogía trifásica).  
**Valor:** redistribución temporal de energía y reducción de picos — no \(W_{net}>0\).

### 2.C — Ascensor de lastre hidrostático + pernos (**línea `kilometro_sim`**, canónica actual)

Documento maestro: `kilometro_sim/MAQUINA_KILOMETRO_LASTRE.md`.

#### Piezas

| Pieza | Función |
|-------|---------|
| **Objeto** (módulo Kilómetro) | Viaja por la guía; con **3 pesos** flota; con **4** se hunde |
| **Recorrido** | Riel + estantería de lastres (stock ALTA / BAJA) |
| **Pesos** | Lastres discretos = “celdas” de la batería |
| **Pernos R / O** | Embrague: peso↔Recorrido o peso↔Objeto (make-before-break) |
| **Fluido** | Boyancia barata para subir el módulo ligero |

#### Umbral 3 flota / 4 se hunde

Un solo lastre conmuta el comportamiento hidrostático — sin bombas grandes de volumen.

Parámetros de simulación de referencia (`sim_enjambre_pesos_impar.py`):

| Símbolo | Valor |
|---------|-------|
| \(m_{peso}\) | 10 kg |
| \(m_{obj,base}\) | 30 kg |
| \(n_{float}/n_{sink}\) | 3 / 4 |
| \(\Delta h\) | 15 m |
| \(N_{stock,alta}\) | 5 |
| \(\eta_{gen}/\eta_{lift}\) | 0,85 / 0,90 |
| \(E_{perno}\) | 1,5 J |
| \(f_{drag}\) | 0,06 |

#### Ciclo modo A (almacén abierto — descarga de batería)

```text
ALTA: objeto n=3 flota, stock listo
  → enganche (perneo) n=4 se hunde
  → bajada Δh + regeneración (solo lastre EXTRA cuenta como PE gastable)
  → entrega en BAJA n=3
  → subida por boyancia + recarga de distancia en estación S
  → stock ALTA −1, BAJA +1
```

Hasta `stock_ALTA=0` → **STOP**.  
Número de ciclos = tamaño del almacén (**no** magia de paridad impar/par).

#### Ciclo modo B (cierre con lift)

Subir cada lastre BAJA→ALTA cuesta más de lo recuperado:

\[
W_{gen}\approx\eta_{gen}\,m g\Delta h,\qquad
W_{lift}=\frac{m g\Delta h}{\eta_{lift}}
\quad\Rightarrow\quad
W_{lift}>W_{gen}\ \text{si }\eta<1
\]

---

## 3. Resultados numéricos de `kilometro_sim` (reproducibles)

### 3.1 Chequeo Gemini + 1ª ley (`check_gemini_y_fisica.py`)

| Modelo | Resultado |
|--------|-----------|
| Gemini **original** (params del chat) | **No completa** el ciclo; se queda en bajada; \(E_{neto}\sim+53\) J es **falso éxito** (PE parcial, sin rearme) |
| Gemini **parcheado** | Completa fases; default \(E_{neto}\sim-466\) J; positivos del sweep = **PE no rearmada** |
| **Ciclo cerrado** 5 vueltas | \(W_{motor}\approx4391\) J, \(W_{gen}\approx3491\) J, \(W_{net}\approx-900\) J, \(\eta_{paid}\approx0{,}80\), \(W_g\sim0\) |
| Sweep 81 combos, 48 cierres válidos | **0** con surplus post-KE &gt; 0; mejor \(\eta_{paid}\sim0{,}84\) |

**Veredicto del propio JSON:**

> *Kilometro = batería/convertidor, no motor perpetuo. Doble kilometro por perneo = buffer de potencia, no generación neta.*

### 3.2 Enjambre lastre 3/4 (`sim_enjambre_pesos_impar`)

| Concepto | Valor aprox. |
|----------|----------------|
| \(m g h\) por peso (10 kg, 15 m) | ~1472 J |
| \(W_{rec}\) honesto / ciclo | ~1176 J |
| Pernos / ciclo | ~6 J |
| Surplus / ciclo **sin** lift | ~1170 J |
| Ciclos hasta STOP (N=5) | **5** |
| Surplus acum. A sin reset | ~5849 J |
| Coste subir 5 pesos | ~8175 J |
| Saldo al **cerrar** almacén | **~−2326 J** |
| Surplus / ciclo **con** lift (modo B) | **~−465 J** |

### 3.3 Doble Kilómetro 90° + perneo (`doble_km_90_perneo.py`)

**Política de regeneración (explícita):**

| Fase | \(\cos\phi\) | Actuación |
|------|--------------|-----------|
| Desfavorable | \(\ge 0\) | Motor o libre — **REGEN OFF** |
| Favorable | \(&lt; 0\) | Regeneración ON |

Con dos masas a 90°:

\[
\tau_{g,total} = -m g R\sqrt{2}\,\cos(\phi+\pi/4)
\]

→ el peak-to-peak gravitatorio **crece \(\times\sqrt{2}\)**; la utilidad es **desfase** y embrague, no cancelar el 1.er armónico (eso es 3×120°).

Métricas de un caso de referencia (pulso 2,5 kW × 0,8 s ≈ 2 MJ pedidos):

| Modo | \(W_{net,elec}\) | Cobertura de pulso | Nota |
|------|------------------|--------------------|------|
| Mono | ~−1524 J | **0 %** | \(\eta_{paid}\sim0{,}58\) |
| Dual rígido | ~−1952 J | ~10 % | mejor que mono en entrega de pulso |
| Dual + perneo | ~−2682 J | **~43 %** | supercondensador mecánico ante picos |

\(\eta_{paid}<1\) siempre. Objetivo de diseño: **\(P_{pico}/P_{media}\)** y cobertura de pulsos, **no** \(W_{net}>0\).

### 3.4 Recorrido reduccionista (`recorrido_reduccionista.py`)

- \(r\) pequeño en tramos desfavorables, grande en favorables.  
- Rearme de radio **en parado** (sin centrífuga) cuesta menos que en giro.  
- Ciclo reduccionista: **gasto eficiente del eje**, no sobreunidad.  
- Gravedad puede ayudar en parado según ángulo (p. ej. retract en top / extend en bottom).

---

## 4. Modelo energético formal (sin trampa)

### 4.1 Campo gravitatorio

\[
\oint \mathbf{F}_g\cdot d\mathbf{r} = 0
\]

en ciclo cerrado de **cada** masa.

### 4.2 Modo almacén abierto

\[
E_{util}\approx N\cdot\eta_{gen}\cdot(1-f_{drag})\cdot m_{peso}\,g\,\Delta h - E_{pernos}
\]

= descarga de PE que **alguien cargó antes** (red, grúa, oleaje, otro proceso).

### 4.3 Cierre

\[
E_{reset}=N\cdot\frac{m_{peso}\,g\,\Delta h}{\eta_{lift}}
\quad\Rightarrow\quad
E_{util}-E_{reset}<0
\]

típico con \(\eta<1\).

### 4.4 Analogía BESS

| BESS químico | Kilómetro lastre |
|--------------|------------------|
| Carga con red | Subir lastres a ALTA |
| Descarga | Enganchar, bajar, generar |
| SOC | Nº de pesos en stock ALTA |
| Round-trip | \(\eta_{gen}\times\eta_{lift}\times(1-drag)\) |

---

## 5. Prototipo físico en curso

### 5.1 Tanque MVP

`PROTOTIPO_TANQUE_PLANOS_Y_BOM.md` + `BOM_prototipo_tanque_v1.csv`  
Calibrar **3 flota / 4 hunde** con masas reales en agua.

### 5.2 v1.1 ESP32 + solenoides

`kilometro_sim/v1_1_ESP32/`:

- FSM: `IDLE_ALTA → ENGANCHE → ESPERA_BAJADA → ENTREGA → ESPERA_SUBIDA → CICLO_OK`  
- Make-before-break en firmware  
- Finales de carrera ALTA/BAJA  
- Log Serial de tiempos  
- Recarga de lastre a ALTA puede seguir **manual** en v1.1  

**Sigue siendo batería de lastre**, no generador perpetuo.

### 5.3 Criterios de aceptación de prototipo

1. Balance medido: \(E_{out}/E_{in}\) con SOC (stock) **constante** (incluye lift).  
2. Nunca publicar solo “julios mientras se vacía el desván” sin reset.  
3. Cierre residual mecánico del modelo &lt; 1 % en sim antes de hardware.  
4. Interlocks: prohibido estado perno 0/0 (peso suelto).

---

## 6. Fluidos densos (línea v15 — diseño de materiales)

| Fluido | \(\rho\) | Uso |
|--------|----------|-----|
| Salmuera 25 % | ~1200 | Referencia barata |
| Fe+aceite 30–60 % Fe | ~2,7–5,1×10³ | Práctico; sedimentación |
| Galio | ~6095 | Costoso |
| Mercurio | 13600 | Solo referencia; **tóxico** |

Objetos de contraste: He / vacío.  
Fuerza de orden \(F_n=(\rho_f-\rho_o)V g\).

---

## 7. Interfaz con ZypyZape y Quijote

```text
Eólica / red / grúa  ──carga──►  Kilómetro (SOC lastre o KE)
                                  │ descarga / picos
                                  ▼
                         servicios de red / microred
                                  │
                    fase opcional → ZypyZape (ACEL/FREN)
Quijote = en pala; Kilómetro = módulo externo. Misma filosofía: orquestar inercia.
```

En un **emplazamiento nuclear apagado** (ver doc 05): bucles de fluido y máquinas rotativas pueden hostear ideas Kilómetro/ZZ como **campus de buffer**, no como magia residual.

---

## 8. Lo que un interesado puede pedir (y obtener)

| Entregable | Estado |
|------------|--------|
| Scripts + JSON de 1ª ley | **Disponibles** en `kilometro_sim` |
| Doc de máquina de lastre | **Disponible** |
| BOM + firmware ESP32 v1.1 | **Diseño** |
| Medida de η en tanque real | **Pendiente** |
| Planta submarina multi-MW | **Conceptual** — no TRL de obra |

---

## 9. Conclusiones técnicas

1. Kilómetro **no** genera energía neta de la gravedad.  
2. Kilómetro **sí** es una batería de lastre / buffer de potencia con contabilidad clara.  
3. `kilometro_sim` eleva el rigor: cierra la 1ª ley y desmonta falsos positivos de chats.  
4. El camino creíble es: **tanque → η medida → enjambre → acoplamiento a red**, no marketing de sobreunidad.  
5. En transición energética (doc 05): aporta **servicios y almacenamiento**, no baseload nuclear.

---

## 10. Referencias de archivo

| Artefacto | Ruta |
|-----------|------|
| Chequeo + veredicto | `kilometro_sim\RESUMEN_CHEQUEO.txt`, `resultados_chequeo.json` |
| Lastre completo | `kilometro_sim\MAQUINA_KILOMETRO_LASTRE.md` |
| Doble 90° | `kilometro_sim\README_DOBLE_KM.md` |
| ESP32 | `kilometro_sim\v1_1_ESP32\` |
| Sims históricas | `repo\kilometre;(soles_bateria)\` o `33x1\repo\...` |
| Veredicto archivo ingeniería | `33x1\repo\teoremasgrok\ingenieria\05_KILOMETRO.txt` |

---

*Documento 03/05 — Dossier técnico VMA v1.1. Solo contenido técnico.*
