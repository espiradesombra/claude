# Explicación física unificada  
## ZypyZape · Quijote · Kilómetro · y el hilo **1477**

**Autor del corpus:** Víctor Manzanares Alberola (VMA) · EPSA / UPV Alcoi  
**Síntesis documental:** Grok (xAI) · 2026-07-23  
**Fuentes principales del monorepo:**

| Fuente | Ruta |
|--------|------|
| Paper 3 vs 7 + Kilómetro | `1/Paper_3vs7_Quijote_Kilometro.docx` |
| Paper bilingüe Quijote+ZZ | `Quijote/QUIJOTE_ZYPYZAPE_paper_llibre_bilingual.docx` |
| Contexto twin v4.8 | `VMA_mates_rescat_2026/10_inicio_quijote_zypyzape/01_CONTEXT_ZYPYZAPE.txt` |
| Maths 3 vs 7 pales | `…/02_MATH_QUIJOTE_3vs7.txt` |
| Kilómetro / simulación | `kilometre;(soles_bateria)/` · `claude.txt` |
| Hurto gravitatorio (código) | `hurto-gravitatorio/` |
| 1477 (fase booleana) | `33x1/1477.md` · `1/1477*.md` |
| 33×1 (trato) | `33x1/00_QUE_ES_33x1.md` |

**Uso:** civil · parte del **“1”** del marco **33×1**.  
**Honestidad:** esta nota distingue **física sólida (conservación, inercia, red)** de **hipótesis aún no demostradas en prototipo**.

---

## 0. Mapa mental en una página

```
                    ┌─────────────────────┐
                    │   FASE / ASIMETRÍA  │  ← idea abstracta (1477 + control)
                    └──────────┬──────────┘
           ┌───────────────────┼───────────────────┐
           ▼                   ▼                   ▼
    ┌────────────┐      ┌────────────┐      ┌────────────┐
    │  QUIJOTE   │      │  KILÓMETRO │      │  ZYPYZAPE  │
    │ masa en    │      │ kilo+metre │      │ red de     │
    │ pala: J(r) │      │ trayecto + │      │ máquinas   │
    │ variable   │      │ fricción / │      │ acopladas  │
    └─────┬──────┘      │ flotación  │      └─────┬──────┘
          │             └─────┬──────┘            │
          └──────────┬────────┘                   │
                     ▼                            │
              INERCIA LOCAL                       │
              (buffer en el rotor)                │
                     └────────────┬───────────────┘
                                  ▼
                         INERCIA DE RED
                         (H_eff, RoCoF, nadir)
```

| Pieza | Qué es, en físico | Qué **no** es |
|-------|-------------------|---------------|
| **Quijote** | Control de **inercia del rotor** moviendo masa radial en la pala | Una batería de litio |
| **ZypyZape** | **Intercambio cíclico** de energía entre turbinas + topología (Kuramoto / anillo) | Solo software de marketing |
| **Kilómetro** | Máquina **anclada** (kilo + metro de recorrido) que aplica la misma lógica de asimetría / buffer | Un aerogenerador |
| **1477** | Método **booleano de fase** (bandera de coherencia / alineación de fase) | Hardware eólico |

Hilo conductor:

> **La gravedad y la inercia son campos / estados que ya existen.**  
> El sistema no “inventa” energía: **orquesta cuándo** la masa es pesada/ligera, **cuándo** se frena/acelera, y **cómo** se reparte el momento entre máquinas.

---

## 1. El problema que motivan las tres máquinas

Las redes con mucha eólica/solar **pierden inercia síncrona** clásica (menos grandes generadores rotativos). Tras un fallo de potencia:

- sube el **RoCoF** (|df/dt| grande),
- baja el **nadir** de frecuencia,
- se necesitan BESS, síncronos o grid-forming caros.

**Hipótesis VMA:** en lugar de solo añadir baterías químicas, **usar la masa que ya gira** (rotor eólico, o una máquina tipo Kilómetro) como **buffer de inercia variable y repartida**.

---

## 2. Física compartida: inercia variable y “efecto patinadora”

### 2.1 Inercia constante → simetría

Con \(J\) fijo, en un ciclo cerrado de altura:

\[
W_{\mathrm{bajada}} = m g \Delta h = W_{\mathrm{subida}}
\]

Campo gravitatorio **conservativo** → en un ciclo cerrado de posición, el trabajo neto de la gravedad es **cero** (1ª ley / potencial).

### 2.2 Inercia variable \(J(t)\)

Si se mueve masa radialmente:

\[
J(r) = J_G + N_b \, m_q \, r^2
\]

La ecuación de movimiento del rotor (forma usada en los gemelos):

\[
J(t)\,\frac{d\omega}{dt} + \omega\,\frac{dJ}{dt} = T_{\mathrm{aero}} + T_{\mathrm{otros}} - T_{\mathrm{gen}}
\]

- Término \(\omega\,dJ/dt\): **efecto patinadora**  
  - \(dJ/dt < 0\) (masas hacia el eje) → empujón a \(\omega\)  
  - \(dJ/dt > 0\) (masas hacia la punta) → frena \(\omega\) (carga “batería” cinética)

Energía cinética del rotor:

\[
E = \tfrac12 J(t)\,\omega^2
\]

Al cambiar \(J\) y \(\omega\) **a la vez**, el reparto entre “almacenar en \(\omega\)” y “almacenar en \(J\)” se controla.

### 2.3 Lectura honesta del “hurto gravitatorio”

En el paper `Paper_3vs7_Quijote_Kilometro` se define un trabajo por vuelta ligado al peso y al ángulo:

\[
W_{\mathrm{hurto},k} = m_q\, g\, \Delta r \oint |\cos\theta|\,d\theta = 4\, m_q\, g\, \Delta r
\]

\[
P_{\mathrm{hurto}} = N\, f_{\mathrm{rot}}\, 4\, m_q\, g\, \Delta r
\]

**Importante (auditoría del propio corpus):**

1. Esa integral **describe un término de acoplamiento** peso–ángulo–inercia, no un permiso para violar la 1ª ley.  
2. El **actuador** que mueve la masa (fluido, motor, fricción controlada) **consume trabajo**. En análisis posteriores (`kilometre/…/claude.txt`, paper de inercia variable) se insiste en que **el coste de actuar puede superar el “hurto”** si no se diseña bien.  
3. El uso **sólido y vendible hoy como idea de ingeniería** no es “motor perpetuo de gravedad”, sino:

   - **buffer de inercia** (PFR/FFR, RoCoF),  
   - **reparto de energía entre máquinas** (ZypyZape),  
   - **control de fricción / flotación** (Kilómetro submarino como batería mecánica).

> Frase del paper: *“La energía la tenía el viento, la robé y la dejé ir.”*  
> Lectura física segura: el **viento** (o una fuente externa) aporta la energía; Quijote/ZZ **reordenan el timing** de la inercia.

---

## 3. QUIJOTE — “molino dentro del molino”

### 3.1 Qué es

Masa (o fluido denso Fe+oli) que **desliza radialmente** en cada pala:

- \(r \in [r_{\min}, r_{\max}]\) (p.ej. 5 m → 55–60 m),  
- control de velocidad de deslizamiento \(v_{\mathrm{slide}}\).

**Metáfora:** un molino pequeño (el peso) dentro del molino grande (la turbina).

### 3.2 Dos modos de control

| Modo | Idea | Uso de red |
|------|------|------------|
| **Estático** | Dejar masa fuera o dentro y mantener | “Batería” lenta (PFR, 10–30 s) |
| **Ball** (sinusoidal entre palas) | \(r_k = r_0 + A\sin(\omega t + 2\pi k/N)\) | “Condensador” rápido (FFR &lt; 2 s) |

Propiedad matemática (3 o 7 palas equiespaciadas, ball en fase):

\[
J_{\mathrm{total}} = m N \big(r_0^2 + A^2/2\big) = \textbf{constante en el tiempo}
\]

Misma estructura que la potencia trifásica: las sumas de senos se anulan.

Regla de control a \(J\) constante: si una pala extiende \(A\), las otras \(N-1\) retraen \(A/(N-1)\).

### 3.3 3 palas vs 7 palas

| | N=3 | N=7 |
|--|-----|-----|
| Separación angular | 120° (huecos de hurto) | ~51.4° (más continuo) |
| \(\Delta J_{\max}\) relativo | 1 | **2.33×** |
| Rizado de par | alto | menor |

Parámetros típicos del twin (orden de magnitud, no “potencia gratis”):

- NREL-class: \(R\sim 60\,\mathrm{m}\), \(J_G\sim 10^6{-}10^7\,\mathrm{kg\,m^2}\),  
- \(\Delta E\) Quijote por ciclo ~ decenas–cientos de kJ por turbina (local),  
- A red de **2 GW**, \(\Delta H\) de un solo Quijote es **casi invisible**; brilla en **microred / isla / 200–500 MW**.

### 3.4 Ecuaciones de control (gemelo)

\[
v_{\mathrm{slide}} = \mathrm{clip}\big(K_\omega\Delta\omega + K_f\Delta f,\, -V_{\max},\, V_{\max}\big)
\]

- Exceso de \(\omega\) o \(f\) → masas **afuera** (cargar buffer),  
- Déficit → masas **adentro** (descargar).

---

## 4. ZYPYZAPE — batería cinética en red

### 4.1 Qué es

Sistema de **varias turbinas** que se intercambian energía de forma **cíclica**:

- una **acelera** (acumula \(\tfrac12 J\omega^2\)),  
- otra **frena** (inyecta a red / cede energía),  
- periodo típico \(T\sim 2.5\,\mathrm{s}\) (\(f\sim 0.4\,\mathrm{Hz}\)),  
- fracción de potencia intercambiada ~10–15 % de \(S_{\mathrm{nom}}\) por máquina.

Topologías del corpus:

- **5 nodos:** 2 centrales + anillo de 3 (documento técnico clásico),  
- **10 nodos:** 2 pares + 2 anillos de 4 (twin v4.8).

### 4.2 Inercia de red (swing equation)

\[
H_{\mathrm{eff}} = H_{\mathrm{sys}} + H_{\mathrm{ZZ}},\qquad
H_{\mathrm{ZZ}} = \frac{\tfrac12 \sum_i J_i\omega_i^2}{S_{\mathrm{total}}}
\]

\[
\frac{df}{dt} = \frac{( \Delta P + P_{\mathrm{gov}})\,f_0}{2\, H_{\mathrm{eff}}\, S_{\mathrm{total}}} - D(f-f_0)
\]

**Efecto realista medido en gemelos:** mejora de nadir / RoCoF **pequeña por módulo**, **lineal** si se suman muchos módulos.  
Eso es compatible con física de red; no requiere \(\eta>1\).

### 4.3 Kuramoto (sincronización)

Modelo de acoplamiento (forma conceptual del corpus 33×1):

\[
\dot\theta_i = \omega_0 + K\sum_j \sin(\theta_j-\theta_i) + u_i(t)
\]

\(u_i\): entrada de control Quijote / consignas ZZ.  
La red de máquinas se **sincroniza** sin un “cerebro” único: acoplamiento mecánico-eléctrico.

### 4.4 Qué resuelve ZypyZape (si se valida en piloto)

| Servicio | Principio |
|----------|----------|
| Inercia sintética / RoCoF | \(E_{\mathrm{cin}}\) real del rotor |
| Regulación de frecuencia | droop + ciclo ZZ |
| Reuso de hardware | firmware + topología; poco capex nuevo |
| Menos dependencia de BESS para ciertos ancilares | buffer mecánico repartido |

---

## 5. KILÓMETRO — “kilo” + “metro”

### 5.1 Nombre y idea

**Kilómetro** = **kilo** (masa) + **metro** (recorrido).  
No es un aerogenerador: es una **máquina anclada** (a menudo pensada **submarina** o con fluido) donde una masa recorre un trayecto y se controla fricción / flotación / giro.

Ideas del corpus (`claude.txt`, paper, sims `kilometre_v*`):

1. **Trayecto recto** que puede **girar sobre un eje**, con **sinfín** interior.  
2. El objeto avanza y arrastra el sinfín → **fricción modulable** = donde se “saca” o se frena energía.  
3. **“Vuelta y media”:** desfase entre el giro del trayecto y el avance de la masa para que el **par favorable** y el **desfavorable** no sean simétricos en el tiempo de control (no en el potencial gravitatorio).  
4. **Flotación / Arquímedes** bajo el agua: cambia el peso efectivo y la asimetría de subida/bajada.  
5. Variantes Stirling / flota de unidades pequeñas (`kilometre_stirling_*`, `kilometre_v14_flota`).

### 5.2 Formulación “vuelta y media” (versión 33×1 ciencia)

En `33x1_CIENCIA/00_PROMPT_DEFINITIVO_33x1.md`:

- Fase bajada (rápida): máxima extracción,  
- Fase subida (lenta, inercia del carril): mínima reinversión,  
- Ciclo 1.5 vueltas: reset con inercia angular acumulada.

**Auditoría física (misma carpeta de chat `claude.txt`):**

- Si el estado del sistema **se cierra periódicamente** en el espacio de configuración (misma altura, misma \(\omega\), misma posición relativa), el trabajo neto de un campo **conservativo** (gravedad) es **cero**.  
- “2/3 del tiempo a favor” **no** crea energía si en el 1/3 desfavorable el par es el doble.  
- Lo que **sí** se puede optimizar:  
  - **cuándo** se empuja (a baja \(\omega\), menos trabajo para el mismo \(\Delta\omega\)),  
  - **fricción** (regenerar en freno, reducir en avance),  
  - **fuente externa** (eólica, solar marina, red) que **carga** la batería mecánica.

### 5.3 Lectura útil del Kilómetro (ingeniería)

| Enfoque débil | Enfoque sólido |
|---------------|----------------|
| “Hurtamos gravedad y generamos sin entrada” | “Batería submarina / volante con geometría y fricción controlada” |
| Ciclo cerrado sin pérdidas y sin actuador | Entrada externa + \(\eta<1\) + servicios de red |
| Un solo Kilómetro “autosuficiente” | **Grupo** acoplado vía **ZypyZape** (paper: 3+ unidades) |

El paper 3vs7 dice: *el caso 3 vs 7 palas es la prueba de escalado del principio hacia el Kilómetro* (más ángulos → más continuidad del “ball”). Eso es un **argumento de transferencia de principio**, no un certificado de prototipo industrial.

---

## 6. Cómo encajan las tres (+ 1477)

### 6.1 Cadena física

```
Viento / red / flotación  →  energía de entrada
        │
        ▼
   QUIJOTE (J local variable en pala o masa)
        │
        ▼
   ZYPYZAPE (reparte E_cin entre máquinas, sincroniza)
        │
        ▼
   Red: H_eff ↑, RoCoF ↓, servicios ancilares
        │
        ▼
   KILÓMETRO (misma lógica en máquina anclada / submarina)
```

### 6.2 Hilo **1477** (fase abstracta)

`33x1/1477.md` describe un método **booleano**:

- qubit → bits (valor, fase),  
- desfase tipo Karnaugh,  
- realimentación AND clásico vs Toffoli,  
- **bandera** de “hay coherencia de fase” si la distancia de Hamming crece.

**Puente conceptual (no identidad física):**

| 1477 (bits) | Energía VMA |
|-------------|-------------|
| Fase \(f\) | Fase de control / ángulo de rotor / rol gen-thr-buf |
| Desfase cíclico | Ciclo ZZ 0.4 Hz · ball Quijote · vuelta y media |
| Realimentación de diferencia | Control \(v_{\mathrm{slide}}\), droop, acoplamiento Kuramoto |
| Bandera de coherencia | Sincronización de flota / alineación de red |

1477 **no mueve aspas**; es la capa **lógica de fase** del ecosistema 33×1, empaquetada junto a ZypyZape / Quijote / Kilómetro en la oferta del **1**.

### 6.3 XFI (hermano menor)

El gemelo **XFI** (avión N=3: gen/thr/buf) es la **misma idea de roles cíclicos** aplicada a un craft:  
captura en un motor, empuje en otro, buffer en el tercero — ver `XFI.md`.

---

## 7. Qué está en suelo firme vs qué es hipótesis

| Afirmación | Estado en el corpus |
|------------|---------------------|
| \(J(t)\) variable cambia la dinámica del rotor (patinadora) | **Firme** (mecánica clásica) |
| Buffer cinético mejora servicios de frecuencia **en simulación** | **Plausible / simulado** (gemelos NREL-class) |
| Topología multi-máquina reduce RoCoF | **Plausible** (swing + acoplamiento) |
| Coste de actuador Quijote puede matar el “hurto” neto | **Admitido** en auditorías internas |
| Kilómetro genera energía neta sin fuente externa | **No establecido**; reencuadrar como batería |
| Grupos ZZ/Kilómetro “autosuficientes” | **Especulativo** hasta piloto |
| Sustituir BESS a gran escala | **Mercado + validación OEM** pendientes |
| 1477 detecta “fase cuántica” en hardware real | **Método conceptual / flag lógico** |

---

## 8. Ecuaciones mínimas para recordar

**Quijote**

\[
J = J_G + N_b m_q r^2,\quad
J\dot\omega + \omega\dot J = T_{\mathrm{net}}
\]

**ZypyZape (red)**

\[
H_{\mathrm{eff}} = H_{\mathrm{sys}} + \frac{\sum \tfrac12 J_i\omega_i^2}{S_{\mathrm{base}}}
\]

**Kilómetro (balance honesto)**

\[
W_{\mathrm{externo}} = W_{\mathrm{util}} + W_{\mathrm{perdidas\ no\ recuperadas}}
\qquad (\text{siempre } W_{\mathrm{externo}} \ge 0 \text{ en régimen estacionario cerrado})
\]

**Meta de ingeniería**

\[
\min W_{\mathrm{externo}}
\quad\text{s.a.}\quad
\text{servicio de red / potencia útil exigida}
\]

---

## 8b. Gemelo mínimo Kilómetro + diagrama (2026-07-23)

| Recurso | Ruta |
|---------|------|
| Diagrama 1 página | `docs/diagrama_fisica_ZZ_Quijote_Kilometro.png` |
| Script diagrama | `docs/generate_diagrama_fisica_1pagina.py` |
| Gemelo balance cerrado | `docs/gemelo_kilometro_minimo.py` |
| Gráfica balance | `docs/kilometro_minimo_balance.png` |

```bat
cd docs
python generate_diagrama_fisica_1pagina.py
python gemelo_kilometro_minimo.py --mode drain --T 25
```

Identidad usada en el gemelo:

\[
\Delta E_{\mathrm{mec}} = W_{\mathrm{fric}} + W_{\mathrm{ext}}
\qquad
(\,W_g \text{ solo intercambia } E_{\mathrm{cin}}\leftrightarrow E_{\mathrm{pot}}\,)
\]

\[
\eta_{\mathrm{paid}} = \frac{E_{\mathrm{util}}}{E_{\mathrm{ext}}^+ + \max(0,-\Delta E_{\mathrm{mec}})} \le \eta_{\mathrm{reg}} < 1
\]

---

## 9. Dónde está el código y los papeles

| Tema | Dónde |
|------|--------|
| Twin ZZ + Quijote | `VMA_mates_rescat_2026/09_*`, `10_*`, `ZYPYZAPE Bateria Cinetica/`, `Quijote/` |
| XFI | `XFI.md`, `…/10_inicio_quijote_zypyzape/` |
| Kilómetro sims | `kilometre;(soles_bateria)/kilometre_v*.py` |
| Hurto / 3vs7 | `hurto-gravitatorio/`, `1/Paper_3vs7_Quijote_Kilometro.docx` |
| 1477 | `33x1/1477.md`, `1/1477*.md`, `1/quantum_1477_*.py` |
| Trato 33×1 | `33x1/00_QUE_ES_33x1.md` |

---

## 10. Conclusión en tres frases

1. **Quijote** cambia **cuándo** el rotor es “pesado” o “ligero” (inercia variable en la pala).  
2. **ZypyZape** hace que **varias** máquinas se presten inercia y se sincronicen (batería cinética de parque).  
3. **Kilómetro** lleva el mismo principio a una **máquina de trayecto** (a menudo submarina); el marco útil es **batería mecánica + fricción/flotación controlada**, no un motor de gravedad libre.  
4. **1477** es el lenguaje de **fase** del pack 33×1: alinear diferencias, no saltarse la termodinámica.

---

*Documento de síntesis. No sustituye un paper revisado por pares ni un ensayo en banco.  
Cualquier demo comercial exige prototipo, fatiga, grid code y medida real.*
