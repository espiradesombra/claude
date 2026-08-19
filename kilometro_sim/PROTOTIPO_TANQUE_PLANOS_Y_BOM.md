# Prototipo Kilómetro — Tanque (planos lógicos + BOM)

**Objetivo:** validar en físico, a escala de taller, la **batería de reset por lastre de perneado**.

**Alcance mínimo (MVP):**

| Elemento | Decisión MVP |
|----------|----------------|
| Desnivel | **1 solo Δh** (fijo) |
| Pesos transferibles | **1 lastre** extra (el 4.º) |
| Umbral | **3 flota / 4 se hunde** |
| Pernos | **2 por lastre** (R = recorrido, O = objeto) |
| Lift de recarga | **Manual** en v1 (no automático) |
| Generación eléctrica | **Opcional** (v1.1); v1 mide solo fuerzas/tiempos/energía mecánica estimada |

**No es:** generador perpetuo, planta industrial, ni “robo de gravedad”.  
**Sí es:** prueba de conmutación de densidad + ciclo de inventario de 1 lastre.

Documentos relacionados:

- `MAQUINA_KILOMETRO_LASTRE.md` — teoría de la máquina  
- `sim_enjambre_pesos_impar.py` — simulación de inventario  

---

## 1. Especificación funcional del MVP

### 1.1 Ciclo de prueba (un lastre)

```text
1. Objeto en cota ALTA, n=3, flota (o flota muy positiva).
2. Enganche make-before-break: lastre ALTA → objeto (n=4).
3. Objeto se hunde y baja guiado Δh.
4. En cota BAJA: make-before-break lastre → recorrido BAJA (n=3).
5. Objeto flota y sube (o se recupera con guía de baja fricción).
6. Operador (v1) devuelve el lastre a cota ALTA a mano → recarga batería.
```

### 1.2 Criterios de éxito (pass/fail)

| # | Criterio | Pass |
|---|----------|------|
| C1 | Con 3 pesos (o masa equivalente) el módulo **no se hunde** en reposo | Sí/No |
| C2 | Con 4 (tras enganchar 1 lastre) el módulo **se hunde** sin ayuda | Sí/No |
| C3 | Transferencia ALTA: make-before-break sin soltar el lastre al fondo del tanque | 0 caídas en 20 ciclos |
| C4 | Transferencia BAJA: igual | 0 caídas en 20 ciclos |
| C5 | Tiempo de conmutación pernos (eléctrico o manual) documentado | medido |
| C6 | Energía estimada bajada vs coste de devolver lastre a mano | tabla |

### 1.3 Parámetros de diseño recomendados (taller)

Escalables; empiezan baratos y medibles.

| Parámetro | Símbolo | Valor propuesto v1 | Notas |
|-----------|---------|---------------------|--------|
| Desnivel útil | Δh | **0,80 m** | tanque ~1,0–1,2 m alto |
| Sección tanque | — | 0,40 × 0,40 m (o Ø 0,45 m) | ver módulo |
| Masa estructural módulo | m_base | 1,5–2,5 kg | plástico / aluminio |
| Masa “3 pesos” total en flotación | m_3 | calibrar | ver §3 |
| Masa lastre transferible | m_L | **0,5–1,0 kg** | acero inox / latón |
| Margen 3→4 | — | +8–15 % densidad aparente | crítico |
| Guía | — | 2 varillas Ø8–10 mm | acero inox o aluminio anodizado |
| Fluido | — | agua dulce 15–25 °C | salmuera opcional v2 |

**Energía potencial de referencia del lastre:**

\[
E = m_L \cdot g \cdot \Delta h
\]

Ejemplo: \(m_L = 0{,}8\,\mathrm{kg}\), \(\Delta h = 0{,}8\,\mathrm{m}\):

\[
E \approx 0{,}8 \times 9{,}81 \times 0{,}8 \approx 6{,}3\,\mathrm{J}
\]

Escala de juguete, pero **suficiente para demostrar el principio** y medir con báscula + cronómetro + (opcional) célula de carga.

---

## 2. Arquitectura física

### 2.1 Vista general (plano lógico)

```text
                    ┌──────────────────────────────┐
                    │     TAPA / PUENTE SUPERIOR   │
                    │  [estación S: enganche ALTA] │
                    │         anclaje perno R_alta │
                    └────────────┬─────────────────┘
                                 │
         nivel agua ~~~~~~~~~~~~ │ ~~~~~~~~~~~~~~~~
                                 │
                    ┌────────────┼─────────────────┐
                    │            │                 │
                    │   VARILLA  │  VARILLA        │  PARED TANQUE
                    │     G1     │    G2           │
                    │            │                 │
                    │      ┌─────▼─────┐           │
                    │      │  MÓDULO   │           │
                    │      │  (objeto) │           │
                    │      │  casco +  │           │
                    │      │  2 puntas │           │
                    │      │  perno O  │           │
                    │      └─────┬─────┘           │
                    │            │                 │
                    │   [estación B: perneo BAJA]  │
                    │         anclaje perno R_baja │
                    └────────────┴─────────────────┘
                    │         FONDO REFORZADO      │
                    └──────────────────────────────┘
```

### 2.2 Módulo (objeto) — despiece lógico

```text
        punta superior (cara de enganche ALTA)
        ┌─────────────────────┐
        │  asiento lastre L   │  ← cavidad / horquilla
        │  [perno O] ────●    │  ← solenoide o pasador manual
        ├─────────────────────┤
        │  casco estanco      │
        │  lastre fijo (calibración 3) │
        │  opcional: logger / LED     │
        ├─────────────────────┤
        │  buje / cojinetes   │  ← deslizan en G1, G2
        │  punta inferior     │  ← cara de entrega BAJA
        └─────────────────────┘
```

**Funciones del módulo:**

1. Casco con **volumen** V tal que con masa m_3 flote y con m_3+m_L se hunda.  
2. **Guía** anti-vuelco (dos varillas).  
3. **Asiento del lastre** con perno O (objeto).  
4. Caras superior/inferior alineables con anclajes R del recorrido.

### 2.3 Lastre transferible L

```text
        ┌──────────┐
        │  cuerpo  │  acero / latón, masa m_L
        │  orificio pasante o orejetas
        │  ● asiento perno O (módulo)
        │  ● asiento perno R (guía alta o baja)
        └──────────┘
```

Requisitos:

- No se oxida fatalmente en agua (inox 316, latón, o acero + recubrimiento).  
- Geometría que **no se atasque** en la horquilla.  
- Dos puntos de enclavamiento claros (O y R).

### 2.4 Recorrido — anclajes R

| Anclaje | Cota | Función |
|---------|------|---------|
| R_alta | superior (bajo la tapa o puente) | retiene L cuando no va al módulo |
| R_baja | inferior (cerca del final de carrera) | aparca L tras la bajada |

En v1 ambos pueden ser **pasadores manuales** (palanca fuera del agua).  
En v1.1: solenoides estancos o actuadores lineales con fuelle.

### 2.5 Make-before-break (secuencia obligatoria)

**Enganche en ALTA (L pasa de guía → módulo):**

```text
1. Módulo en posición ALTA, centrado.
2. Cerrar O (módulo sujeta L)  →  O=1, R=1  (ambos)
3. Abrir R_alta                 →  O=1, R=0
4. Verificar: L viaja con módulo
```

**Entrega en BAJA (L pasa de módulo → guía):**

```text
1. Módulo en posición BAJA, centrado.
2. Cerrar R_baja               →  O=1, R=1
3. Abrir O                     →  O=0, R=1
4. Verificar: L queda en guía baja; módulo con n=3
```

**Nunca** O=0 y R=0 a la vez.

---

## 3. Calibración de flotabilidad (crítico)

### 3.1 Objetivo

Sea \(\rho\) densidad del agua ≈ 1000 kg/m³.

\[
F_b = \rho \cdot V \cdot g
\]

- Con masa \(m_3\): \( m_3 g < F_b \) (flota; margen 5–15 % de empuje sobrante).  
- Con masa \(m_4 = m_3 + m_L\): \( m_4 g > F_b \) (se hunde; margen 5–15 %).

### 3.2 Procedimiento de calibración (sin electrónica)

1. Construir casco cerrado (sin lastre transferible).  
2. Añadir lastre **fijo interno** hasta que flote con ~1–2 cm de francobordo (o apenas sumergido según diseño).  
3. Llamar a esa masa total \(m_3\) (equivale al estado “3 pesos”).  
4. Probar con lastre L enganchado: debe hundirse con decisión (no a flote inestable).  
5. Si no se hunde: subir \(m_L\) o bajar volumen (menos aire).  
6. Si con \(m_3\) ya se hunde: quitar lastre fijo o aumentar volumen.

### 3.3 Tabla de registro (rellenar en lab)

| Ensayo | m_3 (kg) | m_L (kg) | V estimado (L) | ¿Flota n=3? | ¿Hunde n=4? | Notas |
|--------|----------|----------|----------------|-------------|-------------|-------|
| 1 | | | | | | |
| 2 | | | | | | |
| 3 | | | | | | |

---

## 4. Planos lógicos por vistas

### 4.1 Planta (sección horizontal)

```text
        pared
    ┌─────────────────┐
    │  ·G1       ·G2  │   G1, G2 = guías
    │                 │
    │      [MÓDULO]   │
    │                 │
    └─────────────────┘
         0,40 m
```

### 4.2 Alzado (cotas)

```text
    z = 1,10 m ──── tapa / puente, R_alta
    z = 1,00 m ──── nivel agua libre
    z = 0,95 m ──── posición módulo ALTA (centro lastre)
         │
         │  Δh = 0,80 m (centro lastre ALTA → BAJA)
         │
    z = 0,15 m ──── posición módulo BAJA, R_baja
    z = 0,00 m ──── fondo interior
```

(Ajustar números al tanque real; mantener **un solo Δh medido** con cinta métrica entre centros de L en ALTA y BAJA.)

### 4.3 Detalle pernos (esquema)

```text
  GUÍA / ANCLAJE R          LASTRE L           MÓDULO / PERNO O
  ───────────────          ─────────          ────────────────
       [orificio]  ←pinR→  (orejeta)  ←pinO→  [solenoide/pasador]
```

**v1 manual:** pasadores de acero Ø4–6 mm + muelle + cable/bowden hasta fuera del tanque.  
**v1.1:** mini-solenoide push-pull 12 V, carcasa estanca IP68, o actuador lineal corto.

---

## 5. Lista de materiales (BOM) — v1 taller

Precios orientativos EUR (Europa, 2026, orden de magnitud). Ajustar a proveedores locales.

### 5.1 Estructura tanque y guía

| # | Ítem | Spec sugerida | Cant. | Est. € |
|---|------|---------------|-------|--------|
| T1 | Tanque / acuario | 400×400×1200 mm, cristal o PC | 1 | 40–120 |
| T2 | Varilla guía | Inox 316 Ø8–10 mm, L=1,0 m | 2 | 15–30 |
| T3 | Soportes guía superior/inferior | PLA+ o aluminio mecanizado | 4 | 10–25 |
| T4 | Puente superior | Madera contrachapada / Alu 2 mm | 1 | 10–20 |
| T5 | Tope final de carrera inferior | Goma + tornillos inox | 2 | 5 |
| T6 | Junta / sellado pasamuros (si hay cables) | Prensaestopas PG7/PG9 | 2–4 | 8–15 |

### 5.2 Módulo (objeto)

| # | Ítem | Spec sugerida | Cant. | Est. € |
|---|------|---------------|-------|--------|
| M1 | Casco / boya | Botella PET reforzada, o impresión 3D PETG hueca sellada, o PVC | 1 | 5–40 |
| M2 | Cojinetes / buje deslizante | PTFE / nylon para Ø guía | 4 | 8–20 |
| M3 | Horquilla asiento lastre | Impresión 3D PETG / nylon | 1 | 5–15 |
| M4 | Lastre fijo calibración | Plomos de pesca / arandelas inox | lote | 5–15 |
| M5 | Tornillería inox A4 | M3–M5 surtido | 1 | 8–12 |
| M6 | Silicona neutra / epoxy marina | sellado casco | 1 | 8–15 |

### 5.3 Lastre transferible y pernos

| # | Ítem | Spec sugerida | Cant. | Est. € |
|---|------|---------------|-------|--------|
| L1 | Cuerpo lastre | Disco/cilindro inox o latón ~0,5–1 kg | 1 | 10–30 |
| L2 | Pasador perno O | Ø5 mm inox + muelle | 1 | 5–10 |
| L3 | Pasador perno R_alta | idem | 1 | 5–10 |
| L4 | Pasador perno R_baja | idem | 1 | 5–10 |
| L5 | Bowden / cable acero | 1,5 mm, funda | 3 m | 8–15 |
| L6 | Palancas de mando | impreso 3D o aluminio | 3 | 5–15 |

### 5.4 Instrumentación mínima (recomendada)

| # | Ítem | Spec | Cant. | Est. € |
|---|------|------|-------|--------|
| I1 | Báscula 0–5 kg | 1 g resolución | 1 | 10–25 |
| I2 | Cinta métrica / pie de rey | — | 1 | 10 |
| I3 | Cronómetro | móvil | 1 | 0 |
| I4 | Cámara lenta (móvil) | transferencias | 1 | 0 |
| I5 | Célula de carga 5–10 kg + HX711 (opc.) | fuerza en guía | 1 | 15–30 |
| I6 | Arduino Nano / ESP32 (opc. v1.1) | log tiempos | 1 | 10–20 |

### 5.5 Opcional generación / freno (v1.1)

| # | Ítem | Spec | Cant. | Est. € |
|---|------|------|-------|--------|
| G1 | Dinamo / motor DC como generador | 12–24 V | 1 | 15–40 |
| G2 | Polea + hilo inextensible | acoplado al módulo | 1 | 10 |
| G3 | Resistencia de carga + multímetro | medir V,I | 1 | 20–40 |
| G4 | Diodo + condensador (snubber simple) | protección | lote | 5 |

### 5.6 Presupuesto orientativo

| Nivel | Contenido | Total aprox. |
|-------|-----------|--------------|
| **v1 mínima** | Tanque + guías + módulo + lastre + pernos manuales + báscula | **120–250 €** |
| **v1.1** | + solenoides/log + generador + célula carga | **250–450 €** |
| **v2** | Tanque mayor, multi-lastre, enjambre 2 módulos | **500–1200 €** |

---

## 6. Secuencia de fabricación y montaje

### Fase A — Estructura (1–2 días)

1. Colocar tanque nivelado; marcar cota de agua.  
2. Fijar guías G1–G2 verticales y paralelas (tolerancia < 1 mm en 1 m).  
3. Montar puente superior con R_alta.  
4. Montar soporte R_baja a la cota de fin de carrera.

### Fase B — Módulo (1–3 días)

1. Fabricar casco estanco; prueba de inmersión 30 min sin lastre transferible.  
2. Montar bujes; el módulo debe deslizar sin agarrotarse (seco y mojado).  
3. Montar horquilla + perno O.  
4. Calibrar m_3 (§3).

### Fase C — Lastre y transferencias (1–2 días)

1. Fabricar L con orejetas O y R.  
2. Ajustar pasadores: juego < 0,5 mm, enclavamiento positivo.  
3. Ensayar make-before-break en aire 50 veces.  
4. Ensayar en agua 20 veces (criterios C3–C4).

### Fase D — Medición de ciclo (1 día)

1. Medir Δh real (centros de L).  
2. Cronometrar bajada con n=4 y subida con n=3.  
3. (Opc.) medir energía eléctrica si hay generador.  
4. Devolver L a ALTA a mano; anotar esfuerzo/tiempo (recarga de batería).

---

## 7. Protocolo de ensayos

### Ensayo E1 — Flotación estática

- 10 min n=3 en reposo ALTA: no debe tender a bajar.  
- Enganchar L: debe iniciar descenso en < 5 s (o documentar si necesita empujón de despegue por fricción de guía).

### Ensayo E2 — 20 ciclos de transferencia

Para cada ciclo registrar:

| Ciclo | t_enganche_s | t_bajada_s | t_entrega_s | t_subida_s | Fallo (0/1) | Notas |
|------:|-------------:|-----------:|------------:|-----------:|:-----------:|-------|
| 1 | | | | | | |
| … | | | | | | |
| 20 | | | | | | |

**Pass:** 0 fallos de suelta de lastre; 20/20 completados.

### Ensayo E3 — Balance de peaje (honesto)

\[
E_{\text{ref}} = m_L g \Delta h
\]

| Concepto | Cómo se estima |
|----------|----------------|
| Energía disponible lastre | \(E_{\text{ref}}\) |
| Recuperable teórica | \(\eta_{\text{gen}}(1-f_{\text{drag}})E_{\text{ref}}\) |
| Coste recarga manual | tiempo × potencia media operador, o altura×peso si usas polea+dinamómetro |
| Coste pernos | si eléctricos: ∫VI dt |

Conclusión esperada: **sin lift automático, demuestras descarga de inventario**; **con recarga del lastre, el ciclo cerrado no es sobreunitario**.

---

## 8. Seguridad

| Riesgo | Mitigación |
|--------|------------|
| Agua + 12/24 V | Fuente limitada corriente; prensaestopas; GFCI/diferencial en red |
| Aplastamiento de dedos en pernos | Palancas fuera del tanque; no meter manos en transferencia |
| Rotura de tanque | No sobrellenar; apoyo rígido; no subirse al borde |
| Lastre suelto (0/0) | Interlock mecánico o checklist; nunca operar sin R u O |
| Óxido / contaminación | Inox/latón; secar tras uso; no beber el agua del tanque |

---

## 9. Qué fabricar en v1 vs dejar para después

| Incluir en v1 | Dejar para v2+ |
|---------------|----------------|
| 1 Δh, 1 lastre | Multi-lastre / stock N |
| Pernos manuales | Solenoides + PLC |
| Calibración 3/4 | Enjambre de módulos |
| Medición manual | Generación a red / BESS |
| Tanque de mesa | Instalación submarina real |
| Lift manual del lastre | Cinta automática de recarga |

---

## 10. Planos “entregables” de esta carpeta

Este documento **es** el plano lógico + BOM del MVP.

Entregables opcionales siguientes (si se piden):

1. `BOM.csv` — misma lista en hoja de cálculo  
2. Croquis cotas en SVG/PDF  
3. STL de horquilla + palancas (impresión 3D)  
4. Script de log serial para tiempos de ciclo  

### v1.1 — Solenoides + ESP32 (disponible)

Carpeta: `v1_1_ESP32/`

| Archivo | Contenido |
|---------|-----------|
| `README_v1_1.md` | Arquitectura, FSM, BOM electrónica, puesta en marcha |
| `kilometro_v11_esp32.ino` | Firmware Arduino/ESP32 (make-before-break + CSV log) |
| `wiring_v1_1.txt` | Cableado GPIO / MOSFET / FC |
| `BOM_v1_1_extra.csv` | Material electrónico adicional (~80–170 €) |

La recarga del lastre BAJA→ALTA sigue siendo **manual** en v1.1; el firmware automatiza enganche, espera de cotas y entrega.

---

## 11. Resumen ejecutivo del prototipo

```text
TANQUE + 2 GUÍAS
    + MÓDULO (casco calibrado 3 flota)
    + 1 LASTRE (el 4.º)
    + 3 ENCLAVAMIENTOS (O, R_alta, R_baja) con make-before-break
    + ESTACIÓN S siempre arriba (recarga geometría / enganche)
    + RECARGA DE BATERÍA = devolver lastre a ALTA (manual en v1)
```

**Éxito del prototipo:** demostrar 20 ciclos limpios de conmutación y el umbral 3/4.  
**No es fracaso:** que al reponer el lastre el balance energético no sea positivo: **eso confirma que es una batería**, no un motor.

---

## 12. Checklist de compra rápida (v1 mínima)

- [ ] Tanque ~40×40×100+ cm  
- [ ] 2 varillas inox Ø8–10 × 1 m  
- [ ] Material casco estanco + sellador  
- [ ] 4 bujes deslizantes  
- [ ] 1 masa 0,5–1 kg (lastre L)  
- [ ] Lastre fino para calibrar m_3  
- [ ] 3 pasadores + muelles + cable  
- [ ] Báscula y cinta métrica  
- [ ] Móvil para vídeo de transferencias  

**Presupuesto objetivo:** ~150–250 € si ya hay herramientas básicas e impresora 3D.

---

*Prototipo tanque Kilómetro — planos lógicos y BOM v1 — worktree `doble-y-perneado/kilometro_sim`.*
