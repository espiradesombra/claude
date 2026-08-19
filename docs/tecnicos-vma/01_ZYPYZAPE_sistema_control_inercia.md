# ZypyZape  
## Sistema de control y máquinas de inercia — Documento técnico

**Versión:** 1.1 (dossier + límites honestos 33×1 / `teoremasgrok`)  
**Dominio:** parques eólicos existentes, servicios ancilares de frecuencia, inercia sintética  
**Nivel de madurez (TRL ~4–5):** simulación / gemelo digital (calibración NREL 5 MW). **Sin piloto de campo archivado.**

---

## 1. Resumen ejecutivo técnico

**ZypyZape** es un **algoritmo de control de par** para aerogeneradores ya instalados. Convierte grupos de turbinas en una **batería cinética distribuida**: mientras una unidad reduce su velocidad angular y entrega potencia a la red (rol **FREN**), otra reduce el par de generador y acumula energía cinética (rol **ACEL**). El intercambio es cíclico, coordinado y reversible.

### Propuesta de valor técnica

| Aspecto | Descripción |
|---------|-------------|
| Hardware base | Ninguno en la versión software-only: actúa sobre convertidores back-to-back existentes vía firmware / capa de control |
| Activo reutilizado | Momento de inercia del rotor ya instalado \(J\omega^2/2\) |
| Servicio principal | Aumento de **inercia efectiva** de red y mejora de **nadir** / **RoCoF** ante desequilibrios de potencia |
| Topología | Módulos de \(N\) turbinas (referencia: \(N=10\)), emparejamiento y acoplamiento tipo **Kuramoto** |
| Alineación regulatoria (objetivo) | ENTSO-E NC RfG — inercia sintética / grid support |

### Lo que no es

- No es un Megapack electroquímico 1:1 (no almacena MWh de litio).
- No sustituye condensadores síncronos en todos los escenarios de grid code.
- No garantiza, por sí solo, un aumento neto de captación eólica (AEP) en campo: esa hipótesis permanece abierta.

---

## 2. Principios físicos

### 2.1 Energía cinética rotacional

Cada turbina \(i\) almacena:

\[
E_{c,i} = \frac{1}{2} J_i\,\omega_i^2
\]

El control de par del generador \(T_{gen}\) actúa sobre la ecuación de movimiento del rotor:

\[
J\,\frac{d\omega}{dt} = T_{aero} - T_{gen} - T_{roz} + T_{ZZ} + T_{IS}
\]

donde \(T_{ZZ}\) es la componente de par del ciclo ZypyZape e \(T_{IS}\) eventuales contribuciones de inercia sintética emulada en el convertidor.

### 2.2 Leyes de Faraday y Lenz (máquina eléctrica)

- **Faraday:** la tensión inducida escala con la velocidad \(\omega\); la potencia eléctrica entregada depende de \(T_{gen}\cdot\omega\).
- **Lenz:** al aumentar el par de frenado (inyección a red), el rotor se resiste y la “opacidad” aerodinámica efectiva se modula vía \(\lambda = \omega R / v_{viento}\).

### 2.3 Modelo aerodinámico (región MPPT)

\[
P_{aero} = \frac{1}{2}\,\rho\,A\,v^3\,C_p(\lambda,\beta)
\]

Parámetros de referencia del gemelo (familia NREL / modelo 2,5–5 MW):

| Parámetro | Símbolo | Valor de referencia |
|-----------|---------|---------------------|
| Radio rotor | \(R\) | 60–63 m |
| Inercia rotor | \(J_G\) | \(5\times10^6\) kg·m² (modelo 2,5 MW) / \(3{,}54\times10^7\) (NREL 5 MW en gemelo v9.4) |
| Potencia nominal | \(S_{nom}\) | 2,5 MW o 5 MW según gemelo |
| \(C_p\) máximo | \(C_{p,max}\) | ≈ 0,48–0,486 |
| TSR óptimo | \(\lambda_{opt}\) | 7,55 |
| Densidad aire | \(\rho\) | 1,225 kg/m³ |

### 2.4 Ecuación de balanceo de la red (swing equation)

Para una red agregada de potencia base \(S_{total}\):

\[
H_{eff} = H_{sys} + H_{ZZ}
\]

\[
H_{ZZ} = \frac{\frac{1}{2}\sum_i J_i\omega_i^2}{S_{total}}
\]

\[
\frac{df}{dt} = \frac{(\Delta P + P_{gov})\,f_0}{2\,H_{eff}\,S_{total}} - D\,(f - f_0)
\]

A mayor \(H_{eff}\), menor **RoCoF** y, en general, mejor **nadir** de frecuencia ante una pérdida de generación \(\Delta P < 0\).

---

## 3. Arquitectura del sistema de control

### 3.1 Roles por turbina

| Rol | Acción de control | Efecto en \(\omega\) | Efecto en red |
|-----|-------------------|----------------------|---------------|
| **CAPT** | MPPT / par óptimo | Seguimiento de viento | Potencia “base” |
| **ACEL** | Reduce \(T_{gen}\) | \(\omega\) sube → carga cinética | Menos potencia instantánea a red |
| **FREN** | Aumenta \(T_{gen}\) | \(\omega\) baja → descarga cinética | Más potencia instantánea a red |

El ciclo alterna ACEL/FREN de forma coordinada entre pares (o anillos) de turbinas.

### 3.2 Ciclo ZypyZape (parámetros de referencia v4.x)

| Parámetro | Valor | Notas |
|-----------|-------|-------|
| Fracción de par del ciclo | \(P_{ZZ,frac} = 0{,}13\) | 13 % de \(S_{nom}\) por turbina en el tramo activo |
| Periodo de ciclo | \(T_{ciclo} = 2{,}5\) s | \(f_{ciclo} = 0{,}4\) Hz |
| Topología \(N=10\) | 2 pares centrales + anillos de 4 | Modular, escalable por bloques de 25 MW (10×2,5 MW) |
| Modulacíon | Trapezoidal en \(\lambda\) | Reduce ruido en tip-speed ratio frente a onda cuadrada |

Energía característica del ciclo por turbina (orden de magnitud):

\[
E_{ZZ} \approx P_{ZZ,frac}\cdot S_{nom}\cdot\frac{T_{ciclo}}{2}
\approx 0{,}13\cdot 2{,}5\times10^6\cdot 1{,}25 \approx 406\ \text{kJ}
\]

### 3.3 Acoplamiento Kuramoto

Las fases de rotor (o fases de control de ciclo) se acoplan:

\[
\frac{d\theta_i}{dt} = \omega_i + \frac{K}{N}\sum_{j=1}^{N}\sin(\theta_j - \theta_i)
\]

Propósito: sincronizar el “baile” de pares ACEL/FREN, reducir rizado de potencia del módulo y estabilizar la dinámica colectiva.

### 3.4 Control de pitch (calibración NREL)

Control PI de ángulo de pala \(\beta\):

| Ganancia | Valor típico gemelo |
|----------|---------------------|
| \(K_p\) | 2,0–2,5 °/(rad/s) |
| \(K_i\) | 0,40 |
| Rate limit | 8 °/s |
| Rango | \(\beta \in [0^\circ, 30^\circ]\) |

Regiones operativas:

1. Arranque  
2. MPPT (\(\beta = 0^\circ\))  
3. Nominal / limitación (PI de pitch activo)

---

## 4. Máquinas de inercia en el marco ZypyZape

ZypyZape organiza **tres capas** de “máquina de inercia”:

### Capa A — Inercia rotacional nativa (siempre presente)

El rotor ya es un volante de inercia. ZypyZape **gestiona** cuándo se carga y descarga esa inercia respecto a la red.

### Capa B — Inercia sintética de control (software)

Mediante droop y términos de \(d\omega/dt\) o \(df/dt\) en el convertidor, se emula respuesta inercial adicional. Es la vía de despliegue más rápida (firmware).

### Capa C — Inercia variable mecánica (opcional: Quijote)

Masas móviles en las palas varían \(J(t)\). Documentado en el documento 02. En red de 2 GW con masas de pocos kg, \(\Delta H\) es marginal; el valor es local.

### Capa D — Almacenamiento gravitatorio / flotación (opcional: Kilòmetre)

Módulo mecánico externo o submarino como buffer de energía. Documentado en el documento 03. No es necesario para la versión base de ZypyZape.

---

## 5. Modelo de red y resultados de simulación (honestos)

### 5.1 Caso de referencia

| Magnitud | Valor |
|----------|-------|
| \(S_{total}\) | 2 GW |
| \(H_{sys}\) | 4,0 s |
| \(f_0\) | 50 Hz |
| Droop | 5 % |
| Amortiguamiento \(D\) | 0,05 |
| Perturbación | \(\Delta P = -100\) MW |
| Módulo ZZ | 10 turbinas × 2,5 MW = 25 MW |

### 5.2 Resultados gemelo v4.6 / v4.8 (orden de magnitud)

| Métrica | Referencia (sin ZZ) | Con ZZ (1 módulo) |
|---------|---------------------|-------------------|
| \(f_{nadir}\) | ≈ 49,556 Hz | ≈ 49,558 Hz |
| Mejora de nadir | — | **≈ +0,002 Hz por módulo** |
| Escalado | — | **Lineal** con número de módulos (en el modelo agregado) |
| \(\eta\) global del ciclo (sim.) | — | ≈ 98 % (idealizado) |

**Interpretación:** el efecto por módulo en una red de 2 GW es **pequeño pero medible en simulación**. El valor industrial aparece al:

1. Desplegar muchos módulos en parques grandes, o  
2. Aplicar el control en **redes más pequeñas / islas / microredes**, donde \(H_{sys}\) es bajo y \(\Delta H_{ZZ}/H_{sys}\) es significativo.

### 5.3 Hipótesis abiertas (no demostradas en campo)

1. **\(\Delta C_p\) neto > 0** a lo largo del ciclo (ganancia de captación eólica neta).  
   Plausible por asimetrías en \(P\propto v^3\) y pérdidas, **no validada experimentalmente**.
2. Fatiga de tren de potencia y palas bajo ciclos ACEL/FREN a \(0{,}4\) Hz: requiere análisis IEC y OEM.
3. Armónicos, fault ride-through y re-certificación del convertidor tras cambio de firmware.
4. Interacción con protecciones SCADA existentes (anti-islanding, limitación de rampas).

---

## 6. Integración industrial

### 6.1 Ruta de implementación por fases

| Fase | Contenido | Entregable |
|------|-----------|------------|
| F0 | Gemelo digital (actual) | Scripts Python + gráficos de validación |
| F1 | Banco de pruebas 2–5 kW | Medida de RoCoF/nadir emulado, fatiga basica |
| F2 | Piloto 1–2 turbinas adyacentes | Trazas SCADA reales, AEP, eventos de red |
| F3 | Parque 5–10 MW | Optimización de topología y mercados ancilares |
| F4 | Producto (OEM + TSO) | Firmware certificado, API SCADA, O&M |

### 6.2 Interfaces

- **Entradas:** \(f\), \(df/dt\), \(\omega_i\), \(v_{viento}\), \(P_{grid}\), estados de pitch.  
- **Salidas:** referencias de par \(T_{gen,ref}\) (y opcionalmente \(\beta_{ref}\) asíncrono por pala en variantes avanzadas).  
- **Bus:** OPC-UA / MQTT / protocolo nativo del OEM.  
- **Supervisión:** gemelo digital en paralelo (digital twin) para predicción y anti-windup.

### 6.3 Cumplimiento y riesgos de fabricante (OEM)

| Riesgo | Mitigación técnica |
|--------|--------------------|
| Fatiga por ciclos de par | Límites de \(P_{ZZ,frac}\), rampas trapezoidales, conteo de ciclos |
| Garantía AEP | Modo “solo evento de red” vs. ciclo continuo; telemetría A/B |
| Grid code | Ensayos IEC 61400-21, FRT, THD |
| Responsabilidad firmware | Co-desarrollo o contenedor de control validado por OEM |

---

## 7. Relación con Quijote y Kilòmetre

```
                    ┌─────────────────────┐
                    │   Red / TSO / DSO   │
                    └──────────▲──────────┘
                               │ P, f, RoCoF
                    ┌──────────┴──────────┐
                    │      ZypyZape       │
                    │  (control de par)   │
                    └───┬────────────┬────┘
                        │            │
              ┌─────────▼──┐   ┌─────▼──────────┐
              │  Quijote   │   │  Kilòmetre     │
              │ J(t) palas │   │ buffer grav/flot│
              └────────────┘   └────────────────┘
```

- **ZypyZape solo** → despliegue más barato y rápido (software).  
- **+ Quijote** → más control de \(J\) local; coste y certificación de pala.  
- **+ Kilòmetre** → buffer mecánico externo / submarino; otro dominio de prototipo.

---

## 8. Criterios de aceptación propuestos (piloto)

1. Reducción medible de RoCoF medio en eventos de desequilibrio ≥ X % respecto a baseline.  
2. Mejora de nadir ≥ Y mHz (a fijar con el operador según \(S_{total}\) del piloto).  
3. No degradación de AEP superior a umbral contractual (p. ej. ≤ 0,5–1 % anual) en modo continuo, o AEP neutro en modo “solo evento”.  
4. Cumplimiento de límites de armónicos y de fatiga de diseño del OEM.  
5. Gemelo digital con error acotado vs. SCADA (p. ej. \(P\), \(\omega\) dentro de bandas acordadas).

---

## 9. Tabla de límites honestos (cruzar antes de vender)

Fuente: `33x1\repo\teoremasgrok\segunda_vuelta\LIMITES_HONESTOS.txt` e `ingenieria\01_ZYPYZAPE.txt`.

| Afirmación | Estado |
|------------|--------|
| Gemelo NREL 5 MW reproducible; conservación en ciclo de control | **[OK]** |
| Mejora de sincronismo de parque (Kuramoto) | **Estructural / métrica diferencial clara en sim** |
| +1,4–1,5 % potencia de red “absoluta” | **Dudoso** — no usar como titular |
| Mejora nadir ~mHz por módulo en red multi-GW | **Simulación** — real pero pequeña por módulo |
| Gana RoCoF 0–500 ms frente a BESS | **[NO]** — la electrónica es más rápida |
| Sustituye baseload nuclear | **[NO]** — solo servicios de inercia/coordinación |
| Piloto de campo medido | **Pendiente** |

### Rol en desnuclearización / transición (técnico)

Si se reduce generación síncrona (carbón/nuclear), cae \(H_{sys}\). ZypyZape puede **parcialmente** compensar con inercia gestionada en eólica existente y, en emplazamientos apagados, inspirar **campus de máquinas rotativas** (ver documento 05).  
**No** aporta los MWh baseload de un reactor.

---

## 10. Conclusiones técnicas

1. ZypyZape es **físicamente coherente** como gestor de inercia y potencia de corta duración en rotores existentes.  
2. El gemelo digital muestra **mejora pequeña pero lineal** de nadir en red 2 GW; el caso de uso fuerte es **agregación de módulos** o **redes débiles**.  
3. La ruta industrial prioritaria es **software-only + piloto**, dejando Quijote/Kilómetro como fase 2 opcional.  
4. La credibilidad ante terceros depende de **no inflar AEP/RoCoF** y de un plan de fatiga/certificación con OEM.  
5. En transición energética: valor como **servicio auxiliar**, no como “apagador de nucleares”.

---

## 11. Referencias internas de código

| Artefacto | Ubicación típica |
|-----------|------------------|
| Contexto y parámetros | `just run\zypyzape-contexto\01_CONTEXT_ZYPYZAPE.txt` |
| Gemelo paper | `just run\gemelos\gemelo_v94.py` |
| Twin v4.8 + Quijote | `just run\gemelos\zypyzape_twin_v4_8_quijote.py` |
| Validación NREL / gráficos | dossier `zypyzape-contexto\dossier-elewit\` |

---

*Documento 01/04 — Dossier técnico VMA. Solo contenido técnico.*
