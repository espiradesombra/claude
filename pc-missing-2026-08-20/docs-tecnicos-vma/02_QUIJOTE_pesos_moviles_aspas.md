# Quijote  
## Pesos móviles en aspas de aerogeneradores — Documento técnico

**Versión:** 1.1 (límites honestos + enlace a `kilometro_sim` / peaje de actuador)  
**Dominio:** modulación mecánica de inercia del rotor; dinámica local; acoplamiento opcional a ZypyZape  
**Nivel de madurez:** formalismo + gemelos; hardware de pala **conceptual**. Hurto gravitatorio **neto en ciclo cerrado = 0** (actuador + centrífuga dominan).

---

## 1. Resumen ejecutivo técnico

**Quijote** es un mecanismo de **masa desplazable radialmente** a lo largo de cada pala (o canal interno de pala). Al variar la posición radial \(r_q(t)\) de una masa \(m_q\), se modula el momento de inercia del rotor:

\[
J_{total}(t) = J_G + \sum_{k=0}^{N-1} m_{q,k}\, r_{q,k}(t)^2
\]

donde \(J_G\) es la inercia “estructural” del rotor y \(N\) el número de palas.

### Funciones técnicas

| Función | Descripción |
|---------|-------------|
| Batería mecánica local | Almacenar/restituir energía cinética vía \(\Delta J\) a \(\omega\) casi constante o controlada |
| Respuesta rápida (FFR) | Modo “ball” (oscilación radial) → alta \(dJ/dt\) → par de reacción elevado en corto tiempo |
| Respuesta sostenida (PFR) | Modo estático \(r = r_{max}\) o \(r_{min}\) → “bateria mecánica” de menor potencia, mayor duración |
| Hurto gravitatorio (hipótesis) | Aprovechar componente \(mg\) a lo largo de la pala al bajar/subir la masa en fases seleccionadas del ángulo \(\theta\) |

### Relación con ZypyZape

Quijote **no sustituye** el control de par de convertidor. Es una **capa mecánica opcional** que:

- actúa en el dominio de \(J(t)\) y del par de reacción \(-\omega\,dJ/dt\);
- puede sincronizarse con el ciclo ACEL/FREN de ZypyZape;
- aporta más valor en **microredes y rotores individuales** que en redes multi-GW con \(m_q\) pequeñas.

---

## 2. Geometría y grados de libertad

### 2.1 Pala y canal

Modelo de referencia (gemelos NREL 5 MW / v4.8):

| Parámetro | Símbolo | Valor |
|-----------|---------|-------|
| Radio rotor | \(R\) | 60–63 m |
| Radio hub (inicio canal) | \(R_{hub}\) | 5 m |
| Radio tip útil (fin canal) | \(R_{tip}\) | 55–60 m |
| Masa por pala (referencia simple) | \(M_Q\) | 4 kg (modelo v4.8) |
| Fluido denso (modelo paper) | \(\rho_{fl}\) | ≈ 3386 kg/m³ (Fe+aceite) |
| Diámetro de canal | \(D\) | 0,05 m |
| Velocidad radial máx. | \(v_{slide,max}\) / \(DR_{MAX}\) | 0,5 m/s |
| Palas | \(N\) | 3 (estándar) o 7 (estudio comparado) |

Masa equivalente de fluido en un tramo de canal:

\[
m_q = \rho_{fl}\,\pi\left(\frac{D}{2}\right)^2 (r_2 - r_1)
\]

### 2.2 Posición sinusoïdal (modo ball, \(N\) palas)

\[
r_k(t) = r_0 + A\sin\!\left(\omega_{rot}t + \frac{2\pi k}{N}\right),\quad k=0,\ldots,N-1
\]

con \(r_0 = (r_{max}+r_{min})/2\) y \(A = (r_{max}-r_{min})/2\).

---

## 3. Propiedad fundamental: \(J_{total}\) constante en el ball equilibrado

### 3.1 Demostración (esquema)

\[
J_{total}(t) = m\sum_{k=0}^{N-1} r_k(t)^2
= m\sum_k\left(r_0 + A\sin\phi_k\right)^2
\]

\[
= m\left[ N r_0^2 + 2 r_0 A\sum\sin\phi_k + A^2\sum\sin^2\phi_k \right]
\]

Para fases equiespaciadas \(\phi_k = \theta + 2\pi k/N\), \(N\ge 2\):

\[
\sum_{k=0}^{N-1}\sin\phi_k = 0,\qquad
\sum_{k=0}^{N-1}\sin^2\phi_k = \frac{N}{2}
\]

Por tanto:

\[
\boxed{J_{total} = m\,N\left(r_0^2 + \frac{A^2}{2}\right) = \text{constante en } t}
\]

**Analogía electrotécnica:** la potencia trifásica instantánea es constante; aquí el “cruzado” se anula igual que en sistemas simétricos.

### 3.2 Factor geométrico tipo \(\sqrt{3}\)

\[
\text{Factor}_N = 2\sin\!\left(\frac{\pi}{N}\right)
\]

| \(N\) | \(\text{Factor}_N\) |
|-------|---------------------|
| 3 | \(2\sin 60^\circ = \sqrt{3} \approx 1{,}732\) |
| 6 | 1,000 |
| 7 | \(2\sin(\pi/7) \approx 0{,}868\) |
| 12 | ≈ 0,518 |

Aparece en relaciones de amplitud líder/seguidoras al imponer \(J\) constante bajo control no sinusoïdal.

### 3.3 Regla de control “líder / retractor”

Para mantener \(J\) constante al extender una pala líder en \(+\Delta A\):

- **1 pala** extiende \(A\);
- **\(N-1\) palas** retraen \(A/(N-1)\) cada una.

Verificación de primer orden en \(\Delta J \propto 2r\,\Delta r\):

\[
\Delta J_{ext} = m A (2r),\qquad
\Delta J_{ret} = (N-1)\, m \frac{A}{N-1}(2r) = m A (2r)
\quad\Rightarrow\quad \Delta J_{net}=0
\]

Ejemplos:

- \(N=3\): 1 extiende \(A\), 2 retraen \(A/2\).  
- \(N=7\): 1 extiende \(A\), 6 retraen \(A/6\).

---

## 4. Dinámica del rotor con \(J\) variable

Ecuación de movimiento con inercia variable:

\[
\frac{d}{dt}\big(J\omega\big) = T_{net}
\quad\Rightarrow\quad
J\frac{d\omega}{dt} = T_{net} - \omega\frac{dJ}{dt}
\]

El término \(-\omega\,dJ/dt\) es un **par de reacción**:

| \(dJ/dt\) | Efecto | Interpretación |
|-----------|--------|----------------|
| \(>0\) (masas hacia la punta) | Par negativo extra | “Cargar” la batería mecánica; frena el rotor |
| \(<0\) (masas hacia el hub) | Par positivo extra | “Descargar”; acelera o alivia al generador |

Para \(N\) palas con la misma \(v_{slide} = dr/dt\):

\[
\frac{dJ}{dt} = N\, m_q\, 2\, r\, v_{slide}
\]

### 4.1 Fuerza centrífuga

\[
F_{cent} = m_q\,\omega^2\, r
\]

Siempre hacia afuera. El actuador debe **vencer** (o aprovechar) esta fuerza al retraer o extender. Es el coste mecánico dominante del Quijote a \(\omega\) de aerogenerador.

Orden de magnitud (v4.8): \(F_{cent,medio} \approx 480\) N/pala con \(M_Q=4\) kg y \(r\sim 50\) m a \(\omega\sim 1{,}14\) rad/s.

### 4.2 Límite físico posicional (observación de diseño)

En un paso de tiempo, el desplazamiento radial no puede “saltar” de forma inconsistente con la cinemática del canal. En gemelo v9.4 se impone, entre otras, la cota:

\[
|dr_q| \le \omega\, r_q \quad\text{(límite posicional / cinemático documentado en código)}
\]

junto con \(DR_{MAX}\) de actuador (p. ej. 0,5 m/s).

---

## 5. Ball vs estático: dos modos de servicio

Sea \(\Delta r = r_{max}-r_{min}\) y \(\Delta(r^2) = r_{max}^2 - r_{min}^2\).

Energía cinética disponible por cambio de \(J\) (a \(\omega\) fijo):

\[
\Delta E = \frac{1}{2}\,\Delta J\,\omega^2 = \frac{1}{2}\, N\, m_q\, (r_{max}^2 - r_{min}^2)\,\omega^2
\]

**Importante:** \(\Delta E\) es **idéntico** en:

- **Caso A — Estático:** una maniobra a \(r_{max}\), espera, retorno a \(r_{min}\).  
- **Caso B — Ball:** oscilación continua entre \(r_{min}\) y \(r_{max}\).

La diferencia es **temporal**:

| Modo | Analogía | \(P\) instantánea | Duración útil | Desgaste actuador | Servicio de red |
|------|----------|-------------------|---------------|-------------------|-----------------|
| Ball | Condensador mecánico | Alta (\(dJ/dt\) grande) | Segundos (FFR &lt; 2 s) | Continuo | Fast Frequency Response |
| Estático | Batería mecánica | Baja | Decenas de s (PFR 10–30 s) | Mínimo | Primary Frequency Response |

\[
J_{medio,ball} = N m_q\left(r_{mid}^2 + \frac{A^2}{2}\right),\quad
r_{mid}=\frac{r_{max}+r_{min}}{2},\; A=\frac{\Delta r}{2}
\]

Potencia pico del ball (orden de magnitud):

\[
P_{pic} \approx N\, m_q\, 2\, r_{mid}\, v_{slide,max}\,\omega^2
\]

---

## 6. Comparativa 3 palas vs 7 palas

A igual \(m_q\) por pala e igual \([r_{min}, r_{max}]\):

| Magnitud | Escalado | Ratio 7p / 3p |
|----------|----------|---------------|
| Número de quijotes | \(N\) | 2,333× |
| \(\Delta J_{max}\) | \(\propto N\) | 2,333× |
| \(\Delta E_{max}\) | \(\propto N\) | 2,333× |
| \(P_{pic}\) ball | \(\propto N\) | 2,333× |
| Zonas muertas angulares (heurística paper) | — | 3p ≈ 60°; 7p &lt; 25° |
| Pesos activos simultáneos (ball 7p) | — | Típicamente 3–4 vs 1–2 en 3p |

### Ejemplo numérico (valores del formalismo 3vs7)

- \(m = 4\) kg/pala, \(\omega = 1{,}14\) rad/s, \(r_{min}=5\) m, \(r_{max}=60\) m  

\[
\Delta E_{3p} \approx 28\ \text{kJ},\qquad
\Delta E_{7p} \approx 65\ \text{kJ}
\]

Fracción respecto al ciclo ZZ típico (~406 kJ/turbina):

\[
\frac{\Delta E_{7p}}{E_{ZZ}} \approx 16\%\ \text{por turbina (sin escalar flota)}
\]

---

## 7. Control y “timing” del viento

### 7.1 Control simple (v4.8)

\[
v_{slide} = \mathrm{clip}\big(K_\omega\Delta\omega + K_f\Delta f,\; -V_{max},\; +V_{max}\big)
\]

- Hacia afuera → carga (\(J\uparrow\)).  
- Hacia adentro → descarga (\(J\downarrow\)).

Resultado típico simulado con viento fuerte: \(r_{medio}\) tiende a la punta (~49 m) por centrífuga.

### 7.2 Regla de timing (observación de operación)

Mientras se **carga** (masas hacia afuera):

- el término \(-\omega\,dJ/dt\) **resta** par al rotor;
- la potencia a red puede bajar si no se compensa el generador.

Por tanto:

| Condición de viento | Cargar \(J\) | Descargar \(J\) |
|---------------------|--------------|-----------------|
| Viento fuerte / sobrante de potencia | Preferible | Según red |
| Viento flojo | Evitar (roba lo poco que hay) | Preferible si se necesita par |
| Anticipación \(dv/dt &lt; 0\) | Retractar **antes** de que \(\omega\) colapse | Maximiza \(P_q \propto \omega^2 |dJ/dt|\) |

### 7.3 Política angular “hurto gravitatorio” (reglas de paper)

Una política de referencia por ángulo de pala \(\theta_k\):

\[
r_k =
\begin{cases}
r_{max} & \text{si }\cos\theta_k &gt; 0 \quad\text{(pala en fase “bajando”)}\\
r_{min} & \text{si }\cos\theta_k &lt; 0 \quad\text{(pala en fase “subiendo”)}
\end{cases}
\]

Trabajo gravitatorio idealizado por pala y vuelta (esquema paper):

\[
W_{hurto,k} = m_q\, g\, \Delta r \oint |\cos\theta|\,d\theta = m_q\, g\, \Delta r \cdot 4
\]

**Advertencia de integridad:** este \(W_{hurto}\) es un término de trabajo de campo gravitatorio sobre la masa. El **actuador** debe suministrar el trabajo contra fricción y, sobre todo, **contra la centrífuga**. Auditorías de gemelo (1ª ley) muestran que el **coste de actuador puede superar** el trabajo gravitatorio recuperable; el factor depende de \(m_q\), \(\omega\), \(\Delta r\) y de la velocidad de actuación. Por tanto:

> Quijote **no se presenta aquí como motor de energía neta gravitatoria**.  
> Se presenta como **modulador de inercia y de par de reacción**, con balance energético del actuador a cerrar en cada diseño.

---

## 8. Balance energético del actuador (cierre de 1ª ley)

Para un intervalo \(\Delta t\):

\[
\Delta E_{cin,rotor} + \Delta E_{pot,masas} + E_{fric} + E_{elec,gen}
= E_{viento} + E_{actuador,net} + \cdots
\]

En simulaciones con buffer hidráulico (gemelo v9.4):

- canal + fluido Fe+aceite;
- acumulador de presión (p. ej. 20 bar) y generador hidráulico \(\eta\approx 0{,}85\);
- válvulas de retención;
- pitch asíncrono por pala con términos \(K_p\cos + K_\omega\sin - K_Q\,dJ/dt\).

Criterio de viabilidad del “hurto”:

\[
P_{net} = P_{hurto} - P_{actuadores}
\]

Si \(P_{net} &lt; 0\) de forma sistemática → el valor del sistema es **control de \(J\)**, no generación neta.

Resultados de orden de magnitud documentados en contexto v4.8 (módulo 10 turbinas, 5 con Quijote, \(M_Q=4\) kg):

| Métrica | Valor |
|---------|-------|
| \(\Delta E_{quijote}\) | ≈ 96 kJ (≈ 23,7 % del ciclo ZZ por turbina en ese escenario) |
| \(\Delta H_{ZZ}\) en red 2 GW | ≈ 0,00014 s → **invisible a escala multi-GW** |
| Conclusión operativa | Valor **local** (rotor / microred 200–500 MW / islas) |

---

## 9. Diseño mecánico conceptual

### 9.1 Opciones de actuador

| Opción | Pros | Contras |
|--------|------|---------|
| Masa sólida en guía lineal | Simple, predecible | Desgaste, sellado, un solo \(m_q\) |
| Fluido denso en canal (Fe+aceite) | Masa distribuida, bombeable a buffer | Sedimentación, viscosidad, fugas |
| Buffer hidráulico + acumulador | Recupera parte de la energía de actuación | Complejidad, mantenimiento en pala |

### 9.2 Requisitos de pala

1. Canal estructural compatible con fatiga IEC 61400-1.  
2. Centro de masa y balanceo dinámico del rotor (no introducir desequilibrios no controlados).  
3. Seguridad: bloqueo en \(r_{hub}\) o \(r_{tip}\) ante fallo de control.  
4. Hielo, rayos, humedad y mantenimiento en altura.

### 9.3 Integración eléctrica / SCADA

- Señales: \(r_{q,k}\), \(P_{acc}\), \(Q_{bomba}\), alarmas de fuga.  
- Coordinación con pitch y con el ciclo ZypyZape (misma trama de tiempo o desfase optimizado).  
- Modo degradado: Quijote bloqueado, turbina en control convencional.

---

## 10. Plan de validación técnica

| Etapa | Objetivo | Métrica de éxito |
|-------|----------|------------------|
| Banco rotativo de 1 pala | Medir \(F_{cent}\), fricción, \(dJ/dt\) real | Error &lt; 10 % vs modelo |
| Rotor de laboratorio 3 palas | Verificar \(J_{total}\) constante en ball | \(\sigma_J / J &lt; \varepsilon\) |
| Auditoría 1ª ley | Cerrar balance actuador vs gravitatorio | \(P_{net}\) medido y trazable |
| Piloto en góndola / pala real | Fatiga, SCADA, eventos de red | Sin incidencias de seguridad; efecto local en \(\omega\) |

---

## 11. Tabla de límites honestos

| Afirmación | Estado |
|------------|--------|
| \(J_{total}\) constante en ball sinusoïdal equiespaciado | **[OK]** matemática |
| \(\Delta E\) ball = \(\Delta E\) estático (misma ventana \(r\)) | **[OK]** |
| Hurto gravitatorio con excedente neto ciclo cerrado | **[REFUTADO / W=0]** si se cierra el estado |
| Coste de actuador (centrípeta) vs \(mg\Delta r\) | A menudo **actuador gana** (más caro) |
| RoCoF 0–500 ms mejor que BESS | **[NO]** |
| Valor en microred / rotor individual | **Plausible** |
| Hardware en pala comercial | **No archivado** |

### Puente con `kilometro_sim`

La lección del **recorrido reduccionista** y del rearme en parado aplica a Quijote:  
**no pelear la centrífuga** cuando se puede rearmar \(r\) con \(\omega\approx0\).  
El cuello real (también en molinos “de la Tierra” vs diseños gigantes) es:

\[
\text{par de perturbación (viento/gravedad)} \quad\text{vs}\quad \text{par de actuador disponible}
\]

no un truco de energía neta.

---

## 12. Conclusiones técnicas

1. La matemática del ball equilibrado es **sólida**: \(J_{total}\) puede mantenerse constante con fases equiespaciadas (análogo trifásico).  
2. Ball y estático almacenan la **misma** \(\Delta E\); se eligen según el **perfil temporal** del servicio (FFR vs PFR).  
3. Escalar de 3 a 7 palas escala \(\Delta J\) y \(\Delta E\) como \(N\), con menos zonas muertas.  
4. En redes grandes, con masas pequeñas, Quijote es **marginal para \(H_{sys}\)**; su caso de uso es **dinámica de rotor y microredes**.  
5. El “hurto gravitatorio” debe presentarse con **balance de actuador cerrado**; no como energía neta garantizada.  
6. En desnuclearización (doc 05): Quijote es **I+D de inercia variable**, no sustituto de baseload.

---

## 13. Referencias internas de código

| Artefacto | Ubicación típica |
|-----------|------------------|
| Matemática 3 vs 7 | `just run\zypyzape-contexto\02_MATH_QUIJOTE_3vs7.txt` |
| Reglas paper + simulación | `just run\hurto-gravitatorio\gemell_quijote_paper_rules.py` |
| Gemelo completo v9.4 | `just run\gemelos\gemelo_v94.py` |
| Twin v4.8 integrado ZZ+Quijote | `just run\gemelos\zypyzape_twin_v4_8_quijote.py` |

---

*Documento 02/04 — Dossier técnico VMA. Solo contenido técnico.*
