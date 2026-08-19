# Avances matemáticos  
## MDC, geometría, ofuscación por decimales de π y e — Documento técnico

**Versión:** 1.1 (alineado con pack 33×1 y `LIMITES_HONESTOS`)  
**Dominio:** teoría de números aplicada, factorización exploratoria, cifrado geométrico determinista, hashing/stream K3  
**Nivel de madurez:** implementaciones ejecutables (Python/C), demos AntiPC, motores de alta precisión `Decimal`.  
**Uso declarado:** civil (ver `02_USO_CIVIL.txt`). **No** ruptura de RSA industrial. **No** P=NP.

---

## 1. Resumen ejecutivo técnico

Este bloque reúne cuatro líneas matemáticas interconectadas del corpus VMA:

| Línea | Nombre | Función |
|-------|--------|---------|
| A | **MDC** — Método de Descomposición por Convergencia | Exploración de factores / fase / “diente de sierra” sobre semiprimos y candidatos |
| B | **Geometría** — restos, hipotenusa, Thales, polígonos | Visualización y motor de cifrado por perímetro acumulado |
| C | **Ofuscación π / e** | Aproximaciones de alta precisión usadas como capas de mezcla en el cifrado de fase |
| D | **K3 / Aleatorovix** | Motor de stream, hash y organismo de exploración con entropía y memoria MDC |

Todas comparten un principio de diseño:

> **Convergencia controlada + aritmética exacta o de precisión arbitraria**  
> en lugar de depender solo de heurísticas opacas no auditables.

---

## 2. MDC — Método de Descomposición por Convergencia

### 2.1 Idea central

Dado un entero \(n\) (en demos: semiprimo pequeño \(n=p\cdot q\)), se busca una representación en **dos trenes**:

\[
\boxed{(2x+3)(2y+3)=n}
\]

Si existen enteros \(x,y \ge 0\) que satisfacen la igualdad, se recuperan factores alineados con la rejilla impar desplazada:

\[
p' = 2x+3,\qquad q' = 2y+3
\]

(ajustando signos/paridad según la familia de candidatos; en la práctica se combina con la criba **6k±1**).

### 2.2 Capas operativas (AntiPC / plugins)

| Plugin / capa | Rol |
|---------------|-----|
| `MDC_FACTOR` | Búsqueda de pares \((x,y)\) o factores compatibles (**toy**, \(n\) pequeños) |
| `MDC_REGLA` | Escala logarítmica / planificador tipo “Battery Mode” |
| `MDC_PHASE` | Extrae una **fase** del payload para encadenar con K3 → hash |
| Memoria MDC v6 | Heat-map de sectores \(m\); evita re-visitar zonas “calientes” |

### 2.3 Perfil “diente de sierra” \(d(m)\)

Definición operativa (implementación Aleatorovix / `mdc_memoria.py`):

\[
d(m) = \mathrm{frac}\!\left(\frac{n}{2(2m+3)}\right)
= \frac{n \bmod \big(2(2m+3)\big)}{2(2m+3)}
\]

con aritmética exacta (`fractions.Fraction`) para evitar error de flotante.

Interpretación:

- \(d(m)\) cercano a 0 o a valores “especiales” → candidatos a factor / alineación modular.  
- Desfase \(|d-0{,}5|\) y aceleración \(\Delta d\) alimentan el **mapa de calor** de sectores.

### 2.4 Resonancia con inercia (puente MDC ↔ ZypyZape)

En el organismo Aleatorovix, la acción de tipo “inercia” ajusta el salto de exploración:

\[
k_{inercia}\ \text{crece con el “calor” del sector},\quad
\text{salto} \propto k_{inercia}
\]

Es una **metáfora operativa controlada**: igual que ZypyZape acelera/frena turbinas, MDC acelera/frena el índice \(m\) en la búsqueda. No implica física de red eléctrica dentro del factorizador.

### 2.5 Pipeline industrial (esquema)

```
Fichero  →  K3_FILE  →  K3_HASH  →  MDC_PHASE
Número   →  MDC_FACTOR  +  MDC_REGLA
```

Comandos de referencia:

```bat
cd C:\Users\cuent\Desktop\antipc
antipc.cmd mdc analyze --n 1147
antipc.cmd mdc visual
```

Demo canónica: \(1147 = 31 \times 37\).

### 2.6 Honestidad de alcance

| Afirmación correcta | Afirmación incorrecta (evitar) |
|---------------------|--------------------------------|
| MDC es una **herramienta de exploración** y arquitectura por capas | “Rompe RSA-2048” |
| Útil en demos, docencia, pipelines de integridad | Sustituto de GNFS / ECM en criptografía real |
| Integra fase + hash + memoria | Prueba de primalidad general completa por sí solo |

---

## 3. Geometría: restos, visualización y motor de perímetro

### 3.1 Geometría de primos (lado teoría de números)

Del corpus *Cribas, cotas y estructuras modulares*:

1. **Candidatos 6k±1** para primos &gt; 3.  
2. Intervalos \(I(n)\), longitudes \(L(n)\), marcados \(m(n)\approx(e-2)\sqrt{n}\).  
3. **Criterio operativo:** si \(L(n)-m(n)\ge 2\), se estima espacio para al menos dos primos en \(I(n)\).  
4. **Sistema de restos** \((x,a,b)\) y método gráfico con **hipotenusa** (visual pitagórico del MDC).  
5. **MRAUV:** modelo de densidad de primos por tramos (cinemática de \(\pi(x)/x\)).  
6. **Estructura Sofí** (Sophie Germain modular) y resultados de infinitud en clases \(L_i\).

Estas piezas alimentan cribas (desmemoriada, híbrida, modular) y estimadores de densidad (Criva).

### 3.2 Geometría de cifrado: motor binario → perímetro

Estado inicial compartido (**clave**):

| Componente | Ejemplo |
|------------|---------|
| Secuencia de Thales (escalas) | \([3,5,8,13,21]\) |
| Tipos de figura | equilátero / isósceles / escaleno |
| Puntos base por figura | \([6,12,18]\) |

#### Cálculo de lados (determinista)

\[
\mathrm{base} = \mathrm{puntos}\times 1{,}5
\]

| Figura | Lados (× escala Thales) |
|--------|-------------------------|
| Equilátero | \([\mathrm{base}/3,\ \mathrm{base}/3,\ \mathrm{base}/3]\) |
| Isósceles | \([\mathrm{base}/4,\ \mathrm{base}/4,\ \mathrm{base}/2]\) |
| Escaleno | \([0{,}25,\ 0{,}35,\ 0{,}40]\times\mathrm{base}\) |

#### Criterio de aportación por pares de bits

Se lee la cadena binaria en **pares** (bit huérfano final → se añade `0`):

| Par | Efecto en perímetro | Avance del lector |
|-----|---------------------|-------------------|
| `10` | \(+\) lado[0] | +1 bit |
| `11` | \(+\) lado[1] | +1 bit |
| `00` | \(+\) lado[0] + lado[1] | +1 bit |
| `01` | **sin aporte** (paradoja) | **+2 bits** |

El bit `01` introduce un **vacío de aportación** que rompe análisis de frecuencia ingenuos: el lector avanza, el perímetro no crece.

#### Salida

\[
P_{\mathrm{colapso}} = \sum_{\text{pasos}} \Delta L_i \in \mathbb{Q}\ \text{o}\ \texttt{Decimal}
\]

**Perímetro de colapso:** número real con muchos decimales que resume el camino geométrico.

### 3.3 Reversibilidad

| Condición | ¿Reversible? |
|-----------|--------------|
| Mismo estado inicial (Thales, figuras, puntos) + `Decimal` alta precisión | **Sí** (backtracking / resta de aportaciones) |
| Uso de `float` IEEE-754 | **No** (redondeo destruye el camino) |
| Clave incompleta | **No** práctica |

Requisito de implementación:

```python
from decimal import Decimal, getcontext
getcontext().prec = 100  # o 150–1000 en modo industrial
```

---

## 4. Ofuscación por decimales de π y e

### 4.1 Motivación

π y e son:

- **trascendentes** y de expansión decimal (o en base \(b\)) no periódica;
- generables por **algoritmos deterministas** (Taylor, producto de polígonos, series);
- útiles como **capas de mezcla** cuando se calculan con la **misma precisión y la misma semilla de parámetros** en emisor y receptor.

No se trata de “esconder un secreto en los dígitos públicos de π” (que cualquiera puede calcular), sino de:

1. construir una **aproximación parametrizada** \(\tilde\pi(p_1,p_2,N_{iter})\) y \(\tilde e(N)\);  
2. combinarlas con el **perímetro geométrico** y una semilla de fase;  
3. generar un **chorro de bits** (Aleatorovix) para XOR con el mensaje.

### 4.2 Aproximación de π por doblado de lados (implementación de referencia)

Esquema (archivo `gemini-code-1784158392232.py`):

1. Tomar \(n_{lados} = p_1\cdot p_2\).  
2. Ángulo \(= \tilde\pi_0 / n_{lados}\).  
3. Cuerda \(= 2\sin(\mathrm{ángulo})\) con **Taylor en `Decimal`**.  
4. Perímetro \(= n_{lados}\cdot\mathrm{cuerda}\); \(\tilde\pi = \mathrm{perímetro}/2\).  
5. Iterar doblando lados:

\[
n \leftarrow 2n,\qquad
\ell \leftarrow \frac{\ell}{\sqrt{2+\sqrt{4-\ell^2}}}
\quad\text{(forma equivalente en código)}
\]

Parámetro típico: \(N_{iter}=15\), precisión 150 decimales.

### 4.3 Aproximación de e convergente

Serie/producto iterativo en `Decimal` (50 términos en la demo):

\[
\tilde e = 1 + \sum_{v=1}^{N} t_v,\quad
t_{v} = t_{v-1}/c_v
\]

con coeficientes \(c_v\) definidos por la regla del script (rama par/impar).  
El punto no es competir con la serie factorial clásica en velocidad, sino **fijar un camino numérico reproducible** como capa de ofuscación.

### 4.4 Cadena de ofuscación (fase)

Flujo conceptual:

```
semilla_bits
    → perímetro geométrico (Thales + figuras)
    → mezcla con π̃(p1,p2,N)
    → mezcla con ẽ(N)
    → semilla de fase Decimal
    → Aleatorovix: chorro de bits (acordeón MSL/LSL)
    → XOR con payload
```

**Descifrado:** mismo camino forward con la clave, o **backtracking** sobre el residuo decimal cuando el diseño lo permite.

### 4.5 Límites de seguridad (declaración para terceros)

| Aspecto | Estado |
|---------|--------|
| Determinismo y auditabilidad | Alto (código abierto del motor) |
| Reversibilidad con clave | Sí (con `Decimal`) |
| Resistencia a adversario de alto nivel | **No afirmada** |
| Sustituto de AES/ChaCha/RSA | **No** |
| Uso adecuado | Integridad experimental, demos 33×1, educación, watermarking de fase |

---

## 5. Motor K3 (variante industrial de stream)

### 5.1 Acordeón en enteros (`k3_core.c`)

Estado `K3Motor { v, L }`:

```c
m->L += (m->v + 1);
m->v += 2;
if ((5 * m->L) <= (2 * m->v + 1))
    m->L += (m->v * 2);
stream = (m->L ^ m->v) * 0x9E3779B97F4A7C15ULL;
data[i] ^= (uint8_t)(stream & 0xFF);
```

- Constante dorada en forma entera `0x9E3779B97F4A7C15` (dispersión tipo splitmix).  
- Parámetros de proyecto **33×1:** `base=33`, `rel=1` en el launcher.  
- Operaciones ARX / XOR: adecuadas a implementación embebida.

### 5.2 Precauciones de ingeniería

1. El launcher de referencia **puede borrar el archivo original** tras cifrar: usar siempre copias de prueba.  
2. `mlock`/`memset` orientados a Linux; en Windows preferir WSL o reimplementación Python.  
3. Validar solo en entornos controlados; no exponer claves en repositorios.

### 5.3 Teorema de amplificador de fase (K=3 XOR) — resumen

Documento formal: `THEOREM_PhaseAmplifier_K3_XOR.md`.

- Estado de 4 bits \((v_1,f_1,v_2,f_2)\).  
- `shift` tipo Gray sobre valores; fases fijas.  
- Suma de tres caminos Toffoli XOR-eados.  
- **Invariante:** \(f_2\) se conserva → partición del espacio de estados.  
- Uso: diseño de realimentación de fase booleana, no cifrado de red por sí solo.

---

## 6. Aleatorovix como organismo de exploración

No es un PRNG clásico aislado: es un **organismo de software** que:

1. Toma entropía de bajo nivel (nanos de reloj, pila, heap; ping opcional).  
2. Aplica **máscara Lila** (curvatura + campanas + intérprete mutante).  
3. Decide acciones:
   - 0 — Criba desmemoriada (olvido en RAM),  
   - 1 — Salto 6k+1,  
   - 2 — Salto 6k−1,  
   - 3 — Resonancia inercial MDC / ZypyZape.  
4. Borra rastro operativo tras el ciclo (desmemoria).

Relación con MDC: la memoria de sectores calientes/fríos guía **dónde** no insistir en la búsqueda.

---

## 7. Integración en el ecosistema VMA

```
┌──────────────────────────────────────────────┐
│              Capa matemática                 │
│  Cribas 6k±1 · Criva · MRAUV · Sofí · MDC    │
└───────────────────┬──────────────────────────┘
                    │ fase / d(m) / primos
┌───────────────────▼──────────────────────────┐
│           Capa cifrado / hash                │
│  Geometría perímetro · π̃ · ẽ · K3 · XOR      │
└───────────────────┬──────────────────────────┘
                    │ plugins / integridad
┌───────────────────▼──────────────────────────┐
│              AntiPC / runtime                │
│  K3_HASH · MDC_* · telemetría · demos        │
└───────────────────┬──────────────────────────┘
                    │ metáfora de control
┌───────────────────▼──────────────────────────┐
│     ZypyZape (resonancia inercial en soft)   │
└──────────────────────────────────────────────┘
```

---

## 8. Requisitos de reproducción (laboratorio software)

| Componente | Requisito |
|------------|-----------|
| Python | ≥ 3.10, `decimal` stdlib |
| AntiPC | `antipc.cmd`, plugins MDC |
| Geometría | `encriptacionGeometrica\` |
| Precisión | `getcontext().prec ≥ 100` |
| C (opcional) | `gcc` + WSL/Linux para `k3_core` |
| Datos de prueba | \(n=1147\); bits de demo `0101101011000010` |

### Checklist de auditoría

1. Ejecutar `mdc analyze --n 1147` → factores 31×37.  
2. Cifrar/descifrar un payload corto con motor geométrico + `Decimal`.  
3. Verificar que `float` rompe la reversibilidad (prueba negativa).  
4. Generar \(\tilde\pi\) dos veces con mismos \((p_1,p_2,N)\) → igualdad bit a bit en `Decimal`.  
5. No ejecutar launcher K3 sobre archivos sin backup.

---

## 9. Métricas de calidad matemática (propuestas)

| Métrica | Objetivo |
|---------|----------|
| Reproducibilidad | 100 % con misma semilla y `prec` |
| Error de cierre MDC toy | 0 en demos (\(n\) factoriza) |
| Falsos positivos de fase | Documentar tasa en corpus de prueba |
| Tiempo MDC vs \(n\) | Curva empírica (no polinomial garantizada) |
| Distancia de chorro Aleatorovix | Tests de frecuencia / rachas (batería ligera) |

---

## 10. Límites honestos (mates)

| Afirmación | Estado |
|------------|--------|
| Estructuras Sofí / infinitudes en archivo de teoremas | **[OK en corpus]** — revisión externa pendiente |
| Goldbach completo | **[CONJETURA]** |
| MRAUV ⇒ Goldbach teorema | **[HEUR]** |
| MDC / Método-V = primalidad general | **REFUTADO / toy** |
| P vs NP resuelto por el pack | **NO archivado** |
| AntiPC benchmarks | **[OK]** si entorno reproducible |
| K3 / geometría como AES de producción | **No afirmado** |

**Conexión con energía:** la “resonancia inercial” MDC↔ZypyZape es **metáfora de control de búsqueda**, no física de red.  
El pack mates **no desnucleariza**; aporta integridad, educación y software civil del “1” auditable.

---

## 11. Conclusiones técnicas

1. **MDC** aporta un marco de **exploración modular y por fase**, integrable en pipelines; es **toy** en factorización.  
2. La **geometría de perímetro + Thales** produce cifrado determinista **reversible con clave y `Decimal`**.  
3. Las capas **π y e** son **aproximaciones parametrizadas** de alta precisión, no magia de dígitos públicos.  
4. **K3** es stream compacto alineado con 33×1, con alcance criptográfico **no** nation-state.  
5. Valor para terceros: **auditabilidad y capas** — no titulares de problemas abiertos “resueltos”.

---

## 12. Referencias internas de código y texto

| Artefacto | Ubicación típica |
|-----------|------------------|
| Skill / regla MDC | plugins `mdc_*` en `antipc\` / `repo\antipc\` |
| Memoria MDC v6 | `just run\aleatorovix\nucleo\mdc_memoria.py` |
| Cribas y cotas | `just run\archivos-vma\cribas_cotas_vma.txt` |
| Motor geométrico (teoría) | `encriptacionGeometrica\encriptacion.txt` |
| π / e / Aleatorovix | `encriptacionGeometrica\gemini-code-1784158392232.py` |
| K3 C + launcher | `encriptacionGeometrica\k3_core.c`, `k3_launcher.py` |
| Teorema fase K=3 | `just run\teoremas\THEOREM_PhaseAmplifier_K3_XOR.md` |

---

*Documento 04/04 — Dossier técnico VMA. Solo contenido técnico.*
