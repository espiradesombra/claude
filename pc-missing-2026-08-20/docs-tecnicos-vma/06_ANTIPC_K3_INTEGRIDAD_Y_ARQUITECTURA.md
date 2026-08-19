# AntiPC, K3 e integridad del pack “1”  
## Capa software industrial + cómo auditar lo ofrecido — Documento técnico

**Versión:** 1.0  
**Por qué es el documento 06:** en el catálogo de ingeniería VMA, **K3/AntiPC** tienen el **TRL más alto cercano a producto software** (≈5–6 / 4–5) y son lo que un tercero puede **ejecutar y hashear hoy**, sin pala ni tanque. Sin esta capa, los informes 01–05 no tienen **aval de integridad** ni puente cómputo↔energía.  
**Uso:** civil únicamente (`02_USO_CIVIL.txt`). K3 **no** se vende como sustituto de SHA-3/AES de grado gobierno.

---

## 1. Resumen ejecutivo

| Pieza | Qué es | TRL (honesto) | Estado |
|-------|--------|---------------|--------|
| **AntiPC** | Runtime / arquitectura de cómputo en red local: menos copias, reuse de conocimiento, hubs UDP | 4–5 | Benchmarks A→E reproducibles |
| **K3 / k3hash** | Huella industrial, dedup, streaming, marca de estructura | 5–6 | DLL/Python; **no** criptografía nation-state |
| **Manifiesto SHA-256** | Integridad del pack 33×1 / monorepo | Operativo | Regenerable en disco |
| **MDC plugins** | Fase / regla / factor toy en el mismo pipeline | 4 (demo) | Ver doc 04 |
| **XFI** | Gemelo experimental ZZ en craft N=3/4 | 2–3 | Conceptual / simulación |

**Propuesta de valor a un interesado técnico (software):**

1. Demostrar **métricas medibles** (copias memcpy, pkt/s, η de pipeline).  
2. Encadenar **fichero → K3 → MDC_PHASE → telemetría**.  
3. Verificar que el “1” que se ofrece **es el que está en disco** (hashes).  
4. Enlazar el mismo stack a **señales de control** tipo ZypyZape (UDP / Grafcet), sin pretender magia energética.

---

## 2. Por qué no bastaba con el doc 04

El documento 04 cubre **mates** (MDC, π/e, geometría).  
El 06 cubre **producto software y auditoría**:

```text
04  =  qué matemáticas hay y con qué límites
06  =  qué se ejecuta, cómo se mide, cómo se firma el pack, cómo se integra
```

Un revisor industrial suele preguntar primero:

- ¿Puedo correr algo en 15 minutos?  
- ¿Qué mejora mides?  
- ¿Cómo sé que no me cambias el ZIP mañana?

Eso responde el **06**, no el formalismo de Sofí.

---

## 3. Arquitectura unificada AntiPC (5 capas)

Fuente: `05-docs-tecnicos/ARQUITECTURA_UNIFICADA.txt`.

```text
L4  INDUSTRIAL     telemetría, perfiles, dedup, informes JSON/CSV
L3  RED            hubs UDP, WORK / RESULT / DONE, auth HMAC de puerta
L2  FLOW KERNEL    referencias, Knowledge Resolver, EventBus, Scheduler, plugins
L1  PLUGINS        K3_HASH / K3_FILE / K3_DEDUP · MDC_PHASE / REGLA / FACTOR · auth
L0  FUNDACIÓN      k3hash.c  ↔  hash_engine.py (DLL o puro Python)
```

### Plugins L1 (registro típico)

| PluginID | Función |
|----------|---------|
| `K3_HASH` | Hash streaming de payload |
| `K3_FILE` | Huella de fichero (tamaño + hash) |
| `K3_DEDUP` | Agrupación de duplicados |
| `MDC_PHASE` | Fase / curvatura modular desde payload |
| `MDC_REGLA` | Escala log / planificador |
| `MDC_FACTOR` | Factorización **toy** (números pequeños) |
| Network auth | Challenge HMAC antes de aceptar hub |

### Flujos

**Auditoría industrial**

```text
ruta → K3_FILE → KnowledgeResolver → (reuse | hash)
     → MDC_PHASE → telemetría
```

**Red distribuida**

```text
request_network_permission() → HMAC OK
start_network(hubs=N) → hubs + AUTH
dispatch_remote_hash(payload) → WORK → hub K3 → RESULT
```

**Honestidad de diseño:** la ALU no queda “mágicamente libre”; se reduce **coordinación y recopia** (memcpy de usuario, colas bloqueantes). El trabajo de hash **se hace** (local o en hub).

---

## 4. Benchmarks AntiPC (resultados de referencia)

### 4.1 Arquitecturas A→E (julio 2026, localhost, Python 3.12)

Config: 1000 paquetes, payload 64 B, 4 hubs, ~15 % duplicados.

| Arquitectura | Tiempo | Throughput | ALU rel. |
|--------------|--------|------------|----------|
| A — Convencional (memcpy + cola) | 0,605 s | 1654 pkt/s | 100 % |
| B — Lock-free ring | 0,579 s | 1726 pkt/s | 95,8 % |
| C — Distribuido 4 hubs (hilos) | 0,559 s | 1788 pkt/s | 92,5 % |
| D — Grafcet + existencia + redundancia | 0,908 s | 1102 pkt/s | 150 % (*) |
| E — UDP real (4 hubs) | **0,229 s** | **4361 pkt/s** | **37,9 %** |

(*) D es más “lento” en ese lote pero acumula **cache hits K(N)** (reuse de conocimiento) — otra métrica.

**Conclusión de referencia:** E descarga el maestro; el hash ocurre en hubs. Mejora de tiempo vs A ≈ **62 %** en ese entorno.  
**Nota:** valores absolutos dependen de CPU; lo válido es la **proporción en el mismo equipo**.

### 4.2 Viabilidad UDP ZypyZape (demo red, puerto 3333)

Captura `zypyzape_viability.json` (2 s, 4 hubs):

| Modo | Packets | `memcpy_user_copies` | Tiempo |
|------|---------|----------------------|--------|
| A convencional cola | 2845 | **5690** | ~2,00 s |
| B AntiPC slot ring | 2904 | **2904** | ~2,00 s |

Interpretación: mismo orden de paquetes; **~49 % menos copias de usuario** en B (5690 → 2904).  
Eso es el tipo de número que un integrador edge puede auditar.

### 4.3 Reproducir

```bat
cd C:\Users\cuent\Desktop\33x1\antipc
REM o rutas equivalentes del monorepo / sale-it
scripts\03_benchmark_udp.bat
python benchmark.py --packets 8000 --payload 256 --hubs 8
python udp_benchmark.py --packets 2000
```

Comandos pack 33×1 (`03_COMANDOS.txt`):

```bat
python cli.py hash --text "Hola"
python cli.py mdc analyze 143 --proper
python cli.py network demo --hubs 4 --duration 2
```

---

## 5. K3 — huella industrial (qué es y qué no es)

### 5.1 Problema que resuelve

- Fingerprint de buffers / ficheros / telemetría.  
- Deduplicación e indexación.  
- Marca estructural (proyecto 33×1: `base=33`, `rel=1` en variantes).  
- Enlace con phase-boolean / teorema K=3 XOR (capa formal, doc 04).

### 5.2 Lo que **no** es

| Afirmación incorrecta | Corrección |
|-----------------------|------------|
| “K3 es más seguro que AES-GCM” | **No afirmado**; no es el posicionamiento |
| “Rompe o sustituye SHA-3 en estándares” | **No** |
| “Listo para secreto de Estado” | **No** — uso civil, huella y dedup |
| “Cifrado de red por sí solo” | Es **hash / stream experimental**; ver geometría π/e en doc 04 para motores de fase |

### 5.3 Implementaciones en disco

| Ruta típica | Rol |
|-------------|-----|
| `antipc\hash_engine.py` / `k3hash.dll` | Runtime |
| `encriptacionGeometrica\k3_core.c` | Stream acordeón C |
| hashtool-extract / hashtool-work | CLI histórico HASHTOOLCODE |
| Teorema fase | `teoremas\19_THEOREM_PhaseAmplifier_K3_XOR.md` |

---

## 6. Integridad del pack “1” (manifiesto)

### 6.1 Qué es

El **1** del marco 33×1, en sentido técnico, es el **monorepo civil** (mates + energía + AntiPC + docs).  
El aval no es un eslogan: es **árbol de archivos + SHA-256**.

| Artefacto | Función |
|-----------|---------|
| `generar_manifiesto.py` | Recalcula hashes de piezas clave |
| `MANIFIESTO_HASHES.json` | Lista archivo → SHA-256 |
| `MANIFIESTO_ARBOL.json` / `.sha256` | Árbol más amplio (si existe) |
| `ROOT_HASH.txt` | Resumen de raíz |
| `02_USO_CIVIL.txt` | Límites de uso (obligatorio citar) |

### 6.2 Procedimiento de auditoría de 15 minutos

```text
1. Recibir carpeta / commit acordado
2. python generar_manifiesto.py
3. Comparar MANIFIESTO_HASHES.json con el entregado firmado/publicado
4. Ejecutar:
     - mdc analyze 1147   (o 143)
     - network demo 4 hubs
     - (opcional) check_gemini_y_fisica.py  → doc 03
5. Leer LIMITES_HONESTOS + este 06 + doc 05
6. Rechazar cualquier demo que borre ficheros sin backup (k3_launcher legacy)
```

### 6.3 Qué demuestra (y qué no) el manifiesto

| Demuestra | No demuestra |
|-----------|--------------|
| Bit-a-bit de los archivos listados en la fecha del hash | Que la física de red esté certificada IEC |
| Que no te cambiaron el ZIP en silencio | Que Goldbach está probado |
| Reproducibilidad de demos software | Piloto de turbina real |

---

## 7. Puente cómputo ↔ energía

```text
                    ┌──────────────────┐
                    │  SCADA / gemelo  │
                    │  ZZ / Quijote    │
                    └────────▲─────────┘
                             │ UDP / eventos / Grafcet
                    ┌────────┴─────────┐
                    │  AntiPC L2–L3    │
                    │  hubs, WORK/RES  │
                    └────────▲─────────┘
                             │
              ┌──────────────┼──────────────┐
              │              │              │
         K3_HASH        MDC_PHASE      telemetría
              │              │              │
         integridad    fase/control     auditoría
```

- **ZypyZape en software de red** (demo viability): prueba de arquitectura de hubs, no de aerodinámica.  
- **Kilómetro / tanque ESP32:** el log Serial del 03 puede, en el futuro, entrar como telemetría L4.  
- **XFI:** misma idea de roles gen/thr/buf en un craft simulado (extensión, no producto eólico).

---

## 8. XFI (nota breve — no es el núcleo del 06)

**XFI** (*Experimental Flight Infinite*): gemelo conceptual de craft con **N=3 o 4** turbinas en roles tipo ZypyZape (generación / thruster / buffer).

```bat
cd ...\VMA_mates_rescat_2026\10_inicio_quijote_zypyzape
python gemelo_xfi_avion_4turbinas.py --n 3
python gemelo_xfi_avion_4turbinas.py --compare
```

| Aspecto | Estado |
|---------|--------|
| Física de vuelo real certificable | **No** |
| Laboratorio de control multi-rotor ZZ | **Simulación** |
| Prioridad frente a AntiPC/piloto eólico | **Baja** en venta industrial |

Se lista aquí para que el índice del “1” no olvide la pieza, sin inflar TRL.

---

## 9. Módulos comerciales honestos (software)

Alineado con la lógica de `porque_no_33x1` (módulos separados, no monstruo único):

| Módulo | Contenido | Comprador típico | Dependencia |
|--------|-----------|------------------|-------------|
| **N₁ AntiPC** | Runtime UDP, benchmarks, reuse | Edge / DC local / integrador | Python 3.10+ |
| **N₂ K3** | Hash, dedup, huella | Integridad industrial / telemetría | DLL o pure-Py |
| **N₃ MDC/mates** | Cribas, demos, GUI | Educación / I+D numérica | vma-methods |
| **N₄ Energía** | ZZ / Quijote / Kilómetro | TSO / OEM / lab | docs 01–03 + piloto |

**Regla:** no vender N₄ con promesas de N₁, ni N₂ como AES.

---

## 10. Límites honestos (capa 06)

| Afirmación | Estado |
|------------|--------|
| Benchmarks A→E y viability UDP medibles | **[OK]** en entorno documentado |
| ~49 % menos memcpy usuario en demo slot-ring | **[OK]** captura JSON de referencia |
| AntiPC “libera ALU del universo” | **[NO]** — reduce coordinación/copias |
| K3 criptografía de Estado | **[NO]** |
| Manifiesto = prueba de física de red | **[NO]** — solo integridad de bits |
| Pack monolítico vendible mañana | **[NO]** — módulos N₁–N₄ |
| Listo para SCADA de parque sin OEM | **[NO]** — falta integración y certificación |

---

## 11. Criterios de aceptación (software)

Un piloto software se considera **exitoso** si:

1. Tres máquinas independientes reproducen la **ratio** E/A (o B/A memcpy) ± banda acordada.  
2. Manifiesto regenerado **coincide** con el publicado.  
3. Pipeline `hash → mdc phase` deja traza JSON auditable.  
4. Cero promesas de overunity o de ruptura RSA en materiales de venta.  
5. `02_USO_CIVIL.txt` entregado y aceptado por escrito.

---

## 12. Conclusiones

1. El **06 es prioritario** frente a más refinamiento narrativo de 01–05: es lo **ejecutable y hasheable** hoy.  
2. AntiPC aporta **métricas de red/cómputo**; K3 aporta **huella**; el manifiesto aporta **confianza de entrega**.  
3. Esta capa **no desnucleariza ni mueve turbinas**; **habilita** control, telemetría e integridad del pack civil.  
4. Junto al doc 05, cierra el círculo: *qué física es cierta* + *qué software se puede enseñar mañana*.

---

## 13. Referencias de archivo

| Artefacto | Ruta típica |
|-----------|-------------|
| Arquitectura unificada | `33x1\05-docs-tecnicos\ARQUITECTURA_UNIFICADA.txt` |
| Benchmarks | `...\RESULTADOS_BENCHMARK.txt` |
| Viability JSON | `33x1\antipc\output\zypyzape_viability.json` |
| Fichas E08/E09 | `teoremasgrok\ingenieria\08_ANTIPC.txt`, `09_K3HASH.txt` |
| Comandos | `33x1\repo\33x1\03_COMANDOS.txt` |
| Uso civil | `...\02_USO_CIVIL.txt` |
| XFI | `VMA_mates_rescat_2026\10_inicio_quijote_zypyzape\gemelo_xfi_avion_4turbinas.py` |

---

*Documento 06/06 — Dossier técnico VMA. Solo contenido técnico.*
