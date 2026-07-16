# AntiPC — Adaptive Network Through Parallel Computing

**Autor:** Víctor Manzanares Alberola · **Proyecto:** 33×1 · **CLI:** `0.14.0-cmd` · **DLL:** `0.10.0-c`

> *"Transforming communication into computation."*  
> La red no transporta solo datos: cada paquete UDP es una unidad de cómputo reutilizable.

---

## Contenido

| Ruta | Descripción |
|------|-------------|
| `src/antipc/` | Motor Python completo (CLI, red, MMO, industrial) |
| `src/antipc/antipc_native.dll` | ALU unificada en C |
| `docs/` | Arquitectura científica y resultados |
| `LEEME.txt` | Guía rápida en español |
| `antipc-port-c/` (Escritorio) | Logs de integración v04–v14 |

---

## Fórmula maestra

AntiPC modela el rendimiento útil de un sistema distribuido como:

```
P_util(N) = N · E(N) + K(N)
```

| Símbolo | Significado matemático | Implementación |
|---------|------------------------|----------------|
| **N** | Número de nodos / hubs de cómputo | `hub_node.py`, `UdpFabric`, shards MMO |
| **E(N)** | Eficiencia de coordinación ∈ [0,1] | Pipelining UDP, lotes `WORK_BATCH`, Grafcet |
| **K(N)** | Conocimiento reutilizable (cache hits) | `KnowledgeResolver`, caché hub K3, 2ª pasada sin ALU |

**Ley operativa:** maximizar `P_util` elevando `E(N)` (menos esperas en red) y `K(N)` (mismos payloads → mismo digest sin recalcular).

Benchmark de referencia (arquitectura E, 1000 pkt, 4 hubs): **62 % menos carga ALU** que la arquitectura A convencional.

---

## Hitos matemáticos — mapa completo

### H0 · Fundación K3 (HASHTOOLCODE)

Motor de huella determinista — base de todo el cómputo AntiPC.

| Concepto | Fórmula / propiedad | Módulo |
|----------|---------------------|--------|
| Hash streaming | Rotaciones + XOR + semilla (`k3hash.c`) | `hash_engine.py`, `K3_HASH` |
| Paridad auditable | `hash_C(x) == hash_Python(x)` ∀x | `verify_k3.py`, `antipc k3 verify` |
| Heavy hash | 4 rondas encadenadas K3(payload) | `k3_heavy_hash`, Grafcet, hubs L3 |
| Hashes redundantes | 3 semillas Grafcet: `0xA5A5A5A5`, `0x5A5A5A5A`, `0x12345678` | `antipc k3 redundant` |
| Hamming | `d(h₁, h₂)` sobre representación hex | `antipc k3 hamming` |
| Similitud Jaccard | `\|S₁ ∩ S₂\| / \|S₁ ∪ S₂\|` con shingles K3 | `antipc k3 simil` |
| Índice invertido | término → documentos vía K3 | `antipc k3 search` |
| Motor acordeón 33×1 | `L += v+1`, `v += 2`; validación `5L ≤ 2v+1` | `antipc k3 stream-xor` |

**Hito de integración:** DLL `0.7.0-c` — suite K3 completa en C.

---

### H1 · MDC — Método de Descomposición por Convergencia

Capa L1 de dominio numérico VMA. Herramienta de **exploración** (no factorización RSA industrial).

| Concepto | Expresión | CLI / C |
|----------|-----------|---------|
| Regla mecánica (semiprimo) | `(2x + 3)(2y + 3) = n` | `antipc mdc factor`, `MDC_FACTOR` |
| Dos trenes | Recorrido coordinado en plano (x, y) con colisiones | `antipc mdc analyze`, `mdc_trains.c` |
| Fase modular | Curvatura derivada de `payload mod m` | `MDC_PHASE`, `analisi_modular` |
| Regla log (Battery) | Escala logarítmica de productos | `MDC_REGLA`, `mechanical` |
| K-sweep clásico | Barrido entero de candidatos D | `antipc mdc ksweep` |
| K-sweep predictivo | Predicción de D con menos evaluaciones | `mdc_ksweep.c` (C) |
| **Invariante jerk (MDC-U)** | Pinça 4+4, salto balístico `m−V/A` | `antipc mdc jerk` |

**Demos canónicas:** `n = 1147 = 31×37`, `n = 10473029`, K-sweep `10000319 → D=1949`.

**Hito de integración:** DLL `0.4.0-c` (trenes + criba híbrida), `0.5.0-c` (K-sweep predictivo).

---

### H2 · Teoría de números (vma-methods → C)

Cribas, densidad de primos y logaritmos por convergencia.

| Concepto | Definición / método | AntiPC |
|----------|---------------------|--------|
| Eratóstenes | Criba clásica O(n log log n) | `antipc criba` |
| Criba híbrida VMA | Fusión desmemoriada + segmentada | `antipc criba --hibrida` |
| **Criba desmemoriada** | AND booleano entre listas replicadas | `antipc criba --desmemoriada` |
| Criba modular 6k±1 | Residuo en clase 6k±1 (teoría VMA) | `antipc criba --modular6k` |
| **Criva** | Estimador iterativo de `π(x)/x` | `criva.c`, bridge Python |
| **Newton rápido** | Convergencia con oráculo MEcuation | `antipc newton -f cuadrados` |
| Familias oráculo | `general`, `cuadrados`, `cubos`, `potencia`, `kp`, `mersenne` | `newton_rapido.c` |

**Hito de integración:** DLL `0.3.0-c` (Eratóstenes + Criva), `0.4.0-c` (híbrida), `0.5.0-c` (Newton), `0.10.0-c` (6k±1).

---

### H3 · Convergencia geométrica (encriptacionGeometrica)

Cómputo por perímetros de polígonos y pares de bits — motor reversible con `Decimal`.

| Par bits | Efecto en perímetro | Avance lector |
|----------|---------------------|---------------|
| `10` | + lado[0] | +1 bit |
| `11` | + lado[1] | +1 bit |
| `00` | + lado[0] + lado[1] | +1 bit |
| `01` | sin aporte (paradoja) | +2 bits |

| Variante | Mecanismo | AntiPC |
|----------|-----------|--------|
| Convergencia binaria | Thales + figuras + `Decimal` (prec ≥ 100) | `antipc geo` |
| Aleatorovix | Perímetro × π ofuscado × e convergente | `antipc geo-masivo` |
| K3 geo stream | `(L ^ v) × 0x9E3779B97F4A7C15` → XOR | `k3_geo_aleatorovix.c` |

**Parámetros 33×1:** `base = 33`, `rel = 1`.

**Hito de integración:** DLL `0.6.0-c` — π/e, Aleatorovix, convergencia masiva.

> Uso civil: integridad y educación. No es cifrado frente a adversarios nation-state.

---

### H4 · Grafcet y matriz de existencia (arquitectura D)

| Concepto | Rol matemático | Implementación |
|----------|----------------|----------------|
| Matriz existencia | Evita recomputar digest ya visto | `GrafcetEngine`, `game_engine` |
| Lotes | Amortiza coste de transición de estados | `batch_size` en pipeline B→D |
| Reuse K(N) | Segunda pasada: 0 ALU si referencia resuelta | `KnowledgeResolver`, microkernel |
| Validación heavy | Integridad distribuida por shard | `game_cluster.validate_heavy_batch` |

**Pipeline:** `B (slot-ring recvinto) → D (Grafcet) → WorldState (MMO)`.

---

### H5 · Red como ALU distribuida (arquitectura E, L3)

La comunicación UDP **es** cómputo: cada `WORK` pide un hash; cada `RESULT` devuelve conocimiento.

| Concepto | Modelo | v0.11 |
|----------|--------|-------|
| Protocolo | `WORK` / `RESULT` + HMAC | `udp_protocol.py` |
| Lotes | `WORK_BATCH` (0x0B) / `RESULT_BATCH` (0x0C) | Reduce overhead UDP |
| Pipelining | Ventana `W=128`, lote `B=24` | `PipelinedDispatcher` |
| Throughput | `items / elapsed` | **~2.56×** vs despacho clásico (256 items, 4 hubs) |
| Balanceo | `route(seq, payload) → hub[shard % N]` | Shards MMO, `register_remote_hubs` |
| Caché hub | `K(N)` local por digest K3(payload) | `hub_process_work`, `hub_drain_ring_burst` |

**Comandos:**

```bat
cd src\antipc
python cli.py network bench --count 512 --hubs 4
python cli.py network demo --duration 3 --hubs 4
python cli.py game cluster --shards 4 --hub-hosts IP1,IP2
```

**Hito de integración:** CLI `0.11.0-cmd` — `network_compute.py`, bench reproducible.

---

### H7 · Gemelos digitales — ZypyZape, Quijote, Kilòmetre

Simulación energética del repo VMA enlazada a AntiPC (arquitectura B + modelos físicos).

| Gemelo | Modelo matemático | Enlace AntiPC |
|--------|-------------------|---------------|
| **ZypyZape** | `J·dω/dt = T_aero − T_gen + T_ZZ + T_IS`; swing equation | `zypyzape/slot_ring.py` → `bd_pipeline`, `game_server` |
| **Quijote** | `J_i(r) = J_G + N_b·m_q·r²`; `J·ω̇ = T_net − ω·J̇`; `Cp(λ,β)` NREL 5MW | Gemelo v4.8 / v5 en `repo/quijote/` |
| **Kilòmetre** | `F_n = (ρ_f − ρ_o)·V·g`; `E_k = ½Iω²`; comparativa fluidos densos | `repo/kilometre;(soles_bateria)/` |

**Rutas canónicas en repo:**

```
repo/antipc/zypyzape/          viability_udp.py, hub_emitter.py, slot_ring.py
repo/quijote/                  zypyzape_twin_v4_8_quijote.py, gemelo_v5_quijote_control.py
repo/kilometre;(soles_bateria)/ kilometre_v15_fluids_densos.py
```

**CLI (requiere numpy + matplotlib para Quijote/Kilòmetre):**

```bat
python cli.py gemelo list
python cli.py gemelo info quijote
python cli.py gemelo run zypyzape --hubs 4 --packets 15000
python cli.py gemelo run quijote --variant v48
python cli.py gemelo run quijote --variant v5
python cli.py gemelo run kilometre
```

**Hito de integración:** CLI `0.12.0-cmd` — `gemelos/bridge.py` (puente sin copiar scripts).

> Los gemelos eólicos son **simulación / TRL variable** — no sustituyen banco de pruebas físico (ver `repo/2026/resumen de dinero y de que es solo teorico .txt`).

---

### H8 · Libros 1–6 — métodos VMA → AntiPC

Corpus matemático del autor mapeado a comandos ejecutables.

| Libro | Título | Métodos clave | AntiPC |
|-------|--------|---------------|--------|
| **1** | Números i numeritos | Criba desmemoriada, 6k±1, Salto≈√n, Goldbach asimetría | `criba --desmemoriada`, `--modular6k` |
| **2** | Números oTra VeZ | Dos primos L−m≥2, MRAUV, criba híbrida, espejo | `criba --hibrida`, `mdc analyze` |
| **3** | Sigo en mis Trece | Sofí, Criva, K=9/24, pitagórico visual | `criba`, `mdc visual` |
| **4** | Encriptación y Energía Verde | Convergencia, Aleatorovix, K3, ZypyZape/Quijote/Km | `geo`, `gemelo run` |
| **5** | Factorización (2v+3) | MDC, dos trenes, K-sweep, jerk, numeritos | `mdc factor/analyze/ksweep/jerk` |
| **6** | Pseudocódigo implicaciones | Newton, Método-V, Sofí 15 impl, plantillas 1–7 | `newton`, `teoremas/21–26` |

```bat
python cli.py libro list
python cli.py libro metodos
python cli.py libro info 5
python cli.py libro run 5 --metodo L5-ksweep
python cli.py libro run 4 --metodo L4-convergencia
```

Fuente teoremas: `repo/teoremas/INDICE_MAESTRO.txt` (fichas 01–30).

---

### H9 · DeepSeek 6 2026 + teoremas CLI

Integración del corpus DeepSeek (jun 2026) y lectura de fichas teoremas desde CMD.

| Fuente | Contenido | AntiPC |
|--------|-----------|--------|
| `repo/deepseekjun26/` | Tabla hitos Libros 1–4, líneas matemáticas | `deepseek info libros-tabla` |
| `repo/filestot l5/` | MDC v15–v23, K-sweep predictivo Python | `deepseek run mdc-u -s mdc-v23` |
| `repo/PY L5/` | Scripts DeepSeek jerk/pinza | `mdc jerk` |
| `repo/teoremas/` | INDICE_MAESTRO + fichas 01–30 | `teorema list \| info` |

```bat
python cli.py deepseek list
python cli.py mdc jerk 1147 --factorize
python cli.py teorema info 9
```

**Hito de integración:** CLI `0.14.0-cmd`.

---

### H6 · Proyecto 33×1

Marco auditable que une tratado político [33] y núcleo técnico [1]:

- Motor K3 con `base=33`, `rel=1`
- Manifiesto SHA256 verificable
- Paquete civil según `repo/33x1/02_USO_CIVIL.txt`

---

## Arquitecturas A–E (benchmark)

| ID | Nombre | Idea matemática | Carga ALU ref. |
|----|--------|-----------------|----------------|
| **A** | Convencional | Cola + memcpy (sin reuse) | 100 % |
| **B** | Lock-free | SlotRing SPSC, 1× copia usuario | — |
| **C** | Distribuido | Partición + caché `K(N)` por hub | — |
| **D** | Grafcet | Existencia + lotes + redundancia | — |
| **E** | UDP real | `P_util(N)` en LAN/localhost | **~38 %** (vs A) |

---

## Línea temporal de integración (DLL + CLI)

| Versión | Hitos matemáticos añadidos |
|---------|---------------------------|
| **0.3.0-c** | K3 hash, MDC factor, Eratóstenes, Criva, geo converge, k3_stream |
| **0.4.0-c** | MDC dos trenes (`mdc_trains`), criba híbrida VMA |
| **0.5.0-c** | Newton rápido + oráculo MEcuation, K-sweep clásico/predictivo |
| **0.6.0-c** | Aleatorovix, π/e ofuscados, geo masivo |
| **0.7.0-c** | Suite K3: file, fingerprint, redundant, heavy, hamming |
| **0.10.0-c** | Criba modular 6k±1; k3simil/k3search Python |
| **0.11.0-cmd** | Pipelining UDP, lotes WORK/RESULT, bench red 2.56× |
| **0.12.0-cmd** | Gemelos ZypyZape / Quijote / Kilòmetre — CLI `gemelo` |
| **0.13.0-cmd** | Libros 1–6 métodos VMA — CLI `libro` |

Logs detallados: `Desktop/antipc-port-c/02_INTEGRADO_v04.txt` … `11_INTEGRADO_v11.txt`.

---

## CLI — referencia por área matemática

```bat
cd src\antipc

:: Núcleo y versión
python cli.py version
python cli.py native status
python cli.py native bench

:: K3 / HASHTOOLCODE
python cli.py k3 hash --text "Hola"
python cli.py k3 verify
python cli.py k3 heavy --text "payload"
python cli.py k3 simil --dir docs --threshold 0.30
python cli.py k3 stream-xor --text "33x1"

:: MDC
python cli.py mdc factor 1147
python cli.py mdc analyze --n 1147
python cli.py mdc ksweep 10000319
python cli.py mdc visual --gui

:: Teoría de números
python cli.py criba --limit 50000 --modular6k
python cli.py newton 121 -f cuadrados

:: Geometría / convergencia
python cli.py geo --demo
python cli.py geo-masivo --bits 256

:: Red y cómputo distribuido
python cli.py network bench
python cli.py network demo
python cli.py game cluster --shards 4

:: Gemelos energéticos (repo VMA)
python cli.py gemelo list
python cli.py gemelo run zypyzape
python cli.py gemelo run quijote --variant v5

:: Libros 1–6 (métodos matemáticos)
python cli.py libro list
python cli.py libro metodos --libro 5
python cli.py libro run 5 --metodo L5-jerk

:: DeepSeek 6 2026 + teoremas
python cli.py deepseek list
python cli.py deepseek run mdc-u -s mdc-v23 --n 1147
python cli.py teorema list
python cli.py teorema info 9

:: Industrial
python cli.py industrial inventory
python cli.py industrial audit --dir C:\datos
```

---

## Capas del stack (`AntiPCStack`)

```
L4  Industrial   — telemetría CSV/JSON, auditoría ficheros
L3  Red          — UdpFabric, hubs, pipelining, HMAC
L2  Flow Kernel  — referencias, resolver, plugins UPS
L1  Dominio      — MDC_*, K3_*, NetworkAuthGate
L0  Fundación    — k3hash.c / antipc_native.dll
```

---

## Requisitos y arranque

- Windows 10/11 (o Linux/macOS con Python 3.10+)
- Python 3.10+, solo biblioteca estándar
- `antipc_native.dll` junto a `cli.py` (compilar con CMake + VS si falta)

```bat
scripts\01_instalar.bat
scripts\21_build_antipc_native.bat
cd src\antipc && python cli.py version
```

---

## Honestidad técnica

| Ámbito | Límite declarado |
|--------|------------------|
| K3 | Huella / dedup / índice — **no** contraseñas ni firma legal |
| MDC | Factorización **toy** ≤ 10 dígitos |
| Geo / Aleatorovix | Motor determinista educativo — no AES competidor |
| Red UDP | Laboratorio LAN — sin TLS; no WAN producción |
| Newton VMA | Puede no converger en 200 iter (p. ej. E=121) — comportamiento fiel C/Python |

---

## Soporte y documentación

| Nivel | Recurso |
|-------|---------|
| Arquitectura | `docs/ARQUITECTURA_UNIFICADA.txt` |
| Resultados bench | `docs/RESULTADOS_BENCHMARK.txt` |
| Inventario vivo | `python cli.py industrial inventory` |
| Integración C | `Desktop/antipc-port-c/LEEME.txt` |

**VMA — AntiPC / 33×1** · El Palomar, 2026