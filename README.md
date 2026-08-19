# espiradesombra / claude

Monorepo **VMA** (Víctor Manzanares Alberola): mates, energía (ZypyZape / Quijote / Kilómetro), AntiPC / K3, 33×1, webs techamv / Aleatorovix.

**Repo:** https://github.com/espiradesombra/claude

---

## Estado del arte (2026-08)

### Trato 33×1

> **[1]** = este repositorio técnico **civil** (auditable)  
> **[33]** = **33 años de paz firmada** por los países

| Doc | Contenido |
|------|-----------|
| [`33x1/00_QUE_ES_33x1.md`](33x1/00_QUE_ES_33x1.md) | Definición del trato |
| [`33x1/06_FOR_THE_WORLD_33x1.md`](33x1/06_FOR_THE_WORLD_33x1.md) | Mensaje público (America/world, Tesla·xAI, inercia civil) |
| [`33x1/02_USO_CIVIL.txt`](33x1/02_USO_CIVIL.txt) | Uso civil obligatorio |
| [`MAPA_PC_AL_GIT.md`](MAPA_PC_AL_GIT.md) | Qué se ha traído del Desktop al git |

### Energía — física y gemelos

| Pieza | Estado | Dónde |
|-------|--------|--------|
| **ZypyZape** | Buffer cinético multi-turbina; sims + docs | `ZYPYZAPE Bateria Cinetica/`, `zypyzape-contexto/`, `docs/tecnicos-vma/01_*` |
| **Quijote** | Masas en pala / J(r); 3 vs 7; coste actuador | `Quijote/`, `hurto-gravitatorio/`, `docs/gemelo_quijote_3vs7_*` |
| **Kilómetro** | Máquina kilo+metro; sims + ESP32 | `kilometro_sim/`, `kilometre;(soles_bateria)/`, `docs/gemelo_kilometro_minimo.py` |
| **Explicación unificada** | Texto + diagrama 1 página | [`docs/EXPLICACION_FISICA_ZYPYZAPE_QUIJOTE_KILOMETRO.md`](docs/EXPLICACION_FISICA_ZYPYZAPE_QUIJOTE_KILOMETRO.md) |
| **Honestidad / IAs** | Hype vs trellat | [`docs/DE_ON_VE_TOT_AIXO.md`](docs/DE_ON_VE_TOT_AIXO.md) |

**Estado:** principios de inercia variable y red **sólidos en simulación**; prototipo físico / piloto de parque **pendiente**. Actuador y η&lt;1 siempre cuentan. No es motor de gravedad perpetuo.

### Matemáticas VMA

| Línea | Dónde |
|-------|--------|
| Cribas (desmemoriada, modular 6k±1, híbrida) | `vma-methods/`, `VMA_mates_rescat_2026/01_cribas/` |
| Criva / MRAUV | `VMA_mates_rescat_2026/02_criva_mrauv/`, `Metodo densidad MRAUV/` |
| MDC (diofántico cinemático) | `VMA_mates_rescat_2026/03_mdc/`, `PY L5/`, AntipC plugins |
| Newton Rápido | `Libro6…`, `Metodo Newton Rápido/`, `predecir log raiz/` |
| Goldbach / Sofí / Fermat | `Goldbach/`, `Sophie Germain/`, pack rescat |
| Teoremas | `teoremas/`, `teoremasgrok/` |

**Estado:** corpus amplio + fixes (fase criba 2p/4p, `siguiente_primo` wheel). Parte toy / exploratoria; otra parte documentada en libros 1–6.

### Software industrial civil

| Pieza | Dónde |
|-------|--------|
| AntiPC (runtime / UDP / KOP) | `antipc/`, `antipc-port-c/`, `inbox/carpetas-de-escritorio-select/06-python-runtime/` |
| K3 / hash / integridad | `vma-k3/`, `encriptacionGeometrica/`, `hashtool-*` |
| **EliminaDuplicadosK3** (GUI + exe) | `hashtool-work/bin/EliminaDuplicadosK3.exe` + `elimina_duplicados_k3.py` |
| Sale-it / packs | `sale-it/`, `inbox/subir/` |

### Webs

| Site | Dónde |
|------|--------|
| **Aleatorovix** (ramo, Leonor/Princesa, leona favicon) | `web-aleatorovix/`, `techamv-aleatorovix-DEPLOY/` |
| **techamv / ZypyZape landing** | `techamv-web/` (index completo + `404.html` easter eggs) |

### Literatura / miscelánea reciente

| | |
|--|--|
| Muñeco de Nieve — órbitas en espiral | [`docs/muneco-nieve/`](docs/muneco-nieve/) |
| Inbox Desktop (docx/txt VMA) | [`docs/inbox-desktop-2026-08/`](docs/inbox-desktop-2026-08/) |
| **Rescate montón viejo** (selectivo) |
| Desconocidos montón (Apiñón, CV…) | [`docs/from-monton-desconocidos-2026-08/`](docs/from-monton-desconocidos-2026-08/) | [`docs/from-monton-viejo-2026-08/`](docs/from-monton-viejo-2026-08/) |
| Staging audit PC | [`pc-missing-2026-08-20/`](pc-missing-2026-08-20/) |

---

## Mapa rápido del monorepo

```
claude/
├── 33x1/                    ← trato político-técnico + inbox cartas
├── docs/                    ← física, gemelos, técnicos-vma, muneco, from-monton-viejo
├── kilometro_sim/           ← sims Kilómetro (Desktop/33x1)
├── kilometre;(soles_bateria)/
├── VMA_mates_rescat_2026/   ← pack mates + XFI rescatado
├── vma-methods/ · antipc/ · vma-k3/ · encriptacionGeometrica/
├── Quijote/ · ZYPYZAPE…/ · hurto-gravitatorio/ · zypyzape-contexto/
├── web-aleatorovix/ · techamv-web/ · techamv-aleatorovix-DEPLOY/
├── Libro1…Libro6/ · teoremas/ · Goldbach/ …
├── inbox/                   ← subir + carpetas-de-escritorio select
├── pc-missing-2026-08-20/   ← staging/audit del inventari PC
└── MAPA_PC_AL_GIT.md        ← origen Desktop → destino git
```

Índice por carpeta (histórico): [`FOLDERS.md`](FOLDERS.md)

---

## Cómo empezar

1. **33×1:** lee `33x1/00_QUE_ES_33x1.md`  
2. **Energía:** `docs/tecnicos-vma/00_INDICE.md` + explicación física  
3. **Mates:** `VMA_mates_rescat_2026/` o `vma-methods/`  
4. **Web ramo:** `web-aleatorovix/index.html` (easter eggs: escribe `leonor` / princesa…)  
5. **Inventari PC:** `MAPA_PC_AL_GIT.md`

### Gemelos toy (docs/)

```bat
cd docs
python gemelo_kilometro_minimo.py --mode drain --T 25
python gemelo_quijote_3vs7_balance.py --mode phase --T 25
python gemelo_grupo_zz3.py --compare-H
```

---

## Avisos

- **Uso civil.** No militar. Ver `33x1/02_USO_CIVIL.txt`.  
- Monorepo **histórico y grande**: hay carpetas duplicadas / chats / WIP. La “fuente de verdad” reciente está en `docs/`, `33x1/`, `kilometro_sim/`, `VMA_mates_rescat_2026/`, `web-aleatorovix/`.  
- No versionar instal·ladors, cracks ni vídeos personals.  
- Si el repo se ve “estropeado” por herramientas locales: prioriza Drive/backup que indiqui l’autor + aquest mapa.

**1477 · 33×1 · make the world great forever** — not by war, by phase alignment and shared knowledge.
