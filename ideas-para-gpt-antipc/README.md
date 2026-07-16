# JUST RUN — Ecosistema VMA / 33×1

**Paquete ejecutable y presentable** del trabajo de I+D de Víctor Manzanares Alberola (EPSA, UPV Alcoi).

Selección curada de matemáticas (MDC, Criva, cribas), energía (Quijote, ZypyZape, hurto gravitatorio), runtime industrial (K3, AntiPC) y marco 33×1. Los originales en `repo\`, `vma\` y `lee arbusto\` **no se modifican** — esto es una copia reproducible lista para ejecutar, auditar o presentar.

| | |
|---|---|
| **Versión paquete** | 2026-07-11 |
| **Archivos** | ~211 (ver `MANIFEST.txt`) |
| **TRL orientativo** | K3 / AntiPC 3–5 · Gemelos simulación 4 · Hurto gravitatorio sim |
| **Idioma docs** | Español (código: ES/EN mixto) |

---

## Arranque en 10 segundos

**Windows — doble clic:**

```
RUN.bat
```

**Terminal:**

```bat
cd "just run"
RUN.bat demo
python vma-run\vma_run.py demo
```

La demo unifica en &lt;1 s las seis piezas del ecosistema:

1. K3 — auditoría industrial  
2. MDC — factorización cinemática  
3. Criva — densidad de primos  
4. Criba modular 6k±1  
5. Amplificador de fase K=3 (XOR)  
6. Hurto gravitatorio — tabla NREL 5 MW  

**Comandos útiles:**

```bat
RUN.bat mdc 10403
RUN.bat k3 --text "hola"
RUN.bat criva --compare 100,1000,10000
RUN.bat hurto
RUN.bat antipc
```

---

## Requisitos

| Componente | Requisito |
|------------|-----------|
| Python | 3.10+ |
| Dependencias base | `pip install -r requirements.txt` → numpy, matplotlib |
| vma-k3 | `cd vma-k3 && pip install -e .` |
| AntiPC | `cd antipc && scripts\01_instalar.bat` |
| Aleatorovix (C) | WSL o MinGW + gcc |
| Libro PDF | `libro-metodos-simples\COMPILAR.bat` (pdflatex) |
| Ejecutable opcional | `bin\vma-run.exe` si está compilado |

---

## Contenido del paquete

| Carpeta | Qué es | Empezar por |
|---------|--------|-------------|
| `vma-run/` | Ejecutable unificado + `lib/` | `RUN.bat`, `CONCEPTO.txt` |
| `gemelos/` | Gemelos digitales Quijote v3–v10 | `gemelo_v6_fisic.py` |
| `hurto-gravitatorio/` | Simulación hurto NREL 5 MW | `gemell_quijote_paper_rules.py` |
| `zypyzape-contexto/` | Contexto, matemática 3 vs 7, dossier Elewit | `01_CONTEXT_ZYPYZAPE.txt` |
| `vma-k3/` | Motor industrial Theorem K3 | `pip install -e .` → `vma-k3` |
| `antipc/` | Runtime AntiPC + UDP benchmark | `INICIO.bat` |
| `gptcomputing/` | Spec fundacional AntiPC (chat → código) | `CONCEPTO.txt`, `MAPA_*.txt` |
| `archivos-vma/` | Corpus matemático libro + ForPkm | `cribas_cotas_vma.txt` |
| `aleatorovix/` | Organismo Lila, criba desmemoriada | `RUN_ALEATOROVIX.bat` |
| `teoremas/` | Teorema K3 XOR, criva.py | `THEOREM_*.md` |
| `33x1/` | Marco científico y pacto territorial | `00_PROMPT_DEFINITIVO_33x1.md` |
| `libro-metodos-simples/` | LoveArt txt + LaTeX LoveEarthHacker | `LOVEART_fisica_o_informatica.txt` |
| `cpu-antipc/` | Contexto C++ AntiPC (extracto) | `CONCEPTO.txt` |

Documentación extendida en español: **`LEEME.txt`**  
Índice maestro (conjeturas, teoremas, métodos, mercados…): **`TABLA_CONTENIDO_VMA.md`**

---

## Ruta de demostración (~15 min)

| Paso | Acción | Tiempo |
|------|--------|--------|
| 0 | `RUN.bat` — demo vma-run | 10 s |
| 1 | Leer `archivos-vma\cribas_cotas_vma.txt` + contexto ZypyZape | 5 min |
| 2 | `hurto-gravitatorio\` → regenerar CSV/LaTeX | 2 min |
| 3 | `gemelos\gemelo_v6_fisic.py` → gráficos | 1 min |
| 4 | `vma-k3` modo industrial | 30 s |
| 5 | `aleatorovix\` criba + MDC v6 | 30 s |
| 6 | `antipc\INICIO.bat` (si hay red local) | opcional |

Scripts de acceso rápido en raíz: `RUN_GEMELOS.bat`, `RUN_ALEATOROVIX.bat`, `RUN_ARCHIVOS_VMA.bat`, `RUN_ANTIPC.bat`

---

## Resultados clave (simulación, reproducibles)

**Hurto gravitatorio** (`hurto-gravitatorio/quijote_results.csv`):

| Modo | P_net | η_hurto | ΔP_grid |
|------|-------|---------|---------|
| 3 pales PAPER | +9.0 kW | 3.0× | +4.16% |
| 7 pales PAPER | +31.1 kW | 3.7× | +4.31% |

Modo CONTROL (centrífuga práctica): aún negativo — brecha documentada como pendiente.

**AntiPC v0.1.0** (knowledge + UDP): ~31% ALU ahorrada con 35% duplicados (ver `gptcomputing/MAPA_gptcomputing_a_codigo.txt`).

---

## Valor y uso previsto

- **Due diligence técnica** — gemelos + CSV + contexto numérico  
- **Pitch institucional** — dossier Elewit en `zypyzape-contexto/dossier-elewit/`  
- **I+D deeptech** — K3, AntiPC, Criva como líneas B2B  
- **Marco 33×1** — pacto territorial y whitepaper 1477  

Estimación conservadora documentada en `LEEME.txt` §4: 500k–3M € opción I+D; K3 piloto 20k–150k €/año posible.

---

## Qué NO incluye este paquete

- Originales completos del escritorio (`Desktop\archivos` ~98k archivos)  
- Carpeta `UPV EPSA` ni material didáctico masivo  
- `_entrada` pesada de vma-k3 (vídeos/zips)  
- AntiPC definitivo C++ (pendiente integración)  
- Piloto hardware 2 turbinas ni certificación IEC  
- Contrato 33×1 redactado (borrador pendiente)  

Próximas piezas identificadas: Goldbach (55 archivos), `kilometre;batería` (38), extracción Copilot.7z — ver `archivos-vma/ANALISIS_FASE2.txt`.

---

## Origen de las copias

| Destino `just run/` | Origen |
|---------------------|--------|
| `gemelos/` | `repo\Quijotee\` |
| `hurto-gravitatorio/` | `repo\zypyzape quijote ballant (0)\` |
| `zypyzape-contexto/` | `claude-main-extract\claude-main\2026\` |
| `vma-k3/` | `vma\` |
| `antipc/` | `lee arbusto\sale-it\` |
| `gptcomputing/` | `lee arbusto\gptcomputing (2).txt` |
| `archivos-vma/` | `Desktop\archivos\` (curado) |
| `aleatorovix/` | `repo\txt l5\` |

---

## Zip unificado

```bat
EMPACAR_ZIP.bat
```

Genera `just-run-unificado.zip`. Valoración honesta: **`VALOR_ZIP_UNIFICADO.txt`** (80k–350k € opción deeptech hoy).

**Continuar AntiPC en ChatGPT:** `gptcomputing/CONTINUACION_PARA_GPT.txt` + adjuntos listados → Sprint 2 (Ledger, Policy, Metrics).

---

## Checklist — próxima subida

Usar esta lista antes de publicar en GitHub, SharePoint, Elewit o envío a terceros.

### Preparar

- [ ] Comprobar `RUN.bat demo` en máquina limpia (Python 3.10+)
- [ ] Revisar que no hay contraseñas, claves MEGA ni `Contraseñas Chrome.csv`
- [ ] Actualizar `MANIFEST.txt` si se añaden archivos
- [ ] Opcional: recompilar `bin\vma-run.exe` si cambió `vma_run.py`
- [ ] Generar ZIP sin `__pycache__`, `.egg-info`, `build/`, `*.pyc`

### Incluir en el mensaje de subida

- [ ] Enlace o adjunto: este `README.md`
- [ ] Autor: Víctor Manzanares Alberola — EPSA UPV Alcoi
- [ ] Aviso: simulación ≠ hardware certificado; hurto CONTROL pendiente
- [ ] Punto de entrada: `RUN.bat` + `LEEME.txt`
- [ ] Para REE/Elewit: `zypyzape-contexto/dossier-elewit/`

### Destinos sugeridos

| Destino | Qué subir primero |
|---------|-------------------|
| **GitHub privado** | Repo completo `just run` (~2–3 MB sin zips) |
| **Elewit / REE** | `zypyzape-contexto/dossier-elewit/` + `hurto-gravitatorio/quijote_results.csv` |
| **Piloto K3 B2B** | `vma-k3/` + `vma-run` demo K3 |
| **Paper matemático** | `archivos-vma/cribas_cotas_vma.txt` + `teoremas/` |
| **SharePoint UPV** | `README.md` + `LEEME.txt` + demo vídeo de `RUN.bat` |

### Versión sugerida para el commit/tag

```
just-run-2026.07.11
```

Mensaje de release corto:

> Paquete ejecutable VMA: demo unificada, gemelos Quijote, hurto gravitatorio (+9/+31 kW sim), K3, AntiPC v0.1.0, corpus matemático archivos-vma, marco 33×1.

---

## Licencia y propiedad intelectual

**Propiedad intelectual:** Víctor Manzanares Alberola.

El plan comercial y el pacto **33×1** contemplan licencia territorial perpetua para adherentes. Este paquete es una **copia de trabajo** para ejecución y presentación; no transfiere derechos automáticamente.

Contacto académico: EPSA — Universitat Politècnica de València (campus Alcoy).

---

## Changelog del paquete

| Fecha | Cambio |
|-------|--------|
| 2026-07-11 | Creación `just run`: vma-run, gemelos, hurto, K3, AntiPC, 33×1 |
| 2026-07-11 | + `gptcomputing/`, `aleatorovix/`, `cpu-antipc/`, libro LaTeX |
| 2026-07-11 | + `archivos-vma/` fase 1–2 (32 archivos, `ANALISIS_FASE2.txt`) |
| 2026-07-11 | + `README.md` para próxima subida |

---

*just run — copia lo valioso, ejecuta, presenta.*