# Índice técnico — repositorio `claude` (VMA / 33×1)

**Autor:** Víctor Manzanares Alberola · **Actualizado:** 2026-07-16

Este índice ordena los paquetes ejecutables y el corpus matemático. El resto de carpetas históricas se conservan como archivo.

---

## Paquetes principales (usar estos)

| Carpeta | Versión / estado | Comandos / uso |
|---------|------------------|----------------|
| **antipc/** | CLI `0.14.0-cmd` · DLL `0.10.0-c` | `cd antipc/src/antipc && python cli.py version` |
| **vma-methods/** | Cribas, Newton, Criva | `vma-methods.cmd` o `python -m vma_methods` |
| **teoremas/** | Fichas 01–30 + formal | `antipc teorema list` (vía antipc) |
| **33x1/** | Manifiesto, uso civil, comandos | `01_INDICE_TECNICO.txt` |
| **docs/integracion-antipc/** | Logs v04–v14 port C | `14_INTEGRADO_v14.txt` |
| **encriptacionGeometrica/** | Convergencia binaria, K3, Aleatorovix | `_demo_convergencia.py`, `k3_core.c` |
| **aleatorovix/** | Organismo v1.0 (entropía + Lila + criba) | `python aleatorovix.py` |
| **docs/aleatorovix-firma-geometrica/** | Hub firma geométrica + Aleatorovix | `INDICE_ALEATOROVIX_GEO.txt` |
| **vma/** | Paquete K3 cross-verifier | `vma/k3/` |
| **graficas y explicaciones/** | Gráficas MDC, MRAUV, AntiPC | `00_INDICE.txt` |

---

## Libros y MDC (corpus)

| Carpeta | Contenido |
|---------|-----------|
| `Libro1 Números i numeritos/` … `Libro6 NewtonRapido…/` | PDFs y fuentes por libro |
| `Libro5 Factorizacion con 2v+3/` | Corpus Libro 5 completo |
| `filestot l5/` · `PY L5/` | MDC v15–v23, K-sweep, DeepSeek Python |
| `deepseekjun26/` | Chats y tablas DeepSeek jun 2026 |
| `lee arbusto/` | Libro 4 convergencias (completo) |
| `teoremas/formal/` · `spacial/` · `ingenieria/` | Extensiones teoremas vma-github-2026 |
| `ideas-para-gpt-antipc/` | Roadmaps AntiPC históricos |
| `just run/` | Scripts gemelos ejecutables |

---

## Gemelos energéticos

| Carpeta | Sistema |
|---------|---------|
| `zypyzape-contexto/` · `ZYPYZAPE Bateria Cinetica/` | ZypyZape |
| `Quijote/` · `Quijotee/` · `zypyzape quijote ballant/` | Quijote |
| `kilometre;(soles_bateria)/` | Kilòmetre |

Ejecutables vía AntiPC: `antipc gemelo list | run`

---

## Hash / K3

| Carpeta | Nota |
|---------|------|
| `vma-k3/` · `hashtool-work/` | Motor K3 y herramientas |
| `antipc/src/antipc/` | Suite integrada `antipc k3` |

---

## Aleatorovix + firma geométrica

| Ruta | Contenido |
|------|-----------|
| `docs/aleatorovix-firma-geometrica/` | Índice, README, fuentes C auditables |
| `encriptacionGeometrica/` | Teoría Tales, demo Decimal, K3 industrial |
| `aleatorovix/` | Organismo canónico (espejo en `just run/aleatorovix/`) |
| `antipc/.../geo_convergence.c` | Convergencia binaria en C |
| `antipc/.../k3_geo_aleatorovix.c` | π/e ofuscado + cifrado masivo |
| `teoremas/19_teorema_phaseamplifier_k3.txt` | Teorema fase K3 XOR |

```bat
cd antipc\src\antipc
python cli.py geo --bits 0101101011000010
python cli.py geo-masivo --text "hola vma"
python cli.py k3 stream-xor --text "libro4-33x1"
python cli.py libro run L4-convergencia
```

---

## Histórico / duplicados (no borrar sin revisar)

- `antipc2/` — snapshot anterior; **canónico = `antipc/`**
- `grok/` · `CLAUDE TXT/` · `copilot 2025/` — chats exportados
- `UNION/` · `predecir log raiz/` — variantes Newton / log

---

## Flujo de actualización desde Escritorio

```powershell
# 1. Sincronizar paquetes
robocopy C:\Users\cuent\Desktop\repo\antipc C:\Users\cuent\Desktop\repo\_clone_claude\antipc /E /XD .git __pycache__ build
robocopy C:\Users\cuent\Desktop\repo\vma-methods C:\Users\cuent\Desktop\repo\_clone_claude\vma-methods /E /XD __pycache__
robocopy C:\Users\cuent\Desktop\encriptacionGeometrica C:\Users\cuent\Desktop\repo\_clone_claude\encriptacionGeometrica /E /XD __pycache__
robocopy "C:\Users\cuent\Desktop\repo\just run\aleatorovix" C:\Users\cuent\Desktop\repo\_clone_claude\aleatorovix /E /XD __pycache__

# 2. Commit + push
cd C:\Users\cuent\Desktop\repo\_clone_claude
git add -A
git commit -m "Actualización VMA"
git push origin main
```

---

## Enlaces

- Repo: https://github.com/espiradesombra/claude
- AntiPC: https://github.com/espiradesombra/claude/tree/main/antipc