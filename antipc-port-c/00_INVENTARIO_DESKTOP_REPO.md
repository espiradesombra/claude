# Inventario Desktop + Repo — qué hay y qué ya está en C

**Fecha:** 2026-07-16 · **Regla:** no mover archivos; solo referencias.

---

## Ya en C (repo — `antipc_native.dll` v0.3.0-c)

| Módulo | Fuente Python/original | Archivo C repo |
|--------|------------------------|----------------|
| K3 hash | `k3hash.c` (HASHTOOLCODE) | `native/k3hash/src/k3hash.c` |
| MDC factor toy | `mdc_lib/factoritzacio.py` | `native/antipc_core/src/mdc_factor.c` |
| Criba Eratóstenes | `vma-methods/classic.py` | `native/antipc_core/src/sieve.c` |
| Criva básica | `vma-methods/criva.py` | `native/antipc_core/src/criva_est.c` |
| Stream K3 33×1 | `encriptacionGeometrica/k3_core.c` | `native/antipc_core/src/k3_stream.c` |
| Geo convergencia | `encriptacionGeometrica/_demo_convergencia.py` | `native/antipc_core/src/geo_convergence.c` |

**DLL desplegada:** `repo\antipc\src\antipc\antipc_native.dll`

---

## Escritorio — carpetas nuevas / activas

| Ruta | Contenido | Estado C |
|------|-----------|----------|
| `encriptacionGeometrica\` | 51 archivos: teoría, `k3_core.c`, borradores DeepSeek/Gemini | Parcial: k3_stream + geo en DLL |
| `encriptacionGeometrica\deepseek_c_*.c` | π Taylor, motor K3 long double | **No integrado** — candidato fase 2 |
| `encriptacionGeometrica\gemini-code-*.py` | Aleatorovix, backtracking Decimal | Python; C posible (double) |
| `.grok\skills\` | encriptacion-geometrica, vma-*, modo-canto | Solo guías agente |
| `monton viejo revisar\gemelo_*.py` | Gemelos aerogenerador (Kuramoto, Cp) | Python/NumPy — C opcional fase 3 |
| `archivos\` | fermat_modular, estructura_sofi, etc. | Revisar caso a caso |

---

## Repo — núcleos matemáticos VMA

| Ruta | Módulos | ¿Pasar a C? |
|------|---------|-------------|
| `repo\vma-methods\` | cribas (3), criva, newton, classic | **Alta** cribas + newton; criva parcial hecho |
| `repo\antipc\src\antipc\` | CLI, Grafcet, UDP, microkernel | **Media** hot loops; **baja** orquestación |
| `repo\antipc\mdc_lib\` | analyze, visual, factoritzacio | **Alta** analyze trains; **baja** GUI/GIF |
| `repo\antipc\referencia\mdc\*.hpp` | C++ header-only MDC | Ya diseñado para C++; port directo |
| `repo\filestot l5\ksweep_predictiu.py` | K-sweep entero predictivo | **Alta** — puro entero, N grande |
| `repo\filesclaude 6-5\` | mdc, criva, cribas, discriminant | Duplicados vma-methods / MDC |
| `repo\goldbach\`, `repo\diamante\` | Goldbach, heurísticas | Media-baja (exploratorio) |
| `repo\33x1\` | Manifiesto, docs | No ALU — queda Python/batch |

---

## C existente disperso (no unificado)

- `repo\vma-k3\c\`, `repo\hashtool-extract\`, `_zip_extract\grok\vma\c\` — copias k3hash/k3_audit
- `repo\just run\aleatorovix\organismo_lila_v99.c` — organismo caótico
- `encriptacionGeometrica\` — 8+ `.c` borradores Jul 2026

**Recomendación:** unificar en `antipc_native` v0.4+ sin mover originales; copiar lógica desde este escritorio.