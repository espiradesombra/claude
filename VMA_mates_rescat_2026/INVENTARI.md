# VMA mates — pack de rescat 2026-07-23

**Origen:** `C:\Users\cuent\Desktop\monton viejo revisar`  
**Destí:** `C:\Users\cuent\Desktop\VMA_mates_rescat_2026`  
**Autor material:** Víctor Manzanares Alberola (VMA)  
**Rescat fet per:** Grok (xAI) en revisió 2ª del PC

Aquest pack **no esborra res del montón original**. Només copia el que és recuperable i útil (codi + documents canònics), sense instal·ladors, APKs ni basura.

---

## Qui ha de mirar què primer

| Prioritat | Carpeta | Contingut |
|-----------|---------|-----------|
| **P0** | `00_sunraman_eratostenes/` | Pseudocodi **sunraman** (criba 6k±1) + docs Eratòstenes |
| **P0** | `01_cribas/` | `cribas.py` unificat + modular/híbrida + OpenMP C |
| **P0** | `03_mdc/` | MDC base + línia `mdc_v18…v23` + híbrids |
| **P1** | `02_criva_mrauv/` | Criva densitat + MRAUV + salts |
| **P1** | `04_newton_rapido/` | Newton Rápido + papers |
| **P1** | `05_sofi_fermat_goldbach/` | Sofí, Fermat, Goldbach v8–v10 |
| **P1** | `07_documents_canonics/` | Tractats VMA + Números i numeritos |
| **P2** | `06_upv_c_historics/` | Experiments C UPV (desmemoriada, híbrida, OpenMP) |
| **P2** | `08_scripts_anexos_VMA/` | Scripts curts d’annexos (còpia de treball) |
| **P2** | `09_gemelos_zypyzape_ultims/` | Gemelos v9/v10 + K3 hub UI |

---

## 00 — sunraman / Eratòstenes

| Fitxer | Notes |
|--------|--------|
| `sunraman.docx` | Original Word amb pseudocodi |
| `sunraman_pseudocodi.txt` | Export text llegible |
| `Criba de Erastoteles.docx` | Notes clàssiques |
| `Eratostenes y jo.docx` | Text personal |
| `monica erast.docx` | Variant |
| `ventajas del metodo hibrido contra erastotes.docx` | Híbrida vs clàssica |

**Fórmules clau (sunraman):**
- primos = `2n+3`
- compuestos = `4mn + 6m + 6n + 9`
- marcatge amb `anterior` / `anteriorM` / `anteriorNY`

**Línia evolutiva:**  
`sunraman` → `criba_modular.py` → `anexoF_criba6kpm1_openmp.c` → `cribas.py` (mòdul complet)

---

## 01 — Cribas (codi)

| Fitxer | Notes |
|--------|--------|
| `cribas.py` | Desmemoriada + modular + híbrida (documentat, millor versió) |
| `criba_modular.py` | Script curt annex |
| `criba_hibrida.py` | Script curt annex |
| `anexoF_criba6kpm1_openmp.c` | C + OpenMP |
| `criba_desmemoriada.txt` | Guion / explainer |
| `Cribas_cotas_y_estructuras_de_primos_VMA.txt` | Text tècnic complet |

---

## 02 — Criva / MRAUV

| Fitxer | Notes |
|--------|--------|
| `criva.py` | Estimador densitat π(x)/x |
| `criva.txt` | Notes llargues |
| `mrauv.py`, `anexoE_L_m_script.py`, `anexoF_mrauv_calibrador.py` | L(n), m(n), calibratge |
| `mrauv_goldbach.py`, `salto_maximo.py` | Goldbach / salts |
| `heuristica_*.py`, `restos_graficos.py` | Heurístiques |
| `antiprimos.txt`, `metodo_calculo_primos.txt` | Notes |

---

## 03 — MDC (Mètode Diofàntic Cinemàtic)

| Fitxer | Notes |
|--------|--------|
| `mdc.py` | Motor unificat documentat |
| `mdc_v18.py` … `mdc_v23.py` | Evolució (preferir **v23**) |
| `mdc_hybrid_v444.py` | Híbrid |
| `mdc_factorizador.py` | App factorització |
| `mdc_benchmark.py`, `benchmark_frontera.py`, `record_mundial.py` | Benchmarks |
| `mdc_fusio_v3.py`, `mdc_libro5.py`, `mdc_estudi_senyal.py` | Variants |

---

## 04 — Newton Rápido

| Fitxer | Notes |
|--------|--------|
| `newton_rapido.py` | Algorisme + oracles MEcuation |
| `paper1_newton_rapido.docx`, `newton_rapido_preprint.docx` | Papers |
| `me_detector.py` | Detector auxiliar |
| `tu_algo.cpp` | Si present (implementació C++) |

---

## 05 — Sofí / Fermat / Goldbach / altres

| Fitxer | Notes |
|--------|--------|
| `sofi_structure.py`, `estructura_sofi_mini.py` | Sophie Germain / L1–L4 |
| `fermat_modular.py` (+ mini) | Alineació modular Fermat |
| `siguiente_primo.py` (+ mini) | Heurística següent prim |
| `Goldbach_Articulo_VMA_v8/v9/v10.docx` | Articles (preferir **v10**) |
| `conjetura_goldbach.xlsx`, `afirmaciones_goldbach.xlsx` | Dades |
| `riemann_deformado.py`, `discriminant.py`, `metodo_restos.py` | Experiments |

---

## 06 — UPV C històrics

Origen: `UPV EPSA\...\computacion paralela\primos\`

- `AlgDesMemoriado.c` — criba desmemoriada
- `9.0metodoHibrido.c` — mètode híbrid
- sèrie `8.*`, `10.*` Comparativa / OpenMP
- `prime_sieve.c`, `fast_sieve.c` — referència Eratòstenes segmentat (extern/clàssic)

---

## 07 — Documents canònics

- `VMA_Primos_Completo_v3.docx` / `v3_1`
- `Cribas_cotas_y_estructuras_de_primos_VMA.docx`
- `tratado_criba_y_primos.docx`
- `Numeros_i_numeritos.pdf`
- Goldbach ja a carpeta 05
- Mersenne, primos seguros, obra científica, col·lecció unificada, Zypyzape doc

---

## 08 — Scripts annexos

Còpia de treball dels `.py` / annexos de l’arrel de `archivos\` (alguns es solapen amb 01/02; serveixen com a backup ràpid).

---

## 09 — Gemelos / K3 (no prims, però VMA)

- `gemelo_v9*`, `gemelo_v10_cp_dinamic.py`
- `zypyzape_minigemelo.py`, `zypyzape_twin_v4_8_quijote.py`
- `k3_ultimo_software_hub.tsx` (UI hub, no motor K3 complet)

---

## El que NO s’ha rescatat (a propòsit)

- APKs, instal·ladors (Ollama, LM Studio, WinToUSB…)
- Software portable / cracks / KMS
- `preogramas\` complet
- Desenes de còpies `(1)(2)(3)` de Word
- Dumps de xat HTML/txt de multi-GB dins claude-extract
- Tot OneDrive UPV no relacionat amb prims

El montón original **segueix intacte** si cal anar a buscar quelcom més.

---

## Rutes originals (referència ràpida)

```
monton viejo revisar\archivos\sunraman.docx
monton viejo revisar\archivos\criba_modular.py
monton viejo revisar\archivos\anexoF_criba6kpm1_openmp.c
monton viejo revisar\archivos\claude-main-extract\claude-main\py\
monton viejo revisar\archivos\claude-main-extract\claude-main\PY L5\
monton viejo revisar\archivos\UPV EPSA\...\computacion paralela\primos\
```

---

## Propers passos suggerits

1. Verificar que `cribas.py` i `mdc_v23.py` executen al teu Python.
2. Enllaçar/copiar mòduls canònics cap als skills `vma-methods` i `vma-mdc`.
3. Si cal, segona passada: patents, Wolfram, Libro4 energia, etc.
