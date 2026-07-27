# Era sunraman / cribas (mapa de archivos)

Rescate y ubicación de material **viejo** de la época sunraman / cribas 6k±1 / Eratóstenes.  
**Conclusión:** lo canónico y usable **ya está en el monorepo**; fuera hay sobre todo copias.

## Cadena evolutiva (cómo se relacionan)

```
sunraman (pseudocodi / Word)
    → criba_modular.py / modular 6k±1
    → anexoF_criba6kpm1_openmp.c  (C + OpenMP)
    → cribas.py  (Desmemoriada + Modular + Híbrida)
    → vma-methods / antipc-port-c / aleatorovix.nucleo.criba
```

Fórmulas sunraman (INVENTARI del pack rescat):

- primos ≈ `2n+3`
- compuestos ≈ `4mn + 6m + 6n + 9`
- marcado con `anterior` / `anteriorM` / `anteriorNY`

## Dónde está lo canónico (usar esto)

| Prioridad | Ruta | Qué hay |
|:---------:|------|---------|
| **1** | `VMA_mates_rescat_2026/00_sunraman_eratostenes/` | `sunraman.docx`, `sunraman_pseudocodi.txt`, docs Eratóstenes |
| **2** | `VMA_mates_rescat_2026/01_cribas/` | `cribas.py`, modular, híbrida, desmemoriada.txt, OpenMP C, cotas VMA |
| **3** | `VMA_mates_rescat_2026/06_upv_c_historics/` | C histórico UPV (p. ej. desmemoriado) |
| **4** | `teoremas/cribas/` | Fichas CRIBA_* (incl. `CRIBA_SUNRAMAN_ERATOSTENES.txt`) |
| **5** | `vma-methods/` | Librería ejecutable (cribas por CLI) |
| **6** | `antipc-port-c/src/criba_*.c` | Port C integrado |
| **7** | `aleatorovix/nucleo/criba.py` | Criba en el organismo |

Inventario detallado: [`VMA_mates_rescat_2026/INVENTARI.md`](../VMA_mates_rescat_2026/INVENTARI.md)

## Prueba rápida (código de esa época, limpio)

```bat
cd VMA_mates_rescat_2026\01_cribas
python -c "from cribas import comparar_cribas; comparar_cribas(500)"
```

Pseudocodi legible:

```bat
type VMA_mates_rescat_2026\00_sunraman_eratostenes\sunraman_pseudocodi.txt
```

## Copias / variantes “de la época” (ya en monorepo, no re-subir)

| Sitio | Notas |
|-------|--------|
| `Nueva carpeta0/`, `Nueva carpeta (3)/` | Word **Suntaman** / cribratge segmentat / comparativa de cribas (ortografía “Suntaman” = misma línea) |
| `coses velles/sunraman.docx`, `diamante/sunraman*.docx` | Copias del Word sunraman |
| `junto/`, `nuevo/`, `Libro3…` | `Cribas_cotas_y_estructuras_de_primos_VMA*.docx` |
| `archivos-vma/`, `just run/archivos-vma/` | guion desmemoriada + cotas txt + código anexo |
| `filesclaude 6-5/cribas.py`, `py/cribas.py` | Versiones de `cribas.py` |
| Escritorio `VMA_mates_rescat_2026` | **Misma** pack 00/01 que en GitHub |

## ¿Falta algo fuera sin subir?

Auditoría rápida: **no** aparece un “sunraman perdido” grande fuera del monorepo.  
Lo del Escritorio / `repo/` es **duplicado** o ya empaquetado en `VMA_mates_rescat_2026/`.

## Criterio

No re-subir más copias de `sunraman.docx`.  
Si hay que enseñar esta época: **pack rescat 00+01** + opcional ficha `teoremas/cribas/`.
