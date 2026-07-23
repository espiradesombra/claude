# VMA mates rescat + XFI (2026-07-23)

Pack net rescatat de `monton viejo revisar` + fixes + gemelo XFI.

**Autor material:** Víctor Manzanares Alberola (VMA)  
**Sessió tècnica:** Grok (xAI)  
**Marc:** [33×1](../33x1/) — ús civil (`../33x1/02_USO_CIVIL.txt`)

## Contingut

| Carpeta | Contingut |
|---------|-----------|
| `00_sunraman_eratostenes/` | Pseudocodi sunraman + docs Eratòstenes |
| `01_cribas/` | Desmemoriada, modular 6k±1 (**fase 2p/4p fix**), híbrida, OpenMP C |
| `02_criva_mrauv/` | Criva densitat + MRAUV |
| `03_mdc/` | MDC + `mdc_v18`…`v23` |
| `04_newton_rapido/` | Newton Rápido + papers |
| `05_sofi_fermat_goldbach/` | Sofí, Fermat, Goldbach, **siguiente_primo fix** |
| `06_upv_c_historics/` | Experiments C UPV |
| `07_documents_canonics/` | Tractats VMA |
| `08_scripts_anexos_VMA/` | Scripts curts |
| `09_gemelos_zypyzape_ultims/` | Gemelos Quijote/ZZ |
| `10_inicio_quijote_zypyzape/` | Context + **XFI N=3/4** |

## Proves ràpides

```bat
cd 01_cribas
python -c "from cribas import comparar_cribas; comparar_cribas(500)"

cd ..\05_sofi_fermat_goldbach
python -c "from siguiente_primo import validate; validate(40)"

cd ..\10_inicio_quijote_zypyzape
python gemelo_xfi_avion_4turbinas.py --n 3
python gemelo_xfi_avion_4turbinas.py --compare
```

## Declaració mundial 33×1

Veure: [`../33x1/05_DECLARACIO_33x1_GROK_WORLD.txt`](../33x1/05_DECLARACIO_33x1_GROK_WORLD.txt)

## Integritat

Regenerar hashes del paquet 33×1:

```bat
cd ..\33x1
python generar_manifiesto.py
```
