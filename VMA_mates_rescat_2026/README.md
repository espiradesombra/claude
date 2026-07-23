# VMA_mates_rescat_2026

**Pack neto de matemáticas y gemelos rescatados** (2026-07-23) del montón del PC + fixes + XFI.

Forma parte del **“1”** del trato **33×1**  
(ver [`../33x1/00_QUE_ES_33x1.md`](../33x1/00_QUE_ES_33x1.md)).

**Autor material:** Víctor Manzanares Alberola (VMA)  
**Sessió tècnica:** Grok (xAI)

## Subcarpetas (cada una con su README)

| Carpeta | Descripción |
|---------|-------------|
| [`00_sunraman_eratostenes`](00_sunraman_eratostenes/) | Pseudocódigo sunraman + docs Eratóstenes |
| [`01_cribas`](01_cribas/) | Cribas VMA (desmemoriada, modular **fix**, híbrida, OpenMP) |
| [`02_criva_mrauv`](02_criva_mrauv/) | Criva densidad + MRAUV + saltos |
| [`03_mdc`](03_mdc/) | MDC + evolución mdc_v18…v23 |
| [`04_newton_rapido`](04_newton_rapido/) | Newton Rápido + oráculos + papers |
| [`05_sofi_fermat_goldbach`](05_sofi_fermat_goldbach/) | Sofí, Fermat, Goldbach, siguiente_primo **fix** |
| [`06_upv_c_historics`](06_upv_c_historics/) | Código C histórico UPV (cribas paralelas) |
| [`07_documents_canonics`](07_documents_canonics/) | Word/PDF canónicos VMA primos |
| [`08_scripts_anexos_VMA`](08_scripts_anexos_VMA/) | Scripts cortos de anexos |
| [`09_gemelos_zypyzape_ultims`](09_gemelos_zypyzape_ultims/) | Gemelos Quijote/ZZ v9–v10 |
| [`10_inicio_quijote_zypyzape`](10_inicio_quijote_zypyzape/) | Contextos + **XFI avión N=3/4** |

## Verificación

```bat
cd 01_cribas
python -c "from cribas import comparar_cribas; comparar_cribas(500)"

cd ..\05_sofi_fermat_goldbach
python -c "from siguiente_primo import validate; validate(40)"

cd ..\10_inicio_quijote_zypyzape
python gemelo_xfi_avion_4turbinas.py --n 3
```

- `INVENTARI.md` — inventario del rescat
- `VERIFICACIO_I_INICIS.md` — tests + inicios Quijote/Zypy
- Uso civil: `../33x1/02_USO_CIVIL.txt`
