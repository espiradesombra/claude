# 01_cribas — tres cribas VMA + comparador

**Pack mates rescat** · Víctor Manzanares Alberola · 33×1 (uso civil)

## Qué hay

| Archivo | Contenido |
|---------|-----------|
| `cribas.py` | Desmemoriada, Modular 6k±1 (**fase 2p/4p fix**), Híbrida + segmentada |
| **`benchmark_cribas.py`** | **Comparador**: VMA + Eratóstenes + rueda 6k±1 + trial; CSV/TXT |
| **`TEORIA_CRIBAS.md`** | Teoría sunraman → modular, fase, híbrida, cómo leer el bench |
| `RUN_COMPARADOR.bat` | One-click (modo `--quick`) |
| `criba_modular.py` / `criba_hibrida.py` | Scripts cortos anexos |
| `anexoF_criba6kpm1_openmp.c` | C + OpenMP |
| `criba_desmemoriada.txt` | Guion desmemoriada |
| `Cribas_cotas_y_estructuras_de_primos_VMA.txt` | Texto técnico cotas |
| `results/` | Salidas CSV/TXT del benchmark (generadas) |

Ancestro sunraman: [`../00_sunraman_eratostenes/`](../00_sunraman_eratostenes/)

## Comparador (recomendado)

```bat
RUN_COMPARADOR.bat
```

```bat
python benchmark_cribas.py --quick
python benchmark_cribas.py --limits 1000,10000,100000 --repeats 5 --no-trial
python benchmark_cribas.py --list
```

API corta:

```bat
python -c "from cribas import comparar_cribas; comparar_cribas(5000)"
```

## Validación

Cada método se contrasta con **Eratóstenes** (lista exacta de primos).  
Si sale `✗`, hay bug de fase o de límites — no ignorar.

## Enlaces monorepo

- Wiki era: [`../../wiki/Era-sunraman-cribas.md`](../../wiki/Era-sunraman-cribas.md)
- Inventario pack: [`../INVENTARI.md`](../INVENTARI.md)
- CLI librería: [`../../vma-methods/`](../../vma-methods/)
