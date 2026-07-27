# Demo oficial — Cribas VMA

Pack de cribas con **comparador**, **teoría** y guion de 2 minutos (mismo formato que Sale-It).

Carpeta: [`VMA_mates_rescat_2026/01_cribas/`](../VMA_mates_rescat_2026/01_cribas/)

| Archivo | Uso |
|---------|-----|
| [LEEME_DEMO_1_PAGINA.txt](../VMA_mates_rescat_2026/01_cribas/LEEME_DEMO_1_PAGINA.txt) | Demo en 1 página |
| [GUION_DEMO_2_MINUTOS.txt](../VMA_mates_rescat_2026/01_cribas/GUION_DEMO_2_MINUTOS.txt) | Guion hablado |
| [TEORIA_CRIBAS.md](../VMA_mates_rescat_2026/01_cribas/TEORIA_CRIBAS.md) | Teoría |
| `RUN_COMPARADOR.bat` | Benchmark one-click |

## En 1 paso

```bat
cd VMA_mates_rescat_2026\01_cribas
RUN_COMPARADOR.bat
```

## Guion (tiempos)

| Tiempo | Bloque |
|--------|--------|
| 0:00–0:25 | Enganche: 6k±1 + fase 2p/4p + sunraman |
| 0:25–0:50 | Tres cribas VMA en una frase cada una |
| 0:50–1:25 | Acción: RUN_COMPARADOR (ok + ms) |
| 1:25–1:50 | Números: todas ✓; mejor VMA = Modular |
| 1:50–2:00 | Cierre honesto: Eratóstenes gana en py; C/OpenMP después |

## Resultados de referencia (PC, 2026-07-27, `--quick`)

| L | π(L) | Mejor VMA | Eratóstenes |
|---|------|-----------|-------------|
| 500 | 95 | Modular ~0.04 ms | ~0.03 ms |
| 2 000 | 303 | Modular ~0.14 ms | ~0.10 ms |
| 10 000 | 1229 | Modular ~0.76 ms | ~0.45 ms |

Todas las VMA **correctas** (lista = Eratóstenes).

## Relacionado

- [Era-sunraman-cribas](Era-sunraman-cribas)  
- [Demo-oficial-Sale-It](Demo-oficial-Sale-It)  
- [Matematicas](Matematicas)  
