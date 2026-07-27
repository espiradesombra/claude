# antipc-port-c

**Port a C** de módulos AntiPC (cribas, MDC trenes, Newton rápido) + inventarios de integración y demos.

Autor: VMA · línea industrial AntiPC dentro del monorepo civil 33×1.

## Qué contiene esta carpeta

| Ruta | Contenido |
|------|-----------|
| `src/*.c` | `criba_hibrida`, `criba_modular6k`, `mdc_trains`, `mdc_ksweep_predict`, `newton_rapido`, `port_demo` |
| `include/antipc_port_v04.h` | Cabecera del port |
| `0x_INTEGRADO_v*.txt` | Logs de integración v04–v14 |
| `00_INVENTARIO_*.md` / `01_MAPA_PORT_A_C.md` | Inventario y mapa del port |
| `06_ANTIPC_NETWORK.txt`, `07_GAME_MMO.txt`, … | Notas de red / MMO / cluster |
| `LEEME.txt`, `bench_port_c.py`, `build_desktop.bat` | Guía, bench, build |
| `k3dedup_demo/` | Mini demo dedup |
| `*.json` (list/create/…) | Plantillas API / tareas |

## Relacionado

- Motor completo Python/DLL: [`../antipc/`](../antipc/)
- Pack mates: [`../VMA_mates_rescat_2026/`](../VMA_mates_rescat_2026/)
- Mapa: [`../MAPA.md`](../MAPA.md)
