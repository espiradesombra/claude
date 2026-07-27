# sale-it — go-to-market AntiPC / VMA

Material de **venta y producto**: ficha, licencia, scripts de demo y runtime de referencia.

## Demo oficial (empieza aquí)

| Archivo | Uso |
|---------|-----|
| **[LEEME_DEMO_1_PAGINA.txt](LEEME_DEMO_1_PAGINA.txt)** | Demo en 1 página + mensaje a terceros |
| **[RESULTADOS_DEMO_PC_2026-07-27.txt](RESULTADOS_DEMO_PC_2026-07-27.txt)** | Números medidos en este PC |
| `scripts/03_benchmark_udp.bat` | Benchmark A–E (ganador típico: E UDP) |
| `scripts/05_verificar_k3.bat` | Paridad K3 |

```bat
cd sale-it\scripts
01_instalar.bat
05_verificar_k3.bat
03_benchmark_udp.bat
```

## Qué contiene esta carpeta

| Pieza | Contenido |
|-------|-----------|
| `FICHA_PRODUCTO.txt`, `PRODUCT_SHEET_EN.txt`, `LICENCIA.txt` | Producto y licencia |
| `INICIO.bat`, `LEEME.txt`, `MANIFEST.txt` | Arranque e inventario |
| `src/`, `scripts/`, `docs/`, `referencia/` | Código, scripts, docs y k3hash |
| `LEEME_DEMO_1_PAGINA.txt` | Hoja de demo oficial |
| `RESULTADOS_DEMO_PC_*.txt` | Informes de corrida |

## Enlaces

- Mapa monorepo: [`../MAPA.md`](../MAPA.md)
- 33×1: [`../33x1/`](../33x1/)
- Motor completo: [`../antipc/`](../antipc/)
