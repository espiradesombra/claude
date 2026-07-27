# Demo oficial — Sale-It (AntiPC)

**Esta es la demo prioritaria del monorepo** cuando hay que enseñar algo usable en minutos.

Carpeta: [`sale-it/`](../sale-it/)  
Hoja de 1 página: [`sale-it/LEEME_DEMO_1_PAGINA.txt`](../sale-it/LEEME_DEMO_1_PAGINA.txt)  
Números de este PC: [`sale-it/RESULTADOS_DEMO_PC_2026-07-27.txt`](../sale-it/RESULTADOS_DEMO_PC_2026-07-27.txt)

## En 3 pasos

```bat
cd sale-it\scripts
01_instalar.bat
05_verificar_k3.bat
03_benchmark_udp.bat
```

## Resultados de referencia (PC VMA, 2026-07-27)

| Prueba | Resultado |
|--------|-----------|
| K3 `verify_k3.py` | 5/5 digests OK |
| A–E UDP (1000 pkt, 4 hubs) | **E gana**: 1591 pkt/s vs A 843; ALU E = 53% de A (~**47% menos**) |
| demo_v010 | 26.6% ALU ahorrada |

## Guion de 2 minutos (voz)

Texto listo para leer o memorizar por bloques:

**[`sale-it/GUION_DEMO_2_MINUTOS.txt`](../sale-it/GUION_DEMO_2_MINUTOS.txt)**

| Tiempo | Bloque |
|--------|--------|
| 0:00–0:25 | Enganche: red + `P_util = N·E + K` |
| 0:25–0:45 | Qué hay en el pack |
| 0:45–1:20 | Acción: K3 + `03_benchmark_udp.bat` |
| 1:20–1:50 | Números: E ~1590 pkt/s, ALU ~53% de A |
| 1:50–2:00 | Cierre honesto + siguiente nivel |

## Qué decir (resumen)

- “Corre el bat 03 y mira la fila **E** (UDP real).”  
- Prototipo **Python** reproducible en localhost; no es DPDK/C++ de producción.  
- Fórmula: `P_util(N) = N · E(N) + K(N)`.

## Relacionado

- Motor completo: `antipc/`  
- Ficha de producto: `sale-it/FICHA_PRODUCTO.txt`  
- Wiki: [AntiPC](AntiPC)
