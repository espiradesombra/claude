# AntiPC

**Adaptive Network Through Parallel Computing** — red que no solo transporta datos: cada paquete puede ser unidad de cómputo reutilizable.

> *Transforming communication into computation.*

## Fórmula operativa

`P_util(N) = N · E(N) + K(N)`

| Símbolo | Idea |
|---------|------|
| N | Nodos / hubs |
| E(N) | Eficiencia de coordinación |
| K(N) | Conocimiento reutilizable (cache / 2ª pasada) |

## Carpetas

| Carpeta | Contenido |
|---------|-----------|
| `antipc/` | Motor principal: CLI, UDP, MMO, DLL ALU, docs |
| `antipc-port-c/` | Port C (cribas, MDC trenes, Newton) + logs v04–v14 |
| `antipc2/`, `cpu-antipc/` | Evoluciones / runtime CPU |
| `sale-it/` | Ficha producto, licencia, demos go-to-market |
| `ideas-para-gpt-antipc/` | Roadmaps históricos y pack ejecutable grande |

## Entrada práctica (recomendado)

**Demo de producto:** carpeta [`sale-it/`](../sale-it/) — ver wiki [Demo-oficial-Sale-It](Demo-oficial-Sale-It).

```bat
cd sale-it\scripts
03_benchmark_udp.bat
```

Motor completo (más grande):

```bat
cd antipc
REM LEEME.txt · INICIO.bat · antipc.cmd
```

README industrial largo: `antipc/README.md`.
