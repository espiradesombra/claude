# Mapa: qué pasar a C y qué dejar en Python

Prioridad por **ganancia CPU** × **facilidad** × **uso en AntiPC**.

---

## Prioridad ALTA — siguiente oleada (prototipos en esta carpeta)

### 1. Criba Modular 6k±1 real (`criba_modular6k.c`)
- **Fuente:** `repo\vma-methods\vma_methods\cribas.py` → `CribaModular6k.run()`
- **Por qué C:** mismo patrón que Eratóstenes ya en C; `bytearray` + saltos por primo = memoria compacta, cache-friendly.
- **Ganancia estimada:** 3–10× vs Python para `limit > 10⁵`.
- **Nota:** la versión VMA con saltos +2p/+4p y `anteriorNY` es fase 1b; la actual en Python ya usa Eratóstenes mod 6 (rápida).

### 2. Criba Desmemoriada (`criba_desmemoriada.c`)
- **Fuente:** `cribas.py` → `CribaDesmemoriada`
- **Por qué C:** bucles anidados sobre patrones booleanos; muchas alloc en Python.
- **Ganancia:** alta para `limit > 5000`; complejidad de implementación media.

### 3. Criba Híbrida ascendente/descendente (`criba_hibrida.c`)
- **Fuente:** `cribas.py` → `CribaHibrida.run()`
- **Por qué C:** dos pasadas + residuos; ideal en `unsigned char[]`.
- **Ganancia:** media-alta; útil en benchmarks vma-methods.

### 4. MDC dos trenes — scan entero (`mdc_trains.c`)
- **Fuente:** `repo\antipc\src\antipc\mdc_lib\mdc_analyze.py`
- **Por qué C:** hoy usa `Fraction` (GCD grande en Python); el núcleo es aritmética entera: `(n-6x-9) % (4x+6) == 0`.
- **Ganancia:** 10–50× en `mdc analyze` para n de 6–10 dígitos.
- **Dejar en Python:** `format_report`, GIF, Tkinter.

### 5. K-sweep predictivo entero (`mdc_ksweep_predict.c`)
- **Fuente:** `repo\filestot l5\ksweep_predictiu.py`
- **Por qué C:** diseñado para enteros puros; evita pérdida float64 en N grandes.
- **Ganancia:** crítica si se escala MDC más allá del toy actual.

### 6. Newton Rápido (`newton_rapido.c`)
- **Fuente:** `repo\vma-methods\vma_methods\newton.py`
- **Por qué C:** ~200 iter × `pow`/`log`; trivial en `double`.
- **Ganancia:** moderada (ya es rápido en Python); útil para GUI sin GIL.

### 7. π ofuscado + e convergente (`k3_geo_math.c`)
- **Fuente:** `encriptacionGeometrica\deepseek_c_20260716_f0b9a1.c`, `gemini-code-1784158392232.py`
- **Por qué C:** Taylor `long double` + backtracking; sustituye `Decimal` lento.
- **Ganancia:** alta para Aleatorovix / cifrado por fase.

### 8. Generador Aleatorovix (`aleatorovix.c`)
- **Fuente:** `gemini-code-1784158392232.py` → clase `Aleatorovix`
- **Por qué C:** `sin` en cada bit; millones de bits en chorro XOR.
- **Ganancia:** muy alta en cifrado masivo.

---

## Prioridad MEDIA

| Módulo | Fuente | C | Motivo quedarse en Python |
|--------|--------|---|---------------------------|
| Criva fractal + Mertens | `criva.py` | Sí (double) | Tablas/compare siguen en Python |
| `pi_count` / trial division | `classic.py` | Sí | Ya cubierto por sieve C |
| Grafcet `_heavy_hash` | `grafcet.py` | Parcial | Bucles K3 → llamar `k3_hash_buffer` C |
| `redundant_hash` × N pkt | `grafcet.py` | Parcial | Orquestación Grafcet en Python |
| `analisi_modular` | `mdc_lib/analisi_modular.py` | Sí | Fase/patron_bits son baratos |
| `regla_calculo` | `mdc_lib/regla_calculo.py` | Sí | `exp`/`log` escala Battery |
| UDP slot ring | `network/bd_pipeline.py` | No* | I/O bound; *C solo si zero-copy extremo |
| k3hash dedup/search | `k3hash/examples/` | Sí | Ya hay C en HASHTOOLCODE; exponer en DLL |

---

## Prioridad BAJA — mantener Python

| Módulo | Motivo |
|--------|--------|
| GUI vma-methods (`gui/app.py`) | Tkinter |
| MDC visual / GIF / matplotlib | Render |
| Microkernel, plugins, ledger | Arquitectura, no ALU |
| WaveProvider (RTT ping) | Red / OS |
| Gemelos `gemelo_v*.py` | SciPy/NumPy, plots, pocos segundos de sim |
| `33x1` manifiesto, docs | Sin cómputo pesado |
| CLI argparse, JSON export | Glue |

---

## Matriz resumen

```
                    GANANCIA CPU
                 baja    media    alta
              ┌────────┬────────┬────────┐
        fácil │ criva+ │ newton │ eratos │  ← HECHO
              │ mertens│ modbits│ geo    │  ← HECHO
              ├────────┼────────┼────────┤
     media    │ regla  │ hibrida│ 6k±1   │  ← PROTOTIPOS AQUÍ
              │ calc   │ desmem │ mdc_tr │
              ├────────┼────────┼────────┤
        difíc │ UDP    │ gemelo │ ksweep │
              │        │ backtr │ aleat. │
              └────────┴────────┴────────┘
```

---

## Orden sugerido v0.4 → v0.5

1. `mdc_trains.c` + `mdc_ksweep_predict.c` (AntiPC `mdc analyze` rápido)
2. `criba_modular6k.c` + `criba_hibrida.c` (vma-methods bench)
3. `newton_rapido.c` (vma-methods tab Newton)
4. `k3_geo_math.c` + `aleatorovix.c` (encriptacionGeometrica industrial)
5. Exponer k3dedup/k3search desde HASHTOOLCODE en la misma DLL

---

## Cómo integrar sin mover archivos

1. Prototipar y bench en `Desktop\antipc-port-c\` (esta carpeta).
2. Cuando un `.c` pase tests, **copiar** (no mover) a `repo\antipc\src\antipc\native\antipc_core\src\`.
3. Añadir símbolo a `antipc_native.h` y `native_engine.py`.
4. Recompilar con `repo\antipc\scripts\21_build_antipc_native.bat`.