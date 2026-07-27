# Teoría — cribas VMA y el comparador

**Autor material:** Víctor Manzanares Alberola (VMA)  
**Pack:** `VMA_mates_rescat_2026/01_cribas`  
**Código:** `cribas.py` · `benchmark_cribas.py`  
**Origen sunraman:** `../00_sunraman_eratostenes/`

---

## 1. Qué compara el benchmark

| id | Método | Familia |
|----|--------|---------|
| `desmemoriada` | Criba Desmemoriada | VMA |
| `modular6k` | Criba Modular 6k±1 (fase 2p/4p) | VMA |
| `hibrida` | Híbrida ascendente + descendente | VMA |
| `hibrida_seg` | Híbrida segmentada | VMA |
| `eratostenes` | Eratóstenes clásico | clásica |
| `rueda6k` | Eratóstenes solo en 6k±1 | clásica |
| `trial` | Trial division | referencia (solo L pequeños) |

**Corrección:** se valida cada lista de primos contra Eratóstenes (`π(L)` y conjunto exacto).

**Tiempo:** mejor de N corridas (tras warm-up), en milisegundos.

---

## 2. Eratóstenes clásico (referencia)

- Array booleano de tamaño `L+1`.
- Para cada primo `p ≤ √L`, marcar múltiplos desde `p²` con paso `p`.
- Complejidad práctica ≈ **O(L log log L)** marcas.
- Memoria Θ(L).

Sirve de **oro** de corrección y de reloj de referencia en el CSV (`speedup_vs_erat`).

---

## 3. Rueda 6k±1

Todo primo `> 3` es **6k−1** o **6k+1**.  
Se puede:

1. No iterar pares ni múltiplos de 3 como candidatos a primo.
2. Al marcar múltiplos de un primo `p ≡ ±1 (mod 6)`, usar solo saltos que caen otra vez en 6k±1.

### Fase de saltos 2p / 4p (crítica)

Para un primo `p > 3`:

- `p ≡ 5 (mod 6)` ⇒ `p² ≡ 1 (mod 6)` y el siguiente múltiplo 6k±1 de `p` está a **+2p**.
- `p ≡ 1 (mod 6)` ⇒ el siguiente adecuado está a **+4p**.

Luego se **alterna** 2p y 4p.

Si la fase está mal (p. ej. empezar siempre en +2p), se **saltan compuestos** (ejemplo histórico: `p=7` y el 77).  
**Fix documentado 2026-07-23** en `CribaModular6k` del pack.

---

## 4. Sunraman → modular (ancestro)

Pseudocodi en `../00_sunraman_eratostenes/sunraman_pseudocodi.txt` y Word `sunraman.docx`.

Ideas:

| Idea | Expresión |
|------|-----------|
| Forma de “primos” en la rejilla | `2n+3` (cubre 6k±1 con reindexación) |
| Forma de compuestos | `(2n+3)(2m+3) = 4mn+6m+6n+9` |
| Anti-remarcado | índices `anterior` / `anteriorM` / `anteriorNY` |

**Línea evolutiva:**

```
sunraman  →  criba_modular / Modular6k  →  anexoF OpenMP C  →  cribas.py
```

---

## 5. Criba Desmemoriada (VMA)

Documento y guion: `criba_desmemoriada.txt`.

Idea: la primalidad se expresa con **patrones booleanos** de periodo ligado a `6·p` (estructura mod 6 × primo).  
En lugar de re-marcar a ciegas todo el array, se:

1. Construye un patrón compacto de longitud `6p`.
2. Se aplica de forma cíclica / por AND lógico con la lista de candidatos.

**Objetivo teórico:** menos reescrituras de memoria y menos remarcados (“desmemoria” de marcas ya hechas).

La implementación en `cribas.py` es **didáctica** (puede no ser la más rápida en Python puro frente a Eratóstenes vectorizado en C); el comparador mide la realidad en tu máquina.

---

## 6. Criba Modular 6k±1 (VMA)

- Solo trabaja la semántica de candidatos 6k±1 (+ 2 y 3).
- Marcado con saltos **2p/4p** y fase correcta.
- En el marco teórico largo (cotas / estructuras) aparecen cotas del tipo fracción de `l²` al marcar una sola vez (`anteriorNY`).

Complejidad de marcado: del orden de recorrer múltiplos en la subrejilla 6k±1 (menor densidad que todos los enteros).

---

## 7. Criba Híbrida (VMA)

1. **Ascenso** hasta `mid = L/2` (o `√L` en variantes).
2. **Residuos** `L mod p` para cada primo pequeño.
3. **Descenso** desde arriba marcando múltiplos hacia `mid`.
4. Cierre de la zona media con una pasada ascendente de seguridad.

**Motivación:** en cribas clásicas, muchos “fallos”/remarcados se concentran al final del intervalo; el híbrido reparte trabajo entre mitad inferior y superior.

**Segmentada:** trocea `[2, L]` en bloques; en cada bloque se usan primos base hasta `√L` (estilo segmentado clásico + ideas híbridas).

---

## 8. Cómo leer el benchmark

| Columna | Significado |
|---------|-------------|
| `#π` | Número de primos devueltos |
| `ms` | Mejor tiempo en milisegundos |
| `vs Erat` | `t_erat / t_metodo` (>1 = más rápido que Eratóstenes) |
| `ok` | Igualdad exacta de la lista de primos con Eratóstenes |

### Interpretación honesta

- En **Python puro**, Eratóstenes bien escrito suele ganar a implementaciones didácticas VMA en tiempo wall-clock.
- El valor de las cribas VMA aquí es **estructura** (6k±1, fase, desmemoria, híbrido) y **correctez** del modelo modular — base para ports C/OpenMP (`anexoF_criba6kpm1_openmp.c`).
- Para reivindicar velocidad industrial: compilar el anexo C / usar `antipc-port-c`, no solo este script.

---

## 9. Comandos

```bat
cd VMA_mates_rescat_2026\01_cribas

REM rápido
python benchmark_cribas.py --quick

REM límites a medida
python benchmark_cribas.py --limits 1000,10000,100000 --repeats 5 --no-trial

REM solo métodos VMA + Eratóstenes
python benchmark_cribas.py --methods desmemoriada,modular6k,hibrida,hibrida_seg,eratostenes

REM listar ids
python benchmark_cribas.py --list

REM API corta heredada
python -c "from cribas import comparar_cribas; comparar_cribas(5000)"
```

One-click: `RUN_COMPARADOR.bat`

---

## 10. Referencias en el monorepo

| Ruta | Contenido |
|------|-----------|
| `00_sunraman_eratostenes/` | Word + pseudocodi sunraman |
| `01_cribas/cribas.py` | Tres cribas + `comparar_cribas` |
| `01_cribas/Cribas_cotas_y_estructuras_de_primos_VMA.txt` | Texto largo cotas |
| `01_cribas/anexoF_criba6kpm1_openmp.c` | C paralelo |
| `teoremas/cribas/` | Fichas formales |
| `wiki/Era-sunraman-cribas.md` | Mapa de copias y canónicos |
| `vma-methods/` | CLI de librería |

**33×1:** parte del **1** civil del monorepo. Uso civil / educativo.
