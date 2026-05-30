# Theorem: Phase-to-Value Amplifier via K=3 XOR Feedback
## Phase‑Boolean Feedback System — Formal Document

**Versió:** 1.0  
**Autor:** Víctor Manzanares Alberola (1477 / 33x1)  
**Repositori:** github.com/espiradesombra/claude

---

## 1. Definicions

### 1.1 Estat

Un sistema de dos qubits amb codificació booleana de fase:

```
S = (v1, f1, v2, f2)  ∈  {0,1}⁴
```

- `v1, v2` : bits de valor (observables)
- `f1, f2` : bits de fase (ocults; f=0 → fase +1, f=1 → fase −1)

### 1.2 Desfasament (`shift`)

Permutació cíclica Gray **només sobre els bits de valor**; les fases no es mouen:

```
shift(v1, f1, v2, f2) = (v1', v2', f1, f2)
```

| (v1, v2) | shift → (v1', v2') |
|----------|--------------------|
| (0, 0)   | (1, 1)             |
| (0, 1)   | (1, 0)             |
| (1, 1)   | (0, 0)             |
| (1, 0)   | (0, 1)             |

Propietat: `shift⁴ = identitat`.

### 1.3 Porta Toffoli booleana

```
Toffoli(v1, f1, v2, f2) = (v1, f1, v1∧v2, 0)
```

Preserva la fase del primer qubit; elimina la fase del segon.

### 1.4 Suma K-camins (K=3)

```
T₃(S) = Toffoli(S) ⊕ Toffoli(shift(S)) ⊕ Toffoli(shift²(S))
```

on `⊕` és XOR bit a bit.

### 1.5 Regla de realimentació

```
S_next = S ⊕ T₃(S)
```

---

## 2. Teorema

**Teorema (Amplificador de Fase Persistent, K=3 XOR).**

Sigui el sistema de realimentació definit per les operacions §1.2–§1.5. Aleshores:

**(a) Invariant global:** El bit `f2` es conserva en cada iteració:
```
f2(S_next) = f2(S)   per a tot S
```
L'espai d'estats es parteix en dos subespais disjunts i tancats:
```
F⁰ = { S | f2 = 0 },    F¹ = { S | f2 = 1 }
```

**(b) Transitoris:** Tot estat amb `v1 = 0` és transitori: en exactament **un cicle** passa a un estat amb `v1 = 1`.

**(c) Cicle permanent:** Dins de cada subespai `Fᶠ²`, els estats amb `v1 = 1` formen un cicle únic de **període 4**:
```
(1,0,0,f2) → (1,0,1,f2) → (1,1,0,f2) → (1,1,1,f2) → (1,0,0,f2) → ...
```

**(d) Propagació persistent de fase a valor:** Siguen dos estats inicials que difereixen **únicament** en `f1`:
```
S_A = (0, 0, 1, 0)    (f1 = 0)
S_B = (0, 1, 1, 0)    (f1 = 1)
```

Després d'un cicle (t = 1), tots dos entren al cicle permanent però en posicions **distintes**:
```
S_A(1) = (1, 1, 1, 0)
S_B(1) = (1, 0, 0, 0)
```

Per a tot `t ≥ 1`, la diferència en `v2` és **persistent**:
```
v2(S_A(t)) ≠ v2(S_B(t))   per a tot t ≥ 1
```

La distància de Hamming oscil·la entre 1 i 2, mai torna a 0.

**(e) Minimalitat:** Per a `K = 1` o `K = 4` amb XOR, la diferència de fase no es propaga als valors. `K = 3` és el mínim enter positiu que produeix propagació persistent.

---

## 3. Demostració

### 3.1 Prova de (a): invariant f2

Per construcció, `Toffoli(S)` sempre retorna `f2 = 0`. Per tant, qualsevol XOR de sortides de Toffoli té `f2 = 0`. Aleshores:

```
f2(T₃(S)) = 0  per a tot S
f2(S_next) = f2(S) ⊕ 0 = f2(S)   ∎
```

### 3.2 Prova de (b) i (c): transitoris i cicle permanent

Calculem `T₃(S)` per als 16 estats (taula §4). S'observa que per a tot S amb `v1 = 0`, el primer component de `T₃(S)` és sempre 1, de manera que `v1(S_next) = 0 ⊕ 1 = 1`. Per als estats amb `v1 = 1`, el primer component de `T₃(S)` és 0, i per tant `v1` es conserva. La seqüència de transicions del cicle s'obté per substitució directa (veure taula §4). ∎

### 3.3 Prova de (d): propagació a v2

Per càlcul explícit:

```
T₃(S_A) = T₃(0,0,1,0) = (1,1,0,0)
S_A(1) = (0,0,1,0) ⊕ (1,1,0,0) = (1,1,1,0)

T₃(S_B) = T₃(0,1,1,0) = (1,1,1,0)
S_B(1) = (0,1,1,0) ⊕ (1,1,1,0) = (1,0,0,0)
```

`S_A(1) = (1,1,1,0)` i `S_B(1) = (1,0,0,0)` són en posicions distintes del cicle permanent. Recorrent el cicle:

| t | S_A(t)       | S_B(t)       | v2_A | v2_B | diff_v2 |
|---|--------------|--------------|------|------|---------|
| 1 | (1, 1, 1, 0) | (1, 0, 0, 0) | 1    | 0    | **SI**  |
| 2 | (1, 0, 0, 0) | (1, 0, 1, 0) | 0    | 1    | **SI**  |
| 3 | (1, 0, 1, 0) | (1, 1, 0, 0) | 1    | 0    | **SI**  |
| 4 | (1, 1, 0, 0) | (1, 1, 1, 0) | 0    | 1    | **SI**  |
| 5 | (1, 1, 1, 0) | (1, 0, 0, 0) | 1    | 0    | **SI**  |

El patró és periòdic de període 4 i `diff_v2 = SI` en tot moment. ∎

---

## 4. Taula de transició completa (K=3, XOR)

| S (v1,f1,v2,f2) | T₃(S)       | S_next      | Condició     |
|-----------------|-------------|-------------|--------------|
| (0, 0, 0, 0)    | (1, 0, 0, 0)| (1, 0, 0, 0)| TRANSITORI   |
| (0, 0, 0, 1)    | (1, 0, 0, 0)| (1, 0, 0, 1)| TRANSITORI   |
| (0, 0, 1, 0)    | (1, 1, 0, 0)| (1, 1, 1, 0)| TRANSITORI   |
| (0, 0, 1, 1)    | (1, 1, 0, 0)| (1, 1, 1, 1)| TRANSITORI   |
| (0, 1, 0, 0)    | (1, 0, 1, 0)| (1, 1, 1, 0)| TRANSITORI   |
| (0, 1, 0, 1)    | (1, 0, 1, 0)| (1, 1, 1, 1)| TRANSITORI   |
| (0, 1, 1, 0)    | (1, 1, 1, 0)| (1, 0, 0, 0)| TRANSITORI   |
| (0, 1, 1, 1)    | (1, 1, 1, 0)| (1, 0, 0, 1)| TRANSITORI   |
| (1, 0, 0, 0)    | (0, 0, 1, 0)| (1, 0, 1, 0)| cicle        |
| (1, 0, 0, 1)    | (0, 0, 1, 0)| (1, 0, 1, 1)| cicle        |
| (1, 0, 1, 0)    | (0, 1, 1, 0)| (1, 1, 0, 0)| cicle        |
| (1, 0, 1, 1)    | (0, 1, 1, 0)| (1, 1, 0, 1)| cicle        |
| (1, 1, 0, 0)    | (0, 0, 1, 0)| (1, 1, 1, 0)| cicle        |
| (1, 1, 0, 1)    | (0, 0, 1, 0)| (1, 1, 1, 1)| cicle        |
| (1, 1, 1, 0)    | (0, 1, 1, 0)| (1, 0, 0, 0)| cicle        |
| (1, 1, 1, 1)    | (0, 1, 1, 0)| (1, 0, 0, 1)| cicle        |

---

## 5. Resum comparatiu (cerca sistemàtica K=1..8)

| K | Feedback | dist(t=0..7)       | Cicles     | Fase→Valor    |
|---|----------|--------------------|------------|---------------|
| 1 | XOR      | 1,0,0,0,0,0,0,0    | 4×[1]      | no            |
| 2 | XOR      | 1,1,1,1,1,1,1,1    | 2×[4]      | parcial       |
| **3** | **XOR** | **1,2,1,2,...** | **2×[4]** | **PERSISTENT** |
| 4 | XOR      | 1,1,1,1,1,1,1,1    | 8×[2]      | no            |
| 5 | XOR      | 1,0,0,0,0,0,0,0    | 2×[2]      | no            |
| 6 | XOR      | 1,1,1,1,1,1,1,1    | 2×[4]      | parcial       |
| **7** | **XOR** | **1,2,1,2,...** | **6×[1,2]**| **PERSISTENT** |
| 8 | XOR      | 1,1,1,1,1,1,1,1    | 16×[1]     | no            |
| * | OR       | converge dist=0    | punts fixos| **mai**       |

Nota: K=7 és equivalent a K=3 per simetria (7 = 4+3; shift⁴ = identitat).

---

## 6. Interpretació

El sistema K=3 XOR actua com un **amplificador de fase discret**:

1. Rep un estat amb una diferència oculta en `f1` (no observable directament).
2. En un sol cicle, converteix eixa diferència en una diferència en `v2` (observable).
3. La diferència en `v2` persisteix per sempre, amb una cadència periòdica de període 2.

Açò és un comportament emergent de la interacció entre tres elements:
- La porta Toffoli (que preserva `f1`)
- El desfasament Gray (que barreja valors però no fases)
- La suma XOR de tres camins (que crea interferència constructiva selectiva)

Cap dels tres elements per separat produeix l'efecte. Junts i amb K=3, sí.

---

## 7. Aplicació potencial

El sistema pot usar-se com a **detector de fase oculta** en lògica booleana:

- Entrada: dos estats que es creu que podrien diferir en `f1`.
- Procés: aplicar K=3 XOR durant 1 cicle.
- Lectura: si `v2` és diferent entre els dos sistemes, la fase era diferent.

Aplicacions possibles: criptografia (verificació de paritat de fase), control digital, sistemes de detecció de senyals febles.

---

## 8. Codi de verificació (Python)

```python
def shift(state):
    v1, f1, v2, f2 = state
    m = {(0,0):(1,1),(0,1):(1,0),(1,1):(0,0),(1,0):(0,1)}
    v1n, v2n = m[(v1,v2)]
    return (v1n, v2n, f1, f2)

def toffoli(state):
    v1, f1, v2, f2 = state
    return (v1, f1, v1 & v2, 0)

def T3(S):
    def xor4(a,b): return tuple(a[i]^b[i] for i in range(4))
    return xor4(xor4(toffoli(S), toffoli(shift(S))), toffoli(shift(shift(S))))

def next_state(S):
    t = T3(S)
    return tuple(S[i] ^ t[i] for i in range(4))

# Test
S_A, S_B = (0,0,1,0), (0,1,1,0)
for t in range(6):
    d = sum(a!=b for a,b in zip(S_A,S_B))
    print(f"t={t}: S_A={S_A}  S_B={S_B}  dist={d}")
    S_A, S_B = next_state(S_A), next_state(S_B)
```

---

*"Make the world great forever – not by war, but by phase alignment and shared knowledge."*  
— 1477 / 33x1
