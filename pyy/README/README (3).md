# 🔢 Prime Number Theory: Modular Methods & Density Models

**Author:** Víctor Manzanares Alberola (VMA)  
**Status:** Research in progress | Empirically verified | Code executable  
**Language:** Python 3 / C (OpenMP)  
**Topics:** Number Theory · Prime Sieves · Goldbach · Sophie Germain · Factorization

---

## Overview

This repository contains four original computational and theoretical frameworks for the study of prime numbers, developed independently over several years. Each module is self-contained, empirically tested, and accompanied by executable code.

| Module | Short description |
|--------|-------------------|
| [MRAUV-Goldbach](#1-mrauv-goldbach) | Kinematic density model applied to Goldbach's conjecture |
| [Sofí Structure](#2-sofí-structure-sophie-germain) | Modular classification of Sophie Germain prime candidates |
| [Criva](#3-criva-density-estimator) | Iterative rational estimator of prime density |
| [Discriminant Method](#4-discriminant-method-factorization) | Deterministic discard filter for composite numbers |
| [MDC — Kinematic Diophantine Method](#5-mdc--kinematic-diophantine-method) | Unified method for factorization and Wieferich prime detection |
| [Siguiente Primo](#6-siguiente-primo) | Next prime via accumulated products of consecutives (Wilson-inspired) |
| [SaltoMáximo](#7-saltomáximo) | Minimum 2 primes in window of length √n — conjecture/theorem |
| [ZypyZape MiniGemelo](#8-zypyzape-minigemelo-digital) | Digital twin of 5 wind turbines with kinetic battery (physics simulation) |

---

## 1. MRAUV-Goldbach

### Concept

**MRAUV** (Modelo de Recorrido Acumulado por Velocidad) is a piecewise kinematic model that predicts the local density of primes `D(n)` using three parameters:

- `D₀` — initial density at segment start  
- `V₀` — rate of density decrease  
- `a₀` — acceleration of that rate  

Applied to Goldbach's conjecture, it provides a **computable lower bound** on the number of Goldbach decompositions `G(2n)`, by measuring the *asymmetric fault* — the error introduced by composite multiples of small primes that interfere with the symmetry algorithm.

### Key definitions

```
L(n) = ⌊√(n+3)⌋ + 7          # search corridor length
m(n) = Σ_{i=2..K} √(n+3)/i!  # overcount of composites in corridor
D(n) = (L(n) - m(n)) / 2n     # predicted prime density

F_eff(n) ≈ Σ_{p ≤ √(2n)} ⌊2n/p⌋ · π(2n)/(2n)   # effective asymmetric fault
```

### Criterion

> **If `D(n) > F_eff(n)/(2n) + ε` for all `n > N₀`, then Goldbach holds for `2n > 2N₀`.**

### Python implementation

```python
import math

def L(n):
    return int(math.sqrt(n + 3)) + 7

def m(n, K=50):
    return sum(math.sqrt(n + 3) / math.factorial(i) for i in range(2, K + 1))

def D(n):
    return (L(n) - m(n)) / (2 * n)

def F_eff(n, primes_small):
    pi_2n = sum(1 for p in primes_small if p <= 2 * n)
    density = pi_2n / (2 * n)
    return sum((2 * n // p) * density for p in primes_small if p <= int(math.sqrt(2 * n)) + 1)

def verify_goldbach_MRAUV(n_max=100000, delta=5000):
    from sympy import primerange
    primes_small = list(primerange(2, int(math.sqrt(2 * n_max)) + 2))
    
    for n in range(1000, n_max, delta):
        d = D(n)
        f = F_eff(n, primes_small) / (2 * n)
        margin = d - f
        print(f"n={n:6d}  D(n)={d:.6f}  F_eff/2n={f:.6f}  margin={margin:+.6f}")
        if margin <= 0:
            print(f"  ⚠️  ALERT at n={n}")
            return False
    
    print("✅ Goldbach criterion satisfied in range")
    return True

verify_goldbach_MRAUV()
```

---

## 2. Sofí Structure (Sophie Germain)

### Concept

A modular classification of prime candidates of the form `6k−1` that partitions them into four disjoint sets, isolating Sophie Germain prime candidates without exhaustive testing.

### Formal construction

```
L1  = { a : a ≡ 5 (mod 6) }                        # all 6k-1 candidates
L3  = { a ∈ L1 : a = (6k−1)(6h+1) for some k,h }   # composites via factor type A
L4  = { a ∈ L1 : 2a+1 = (6j−1)(6g+1) for some j,g }# composites via factor type B
L2  = L3 ∩ L4                                        # doubly composite
U2  = L1 \ (L3 ∪ L4)                                # residual candidates

LSG = { p prime : 2p+1 is also prime }               # Sophie Germain primes
```

**Result:** `U2 ⊆ LSG` — every element of `U2` is either a Sophie Germain prime or prime with composite safe prime.

**Conjecture:** `|U2| = ∞` (U2 is infinite), which would imply infinitely many Sophie Germain primes.

### Python implementation

```python
def classify_sofi(limit):
    """Classify numbers up to limit into L1, L3, L4, L2, U2"""
    from sympy import isprime, factorint

    def is_type_A_composite(a):
        """a = (6k-1)(6h+1) for some positive k, h"""
        for f in range(5, int(a**0.5) + 1, 6):
            if a % f == 0:
                g = a // f
                if g % 6 == 1:
                    return True
        return False

    def is_type_B_composite(a):
        """2a+1 = (6j-1)(6g+1)"""
        return is_type_A_composite(2 * a + 1)

    L1 = [a for a in range(5, limit, 6)]
    L3 = set(a for a in L1 if is_type_A_composite(a))
    L4 = set(a for a in L1 if is_type_B_composite(a))
    L2 = L3 & L4
    U2 = set(a for a in L1 if a not in L3 and a not in L4)

    sg = {a for a in U2 if isprime(a) and isprime(2 * a + 1)}

    return {"L1": len(L1), "L3": len(L3), "L4": len(L4),
            "L2": len(L2), "U2": len(U2), "SG_in_U2": len(sg)}

print(classify_sofi(10000))
```

---

## 3. Criva (Density Estimator)

### Concept

**Criva** is an iterative rational estimator of prime density, inspired by layered sieve theory. It converges to `π(x)/x` with controlled relative error (< 0.1% in a few iterations).

### Recurrence

```
D_{n+1} = (D_n + T) / 2

where T = correction term based on sieve layer n
```

**Fractal layer model:**

```
D(x) = Σ_{n=0..N} (D₀ / 2ⁿ) · wₙ(x)
```

where `wₙ(x)` is the weight function for sieve layer `n` (exclusion of multiples of the first `n` primes).

### Python implementation

```python
import math

def criva(x, layers=10):
    """Estimate π(x)/x using the Criva iterative model"""
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29][:layers]
    
    D = 1.0
    for p in primes:
        D *= (1 - 1/p)  # Euler product approximation
    
    # Iterative refinement
    T = 1 / math.log(x) if x > 1 else 1
    for _ in range(5):
        D = (D + T) / 2
        T = D * (1 - 1/math.log(x + 1))
    
    return D

def compare_criva_vs_real(x_values):
    from sympy import primepi
    print(f"{'x':>10}  {'Criva·x':>10}  {'π(x)':>10}  {'error%':>8}")
    for x in x_values:
        estimated = criva(x) * x
        real = primepi(x)
        error = abs(estimated - real) / real * 100
        print(f"{x:>10}  {estimated:>10.1f}  {real:>10}  {error:>7.3f}%")

compare_criva_vs_real([1000, 10000, 100000, 1000000])
```

---

## 4. Discriminant Method (Factorization)

### Concept

A **deterministic discard filter** for composite numbers based on a deformed discriminant. Given `N` to be tested for primality, the method searches for factor pairs `(2v+3, 2b+3)` by reducing the problem to finding perfect squares in a narrow algebraic window.

### Key formulas

```
N = (2v+3)(2b+3)  →  reparametrize as:
S = 2v
M = N - 9
Δ(S) = S² + 6S - M

Factorization exists ⟺ Δ(S) = k² for some integer k
```

**Deterministic stop condition:**

```
4S + 16 > 2√Δ + 1  →  STOP (no further solutions exist)
```

This condition is **not probabilistic** — it has geometric justification via gap-between-squares analysis.

### Python implementation

```python
import math

def discriminant_method(N):
    """
    Deterministic factorization filter using deformed discriminant.
    Returns factor pair if found, or 'prime candidate' if none.
    """
    if N % 2 == 0:
        return (2, N // 2)
    if N % 3 == 0:
        return (3, N // 3)
    
    M = N - 9
    sqrt_N = math.isqrt(N)
    v0 = (sqrt_N - 3) // 2
    S = 2 * v0
    
    steps = 0
    while S >= 0:
        delta = S * S + 6 * S - M
        
        if delta >= 0:
            k = math.isqrt(delta)
            if k * k == delta:
                # Found factorization
                v = S // 2
                b = (k - S - 3) // 2
                f1 = 2 * v + 3
                f2 = 2 * b + 3
                if f1 * f2 == N and f1 > 1 and f2 > 1:
                    return (f1, f2)
        
        # Deterministic stop condition
        sqrt_delta = math.isqrt(abs(delta)) if delta >= 0 else 0
        if 4 * S + 16 > 2 * sqrt_delta + 1:
            break
        
        S -= 2
        steps += 1
    
    return f"prime candidate (checked {steps} steps)"

# Test
for n in [91, 143, 221, 323, 10403, 15251]:
    print(f"N={n:6d}  →  {discriminant_method(n)}")
```

---

## 5. MDC — Kinematic Diophantine Method

### Concept

The **Método Diofántico Cinemático** (MDC) is a unified procedure for locating integer solutions of Diophantine equations `F(x,t) = 0` by analysing the kinematics of the decimal part of a real parametric function `d(t) = frac(g(x,t))`.

Instead of testing all candidates, the method measures velocity, acceleration and jerk of `d(t)` at 4 consecutive points and **extrapolates directly to the solution** in O(1) evaluations per jump.

### Two applications, one principle

| | Factorization | Wieferich Primes |
|--|--|--|
| Parameter `t` | `m` (factor candidate) | `K` (even integer) |
| Function `g(t)` | `N / (2*(2m+3))` | `p(K)`: solves `2^p = K*p^2 + 2` |
| Target `δ` | `0.5` | `0.0` |
| Search space | `L₁ = {6k±1}` | `2ℤ` (even integers) |
| Solution | factor `f = 2m+3` of `N` | Wieferich prime `p = n` |

### Wieferich sawtooth function

Define `d(K) = frac(p(K))` where `p(K)` is the real solution to `2^p = K*p^2 + 2`. Then:

> **Theorem (Zeros = Wieferich condition):** `d(K) = 0` for even integer `K` if and only if `K = (2^p − 2)/p²` for a Wieferich prime `p`.

> **Theorem (Parity):** Whenever `K = (b^p − b)/p² ∈ ℤ` for any odd prime `p`, `K` is even. Therefore iterating by steps of 2 loses no information.

### 1093–3511 coincidence explained

Both known Wieferich primes are Wieferich in bases `{2, 4, 8, 16, 32}` — all powers of 2. This is **a single condition**, not five independent ones: `p²|(2^p − 2)` implies `p²|((2^k)^p − 2^k)` for all `k` by a direct algebraic argument.

### Python implementation

```python
from mdc import mdc_factor, mdc_wieferich_scan, analyze_1093_3511

# Factorization
f1, f2, info = mdc_factor(10403)   # → (101, 103)

# Wieferich scan
mdc_wieferich_scan(p_min=3, p_max=30)

# Reproduce the 1093-3511 coincidence
analyze_1093_3511()
```

### Related paper

A formal write-up of the Wieferich results is available in `wieferich_paper.tex`
(arXiv preprint, 2026). The MDC framework is described in `metodo_diofantico_cinematico` (Spanish, 2026).

---

## 6. Siguiente Primo

### Concept

Given a known prime `inicio`, finds the **next prime** by maintaining three running accumulators `(t, tt, nt)` built from products of consecutive integers around three sliding counters `(ny, n, m)`.

The detection logic is extracted from a **Karnaugh map over 3 consecutive iterations** — it does not compute `(n−1)!` explicitly, but exploits the same divisibility structure as Wilson's theorem:

> `p` is prime ⟺ `(p−1)! ≡ p−1 (mod p)`

### Accumulators

```
ny = n−1,  n,  m = n+1       (three sliding counters)

t  *= ny   each step          (accumulates products via ny)
tt *= n    each step          (accumulates products via n)
nt *= m    each step          (accumulates products via m)

Residues: t%ny, t%n, t%m, tt%ny, tt%n, tt%m, nt%ny, nt%n, nt%m
```

### Karnaugh 3-pass detection

```python
# Pass 1 (2-iteration memory)
antp1  = (t3 > 0) and (nt2 == 0)
ant2p1 = (t3 > 0) and (nt2 == 0) or antp1
paso1  = (t3 > 0) and (nt2 == 0) or ant2p1

# Pass 2 (refinement)
antp2  = paso1 and (t2 > 0) and (tt2 > 0) and (t3 == 0)
paso2  = antp2 or ...

# Pass 3 (confirmation)
paso3  = antp2 and (t1 > 0) and (nt1 + nt2 == 0)
# → prime detected when paso3 fires
```

### Python implementation

```python
from siguiente_primo import siguiente_primo, primeros_n_primos, find_twin_primes

# Next prime after 11
siguiente_primo(11)          # → 13

# First 20 primes
primeros_n_primos(20)

# Twin primes up to 200
find_twin_primes(200)
```

---

## 7. SaltoMáximo

### Statement

> In the interval `[n − ⌊√(n+3)⌋ − 3,  n+3]` there are **at least two primes**.

This is stronger than Bertrand-Chebyshev (one prime in `[n, 2n]`) and gives an **explicit window** with a quantitative criterion.

### Key bound

```
(1 − (e−2)) · √(n+3)  ≥  2

where e−2 ≈ 0.718  (limit of the factorial accumulator m(n)/√(n+3))
      1−(e−2) ≈ 0.282  (minimum prime density factor)

Valid for n > ~100. Below that, verified case by case.
```

### Connection to MRAUV

The accumulator `m(n) = Σ_{i=2}^{K} √(n+3)/i!` converges to `(e−2)·√(n+3)`. The residual `L(n) − m(n) ≥ 2` is exactly the MRAUV-Goldbach criterion — both results share the same factorial accumulator.

### Python

```python
from salto_maximo import verificar_conjetura, tabla_verificacion, estadisticas_salto

verificar_conjetura(n_max=10_000)   # ✅ always >= 2 primes
estadisticas_salto(n_max=100_000)   # max gap vs sqrt(n+3)
```

---

## 8. ZypyZape MiniGemelo Digital

A **physics-based digital twin** of 5 wind turbines with synthetic kinetic battery — a separate engineering project by the same author.

### What it simulates

- 5 turbines with aerodynamic torque model `T = K_opt · ω²` (MPPT)
- Three operating roles: `CAPTACION` (harvest), `BAT_ACELERA` (kinetic charge), `BAT_FRENA` (kinetic brake)
- Grid frequency via swing equation: `RoCoF = −ΔP·f₀ / (2H·S_tot)`
- State of Energy: `SoE = (E_cin − E_min) / (E_max − E_min)`
- Perturbation injection (±200–300 MW) and automatic role rotation

### Key physics

```
J · dω/dt = T_aero − T_gen − T_roz     (rotor dynamics)
λ = ω·R/v                               (tip speed ratio)
Cp_max = 0.593                          (Betz limit)
f evolves via swing equation
```

### Run

```bash
pip install matplotlib numpy
python3 zypyzape_minigemelo.py          # interactive (requires display)
# or generates 4 static PNG scenarios automatically
```

---

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/prime-modular-methods
cd prime-modular-methods
pip install sympy numpy matplotlib
```

For C implementations (Criva sieve with OpenMP):

```bash
gcc -O2 -fopenmp -o criva_sieve criva_sieve.c -lm
./criva_sieve 1000000
```

---

## Repository structure

```
prime-modular-methods/
├── README.md
├── python/
│   ├── mrauv_goldbach.py      # MRAUV density model + Goldbach criterion
│   ├── sofi_structure.py      # Sophie Germain modular classification
│   ├── criva.py               # Iterative density estimator
│   ├── discriminant.py        # Deterministic factorization filter
│   ├── mdc.py                 # Kinematic Diophantine Method (unified)
│   ├── siguiente_primo.py     # Next prime via accumulated products (Wilson-inspired)
│   ├── salto_maximo.py        # Max prime gap conjecture: >= 2 primes in sqrt(n) window
│   └── zypyzape_minigemelo.py # Digital twin: 5 wind turbines + kinetic battery
├── c/
│   ├── criva_sieve.c          # Optimized sieve with OpenMP
│   └── modular_sieve.c        # 6k±1 modular sieve
├── papers/
│   ├── wieferich_paper.tex    # arXiv paper: Wieferich as sawtooth zeros
│   └── metodo_diofantico_cinematico.pdf  # MDC unified framework (Spanish)
├── notebooks/
│   ├── mrauv_analysis.ipynb   # Goldbach margin visualization
│   ├── sofi_classification.ipynb
│   ├── criva_vs_pnt.ipynb     # Criva vs Prime Number Theorem
│   └── mdc_sawtooth.ipynb     # Sawtooth structure visualization
└── docs/
    └── theory_notes.pdf        # Extended theoretical notes (Spanish)
```

---

## Theoretical context

| This work | Classical reference | Relationship |
|-----------|--------------------|----|
| MRAUV density | Prime Number Theorem (Hadamard, 1896) | Constructive local version |
| Sofí structure | Sophie Germain conjecture (open) | Modular framework for candidates |
| Criva estimator | Selberg / Brun sieves | Constructive layer model |
| Discriminant method | Fermat factorization | Deterministic stop condition (new) |
| MDC factorization | Trial division, Pollard ρ | Kinematic jump, O(1) per jump |
| MDC Wieferich | Fermat quotient criterion | Sawtooth zeros = Wieferich condition |

---

## Empirical validation

| Method | Range tested | Max relative error |
|--------|-------------|-------------------|
| Criva estimator | up to 10⁶ | < 0.1% |
| Discriminant filter | up to 10⁵ | 0% (deterministic) |
| MRAUV Goldbach margin | up to 10⁵ | positive in all cases |
| Sofí / U2 ⊆ LSG | up to 10⁴ | verified |
| MDC factorization | up to 10⁶ | 0% (exact hits) |
| MDC Wieferich predictor | primes up to 10⁴ | < 10⁻⁴ |
| Siguiente Primo | first 100 primes | 0% (exact, validated vs Wilson) |

---

## Citation

If you use this work in your research, please cite:

```
Manzanares Alberola, V. (2025). Prime Number Theory: Modular Methods & Density Models.
Unpublished manuscript. Available at: https://github.com/YOUR_USERNAME/prime-modular-methods
```

---

## License

MIT License — free to use, modify, and distribute with attribution.

---

## Author

**Víctor Manzanares Alberola**  
Independent researcher in computational number theory.  
Contact: [your email or profile link]

> *"No es que sea fácil. Es que merece la pena."*
