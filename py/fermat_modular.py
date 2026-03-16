"""
fermat_modular.py
=================
Férmines (números de Fermat) y alineación modular — VMA (Víctor Manzanares Alberola)

Concepto original:
    Los "Férmines" son los números de Fermat  Fₙ = 2^(2ⁿ) + 1.
    VMA descubrió dos propiedades clave:

    1. Producto menos dos:
           F₀ · F₁ · ... · Fₙ  =  Fₙ₊₁ − 2
       (identidad exacta, demostrable por inducción)

    2. Alineación modular con residuo privilegiado:
       · Base de alineación:  Bₙ = 2 × (F₀ · F₁ · ... · Fₙ)
       · Existe un único rₙ tal que:
             rₙ ≡ 1   (mod 2)
             rₙ ≡ 2   (mod Fₖ)   para todo k ≤ n
         y para todo m ≥ n+1:   Fₘ ≡ rₙ  (mod Bₙ)
       · En base Bₙ, los números de Fermat primos solo pueden
         aparecer en la "primera fila" (residuo rₙ).
         Alargando el productorio, la segunda mitad de la primera
         fila también queda descartada como compuesta.

Estructura del módulo
---------------------
1. fermin(n)                  → Fₙ = 2^(2ⁿ) + 1
2. producto_fermines(n)       → F₀ · F₁ · ... · Fₙ
3. verificar_identidad(n)     → comprueba F₀·...·Fₙ = Fₙ₊₁ − 2
4. base_alineacion(n)         → Bₙ = 2 × productorio
5. residuo_privilegiado(n)    → rₙ por CRT (≡1 mod 2, ≡2 mod Fₖ)
6. alineacion_modular(n,kmax) → verifica Fₘ ≡ rₙ (mod Bₙ) para m>n
7. primalidad_fermines(kmax)  → tabla de primalidad de F₀..Fₖₘₐₓ
8. tabla_alineacion(n)        → imprime residuos en base Bₙ
9. demo()                     → demostración completa en consola

Autoría
-------
Investigación original: Víctor Manzanares Alberola
Asistencia de escritura/código: IA
"""

import math
from functools import reduce


# ─────────────────────────────────────────────────────────────
# 1. FUNCIÓN BASE: Fₙ = 2^(2ⁿ) + 1
# ─────────────────────────────────────────────────────────────

def fermin(n: int) -> int:
    """
    Devuelve el n-ésimo número de Fermat Fₙ = 2^(2ⁿ) + 1.

    Ejemplos:
        F₀ = 3, F₁ = 5, F₂ = 17, F₃ = 257, F₄ = 65537
        F₅ = 4294967297 = 641 × 6700417  (compuesto, Euler 1732)
    """
    if n < 0:
        raise ValueError("n debe ser ≥ 0")
    return 2 ** (2 ** n) + 1


# ─────────────────────────────────────────────────────────────
# 2. PRODUCTORIO F₀ · F₁ · ... · Fₙ
# ─────────────────────────────────────────────────────────────

def producto_fermines(n: int) -> int:
    """
    Devuelve F₀ · F₁ · ... · Fₙ.
    Por la identidad VMA, esto es igual a Fₙ₊₁ − 2.
    """
    resultado = 1
    for k in range(n + 1):
        resultado *= fermin(k)
    return resultado


# ─────────────────────────────────────────────────────────────
# 3. VERIFICACIÓN DE LA IDENTIDAD F₀·...·Fₙ = Fₙ₊₁ − 2
# ─────────────────────────────────────────────────────────────

def verificar_identidad(n_max: int = 6) -> list[dict]:
    """
    Verifica la identidad  ∏_{k=0}^{n} Fₖ = Fₙ₊₁ − 2  para n = 0..n_max.

    Demostración por inducción:
        Base n=0: F₀ = 3 = F₁ − 2 = 5 − 2 ✓
        Paso: si ∏_{k=0}^{n} Fₖ = Fₙ₊₁ − 2,
              entonces ∏_{k=0}^{n+1} Fₖ = (Fₙ₊₁−2)·Fₙ₊₁
                                         = Fₙ₊₁² − 2·Fₙ₊₁
                                         = (2^(2^(n+1)))² + 2·2^(2^(n+1)) + 1
                                           − 2·(2^(2^(n+1)) + 1)
                                         = 2^(2^(n+2)) + 1 − 2
                                         = Fₙ₊₂ − 2  ✓
    """
    resultados = []
    for n in range(n_max + 1):
        prod = producto_fermines(n)
        fn1_menos2 = fermin(n + 1) - 2
        ok = (prod == fn1_menos2)
        resultados.append({
            'n': n,
            'F0_hasta_Fn': prod,
            'F(n+1)-2': fn1_menos2,
            'identidad': ok
        })
    return resultados


# ─────────────────────────────────────────────────────────────
# 4. BASE DE ALINEACIÓN Bₙ = 2 × ∏_{k=0}^{n} Fₖ
# ─────────────────────────────────────────────────────────────

def base_alineacion(n: int) -> int:
    """
    Devuelve Bₙ = 2 × F₀ × F₁ × ... × Fₙ = 2 × (Fₙ₊₁ − 2).

    Bₙ es el módulo de alineación: en base Bₙ, los Férmines
    futuros Fₘ (m > n) caen todos en el mismo residuo rₙ.
    """
    return 2 * producto_fermines(n)


# ─────────────────────────────────────────────────────────────
# 5. RESIDUO PRIVILEGIADO rₙ (CRT)
# ─────────────────────────────────────────────────────────────

def _crt_dos_modulos(r1, m1, r2, m2):
    """CRT para dos congruencias: x ≡ r1 (mod m1), x ≡ r2 (mod m2)."""
    # Extendido de Euclides
    g, s, _ = _egcd(m1, m2)
    if (r2 - r1) % g != 0:
        raise ValueError("Sistema sin solución (módulos no coprimos con residuos incompatibles)")
    lcm = m1 * m2 // g
    x = (r1 + m1 * ((r2 - r1) // g) * s) % lcm
    return x, lcm


def _egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = _egcd(b, a % b)
    return g, y, x - (a // b) * y


def residuo_privilegiado(n: int) -> tuple[int, int]:
    """
    Calcula el residuo privilegiado rₙ tal que:
        rₙ ≡ 1   (mod 2)
        rₙ ≡ 2   (mod Fₖ)   para k = 0, 1, ..., n

    Devuelve (rₙ, Bₙ) donde Bₙ es la base de alineación.

    Para todo m > n: Fₘ ≡ rₙ (mod Bₙ).
    """
    # Construimos el sistema de congruencias por CRT iterativo
    r, m = 1, 2  # rₙ ≡ 1 (mod 2)
    for k in range(n + 1):
        fk = fermin(k)
        r, m = _crt_dos_modulos(r, m, 2, fk)

    Bn = base_alineacion(n)
    rn = r % Bn
    return rn, Bn


# ─────────────────────────────────────────────────────────────
# 6. VERIFICACIÓN: Fₘ ≡ rₙ (mod Bₙ) para m > n
# ─────────────────────────────────────────────────────────────

def alineacion_modular(n: int, m_max: int = 8) -> list[dict]:
    """
    Verifica que Fₘ ≡ rₙ (mod Bₙ) para m = n+1, n+2, ..., m_max.

    Esto confirma la propiedad de alineación: todos los Férmines
    futuros caen en el mismo residuo dentro de la base Bₙ.
    """
    rn, Bn = residuo_privilegiado(n)
    resultados = []
    for m in range(n + 1, m_max + 1):
        Fm = fermin(m)
        residuo = Fm % Bn
        ok = (residuo == rn % Bn)
        resultados.append({
            'm': m,
            f'F{m} mod B{n}': residuo,
            f'r{n} mod B{n}': rn % Bn,
            'alineado': ok
        })
    return resultados


# ─────────────────────────────────────────────────────────────
# 7. TABLA DE PRIMALIDAD DE LOS FÉRMINES
# ─────────────────────────────────────────────────────────────

def primalidad_fermines(k_max: int = 6) -> list[dict]:
    """
    Tabla de primalidad para F₀ ... Fₖₘₐₓ.

    Conocido:
        F₀=3, F₁=5, F₂=17, F₃=257, F₄=65537 → primos
        F₅ = 641 × 6700417 → compuesto (Euler, 1732)
        F₆ en adelante → todos compuestos conocidos hasta F₃₂ (≈2024)
    """
    def es_primo_simple(k):
        """Criba básica (solo útil para k pequeño)."""
        if k < 2: return False
        if k == 2: return True
        if k % 2 == 0: return False
        for d in range(3, int(math.isqrt(k)) + 1, 2):
            if k % d == 0: return False
        return True

    resultados = []
    for k in range(k_max + 1):
        fk = fermin(k)
        # Para k ≤ 4 hacemos test directo; k ≥ 5 son números gigantes
        if k <= 4:
            primo = es_primo_simple(fk)
            metodo = 'test directo'
        elif k == 5:
            primo = False
            metodo = 'Euler 1732: 641 | F₅'
        else:
            primo = False
            metodo = 'compuesto (ningún primo Fermat conocido k≥5)'

        resultados.append({
            'k': k,
            'Fk': fk if k <= 5 else f"2^{2**k}+1 (~{2**(2**k-1):.2e} cifras)",
            'primo': primo,
            'bits': 2**k + 1,
            'nota': metodo
        })
    return resultados


# ─────────────────────────────────────────────────────────────
# 8. TABLA DE ALINEACIÓN EN BASE Bₙ
# ─────────────────────────────────────────────────────────────

def tabla_alineacion(n: int, m_max: int = 8) -> None:
    """
    Imprime la tabla de residuos Fₘ mod Bₙ para m = 0..m_max.
    Muestra en qué "fila" (residuo) cae cada Fermín.
    """
    rn, Bn = residuo_privilegiado(n)
    print(f"\nBase de alineación B{n} = {Bn}")
    print(f"Residuo privilegiado r{n} = {rn}")
    print(f"\n  {'m':>4}  {'Fₘ mod B'+str(n):>20}  {'es rₙ?':>8}  {'primo?':>7}")
    print("  " + "-" * 48)

    primos_conocidos = {0, 1, 2, 3, 4}  # F₀..F₄ son los únicos conocidos

    for m in range(m_max + 1):
        Fm = fermin(m)
        res = Fm % Bn
        es_rn = (res == rn % Bn)
        primo = m in primos_conocidos
        marca = "←" if es_rn else ""
        print(f"  {m:>4}  {res:>20}  {str(es_rn):>8}  {'SÍ' if primo else 'no':>7}  {marca}")


# ─────────────────────────────────────────────────────────────
# 9. DEMOSTRACIÓN COMPLETA
# ─────────────────────────────────────────────────────────────

def demo():
    print("=" * 65)
    print("FÉRMINES (NÚMEROS DE FERMAT) Y ALINEACIÓN MODULAR — VMA")
    print("=" * 65)

    # Valores básicos
    print("\n--- Los primeros Férmines Fₙ = 2^(2ⁿ) + 1 ---")
    for k in range(8):
        fk = fermin(k)
        print(f"  F{k} = 2^{2**k} + 1 = {fk}")

    # Identidad producto
    print("\n--- Identidad: F₀·F₁·...·Fₙ = Fₙ₊₁ − 2 ---")
    for r in verificar_identidad(6):
        n = r['n']
        p = r['F0_hasta_Fn']
        fn1 = r['F(n+1)-2']
        ok = "✓" if r['identidad'] else "✗"
        print(f"  n={n}: productorio = {p:>12}  |  F{n+1}-2 = {fn1:>12}  {ok}")

    # Primalidad
    print("\n--- Primalidad de los Férmines ---")
    for r in primalidad_fermines(6):
        p = "PRIMO" if r['primo'] else "compuesto"
        print(f"  F{r['k']:>1}  {p:>10}  ({r['nota']})")

    # Base y residuo privilegiado
    print("\n--- Base de alineación y residuo privilegiado ---")
    for n in range(4):
        rn, Bn = residuo_privilegiado(n)
        print(f"  n={n}: B{n} = {Bn:>10},  r{n} = {rn}")

    # Alineación modular
    print("\n--- Alineación: Fₘ ≡ r₁ (mod B₁) para m > 1 ---")
    tabla_alineacion(n=1, m_max=7)

    print("\n--- Alineación: Fₘ ≡ r₂ (mod B₂) para m > 2 ---")
    tabla_alineacion(n=2, m_max=7)

    # Interpretación de VMA
    print("\n--- Interpretación (VMA) ---")
    print("  En base B₁ = 2×F₀×F₁ = 30:")
    rn, Bn = residuo_privilegiado(1)
    print(f"  Todos los Fₘ con m≥2 caen en residuo r₁={rn} (mod {Bn}).")
    print("  Los Férmines primos solo pueden estar en esa 'primera fila'.")
    print("  Alargando el productorio (mayor n), la segunda mitad")
    print("  de esa fila también queda descartada como compuesta.")

    print("\n--- Verificación numérica de la alineación para n=2 ---")
    resultados = alineacion_modular(n=2, m_max=7)
    for r in resultados:
        ok = "✓" if r['alineado'] else "✗"
        print(f"  m={r['m']}: F{r['m']} mod B₂ = {list(r.values())[1]:>6}  {ok}")


if __name__ == "__main__":
    demo()
