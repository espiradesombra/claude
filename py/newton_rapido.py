"""
newton_rapido.py
================
Newton Rápido con MEcuations — VMA (Víctor Manzanares Alberola)

Algoritmo:
    Dado E > 0 y base b > 0, b ≠ 1, se busca j tal que b^j = E.

    Cada iteración aplica DOS pasos:

        1. MEDICIÓN (primer medidor)
               j ← (E + E·b^(-j) - 1) / E
           Acumula el historial de cuatro diferencias consecutivas:
               d1, d2, d3, d4

        2. SALTO EXPONENCIAL
               j ← j ^ |−0.85 + 1.8·(d3−d2)/(d3−d1)|
           Ajusta la curvatura de convergencia según la aceleración
           de las diferencias acumuladas en el paso anterior.

    Condición de parada: |j − j_anterior| < 10^(-12)

    Parámetros:
        fE  = 1.8    (factor exponencial)
        cte = -0.85  (constante de curvatura)
        sbp = 1.3    (factor de sobrepaso cuando b^j ≥ √E)

MEcuations como oráculos:
    Si E pertenece a una familia conocida, j₀ se obtiene directamente:
        cuadrados p²        → j₀ = ln(E) / 2
        cubos p³            → j₀ = ln(E) / 3
        potencia pⁿ         → j₀ = ln(E) / n
        k·p (k conocido)    → j₀ = ln(E/k)
        Mersenne 2ⁿ−1       → j₀ = n·ln(2) / ln(b)
        sin estructura      → j₀ = 1.0

Autoría:
    Investigación original: Víctor Manzanares Alberola
    Asistencia de escritura/código: IA
"""

import math


# ─────────────────────────────────────────────────────────────
# PARÁMETROS GLOBALES
# ─────────────────────────────────────────────────────────────

fE   =  1.8     # factor exponencial
cte  = -0.85    # constante de curvatura
sbp  =  1.3     # factor de sobrepaso
TOL  =  1e-12   # tolerancia de convergencia
ITER_MAX = 200  # máximo de iteraciones


# ─────────────────────────────────────────────────────────────
# ORÁCULOS MECUATION
# ─────────────────────────────────────────────────────────────

def oraculo(E: float, familia: str = "general", b: float = 10.0, **kwargs) -> float:
    """
    Devuelve j₀ según la familia numérica de E.

    Parámetros opcionales según familia:
        n  → potencia n (familia "potencia")
        k  → factor conocido (familia "kp")
    """
    lnE = math.log(E)
    lnb = math.log(b)

    if familia == "cuadrados":
        return lnE / 2 / lnb
    elif familia == "cubos":
        return lnE / 3 / lnb
    elif familia == "potencia":
        n = kwargs.get("n", 2)
        return lnE / n / lnb
    elif familia == "kp":
        k = kwargs.get("k", 1)
        return math.log(E / k) / lnb
    elif familia == "mersenne":
        # E = 2^n - 1  →  n ≈ log2(E+1)
        n = math.log2(E + 1)
        return n * math.log(2) / lnb
    else:
        return 1.0


# ─────────────────────────────────────────────────────────────
# NÚCLEO: NEWTON RÁPIDO
# ─────────────────────────────────────────────────────────────

def newton_rapido(E: float, b: float = 10.0, j0: float = None,
                  verbose: bool = False) -> dict:
    """
    Calcula log_b(E) mediante el algoritmo Newton Rápido.

    Pasos por iteración:
        1. Medición   j ← (E + E·b^(-j) - 1) / E
                      acumula d1, d2, d3, d4
        2. Exponencial j ← j ^ |cte + fE·(d3-d2)/(d3-d1)|

    Devuelve dict con: j, iteraciones, error, historial_j
    """
    if j0 is None:
        j0 = 1.0

    j = j0
    j_ref = math.log(E) / math.log(b)   # valor exacto para medir error

    # Historial de cuatro diferencias consecutivas
    d1 = d2 = d3 = d4 = 0.0

    historial = [j]
    j_prev = j

    for it in range(1, ITER_MAX + 1):

        # ── PASO 1: MEDICIÓN (primer medidor) ──────────────────
        try:
            bpow = b ** (-j)
        except OverflowError:
            bpow = 0.0

        j_nuevo = (E + E * bpow - 1.0) / E

        # Sobrepaso: si b^j ≥ √E, amortiguar
        try:
            if b ** j >= math.sqrt(E):
                j_nuevo = j_nuevo / sbp
        except (OverflowError, ValueError):
            pass

        # Acumular historial de diferencias
        diff = j_nuevo - j
        d4 = d3
        d3 = d2
        d2 = d1
        d1 = diff

        j = j_nuevo

        # ── PASO 2: SALTO EXPONENCIAL ──────────────────────────
        denom = d3 - d1
        if abs(denom) > 1e-30 and abs(d3 - d2) > 1e-30:
            exp_val = abs(cte + fE * (d3 - d2) / denom)
            if exp_val > 0 and j > 0:
                try:
                    j = j ** exp_val
                except (ValueError, OverflowError):
                    pass  # mantener j si el salto no es aplicable

        historial.append(j)

        if verbose:
            print(f"  it={it:3d}  j={j:.10f}  d1={d1:.2e}  d2={d2:.2e}"
                  f"  d3={d3:.2e}  d4={d4:.2e}  err={abs(j-j_ref):.2e}")

        # Condición de parada
        if abs(j - j_prev) < TOL:
            break
        j_prev = j

    error = abs(j - j_ref)
    return {
        "j":           j,
        "j_exacto":    j_ref,
        "iteraciones": it,
        "error":       error,
        "historial":   historial,
        "d1": d1, "d2": d2, "d3": d3, "d4": d4
    }


# ─────────────────────────────────────────────────────────────
# INTERFAZ CON ORÁCULO MECUATION
# ─────────────────────────────────────────────────────────────

def log_con_oraculo(E: float, b: float = 10.0, familia: str = "general",
                    verbose: bool = False, **kwargs) -> dict:
    """
    Newton Rápido con punto inicial proporcionado por MEcuation.

    Ejemplo:
        log_con_oraculo(E=289, b=10, familia="cuadrados")
        log_con_oraculo(E=4294967295, b=2, familia="mersenne")
    """
    j0 = oraculo(E, familia=familia, b=b, **kwargs)
    resultado = newton_rapido(E, b=b, j0=j0, verbose=verbose)
    resultado["familia"] = familia
    resultado["j0"] = j0
    return resultado


# ─────────────────────────────────────────────────────────────
# DEMO
# ─────────────────────────────────────────────────────────────

def demo():
    print("=" * 70)
    print("NEWTON RÁPIDO CON MECUATIONS — VMA")
    print("Medición (d1..d4) + Salto Exponencial  |  Sin ajuste multiplicativo")
    print("=" * 70)

    # ── Tabla 1: cuadrados p², base 10 ────────────────────────
    print("\n--- Familia cuadrados p² (base b=10) ---")
    print(f"  {'p':>4}  {'E=p²':>8}  {'log₁₀(E)':>10}  {'iter':>5}  {'error':>10}  {'j':>14}")
    print("  " + "-" * 58)

    primos = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    for p in primos:
        E = p * p
        r = log_con_oraculo(E, b=10, familia="cuadrados")
        print(f"  {p:>4}  {E:>8}  {r['j_exacto']:>10.4f}  "
              f"{r['iteraciones']:>5}  {r['error']:>10.2e}  {r['j']:>14.10f}")

    # ── Comparativa con/sin oráculo ────────────────────────────
    print("\n--- Comparativa: sin oráculo vs con oráculo (cuadrados) ---")
    print(f"  {'E':>10}  {'sin j₀ iter':>12}  {'sin j₀ err':>11}  "
          f"{'con j₀ iter':>12}  {'con j₀ err':>11}")
    print("  " + "-" * 62)

    casos = [4, 9, 25, 49, 121, 169, 289, 361, 529, 841]
    for E in casos:
        r_sin = newton_rapido(E, b=10, j0=1.0)
        r_con = log_con_oraculo(E, b=10, familia="cuadrados")
        print(f"  {E:>10}  {r_sin['iteraciones']:>12}  {r_sin['error']:>11.2e}  "
              f"{r_con['iteraciones']:>12}  {r_con['error']:>11.2e}")

    # ── Todas las familias ─────────────────────────────────────
    print("\n--- Oráculos MEcuation por familia ---")
    casos_fam = [
        (4,           10, "cuadrados", {},       "p=2, E=2²"),
        (8,           10, "cubos",     {},       "p=2, E=2³"),
        (1024,        10, "potencia",  {"n":10}, "E=2¹⁰"),
        (35,          10, "kp",        {"k":5},  "E=5·7, k=5"),
        (4294967295,   2, "mersenne",  {},       "E=2³²−1"),
        (77,          10, "general",   {},       "E=7·11, sin oráculo"),
    ]
    print(f"  {'caso':>25}  {'familia':>12}  {'j₀':>10}  {'iter':>5}  {'error':>10}")
    print("  " + "-" * 70)
    for E, b, fam, kw, desc in casos_fam:
        r = log_con_oraculo(E, b=b, familia=fam, **kw)
        print(f"  {desc:>25}  {fam:>12}  {r['j0']:>10.4f}  "
              f"{r['iteraciones']:>5}  {r['error']:>10.2e}")

    # ── Traza de d1..d4 ───────────────────────────────────────
    print("\n--- Traza de historial d1..d4 para E=121, b=10 (oráculo cuadrados) ---")
    newton_rapido(121, b=10, j0=oraculo(121, "cuadrados", 10), verbose=True)


if __name__ == "__main__":
    demo()
