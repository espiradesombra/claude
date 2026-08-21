"""
MDC - Método Diofántico Cinemático (Libro 5)
Autor: Víctor Manzanares Alberola
Asistencia de escritura: IA

Factorización por diente de sierra + pinza doble + invariante parabólico.

Principios:
  - Sin floats en aritmética crítica (enteros exactos)
  - Buffer de decimales guardado línea a línea en TXT
  - Pinza doble: descenso (→) y ascenso (←)
  - Invariante parabólico: segunda diferencia discreta
  - Tercera medición: detector de error de escala / sobrepaso
  - Salto diferencial: no +1 fijo, sino estimado por velocidad
  - Estructura preparada para N grandes (Python big int nativo)
"""

import os
import sys

# ─────────────────────────────────────────────
# CONFIGURACIÓN
# ─────────────────────────────────────────────
LOG_FILE   = "mdc_buffer.txt"
MAX_ITER   = 500_000
TOL_D2     = 1e-9   # tolerancia segunda diferencia (simetría parabólica)
PINZA_PASOS = 4     # puntos por dirección (mínimo 3, mejor 4)


# ─────────────────────────────────────────────
# 1. NÚCLEO: FRACCIÓN EXACTA
#    d(m) = frac(N / (2*(2m+3)))  →  (num, denom)
#    Cero exacto: num == 0  ⟺  factor encontrado
# ─────────────────────────────────────────────
def frac_exacta(N: int, m: int):
    """
    d(m) = frac(N / (2m+3))
    Devuelve (num, denom) con num = N mod (2m+3).
    num == 0  ⟺  (2m+3) divide N  ⟺  factor encontrado.
    """
    denom = 2*m + 3
    num   = N % denom
    return num, denom


def es_factor_exacto(N: int, m: int):
    """Comprobación directa sin ningún float."""
    num, denom = frac_exacta(N, m)
    trivial = (denom == 1 or denom == N)
    return (num == 0 and not trivial), denom


# ─────────────────────────────────────────────
# 2. BUFFER / LOG
#    Cada punto se escribe como: m,num,denom
#    Formato de línea estable para análisis externo.
# ─────────────────────────────────────────────
def log_init():
    with open(LOG_FILE, "w") as f:
        f.write("m,num,denom\n")

def log_punto(m: int, num: int, denom: int):
    with open(LOG_FILE, "a") as f:
        f.write(f"{m},{num},{denom}\n")


# ─────────────────────────────────────────────
# 3. MEDIR LISTA (buffer como lista de dígitos)
#    Devuelve lista de tuplas (m, num, denom)
#    Dirección:  paso = +1 (descenso) o -1 (ascenso)
# ─────────────────────────────────────────────
def medir_lista(N: int, m0: int, pasos: int, paso: int) -> list:
    lista = []
    m = m0
    for _ in range(pasos):
        if m < 1:
            break
        num, denom = frac_exacta(N, m)
        lista.append((m, num, denom))
        log_punto(m, num, denom)
        m += paso
    return lista


# ─────────────────────────────────────────────
# 4. VALOR DECIMAL DE UN PUNTO (float local)
#    Solo se usa para métricas comparativas,
#    nunca para la decisión de factor.
# ─────────────────────────────────────────────
def val(p) -> float:
    _, num, denom = p
    return num / denom


# ─────────────────────────────────────────────
# 5. SEGUNDA DIFERENCIA DISCRETA
#    Para una parábola exacta: Δ²y = constante
#    Usamos 3 puntos consecutivos de la lista.
# ─────────────────────────────────────────────
def segunda_diferencia(lista: list, inicio: int = 0) -> float:
    y0 = val(lista[inicio])
    y1 = val(lista[inicio + 1])
    y2 = val(lista[inicio + 2])
    return y2 - 2*y1 + y0


# ─────────────────────────────────────────────
# 6. TERCERA MEDICIÓN: DETECTOR DE SOBREPASO
#    Si al añadir el punto 3, la segunda diferencia
#    cambia bruscamente → escala mal calibrada.
#    Devuelve True si el paso es correcto.
# ─────────────────────────────────────────────
def tercera_medicion_ok(lista: list, tol: float = 0.10) -> bool:
    if len(lista) < 4:
        return True
    d2_01 = segunda_diferencia(lista, 0)   # puntos 0,1,2
    d2_12 = segunda_diferencia(lista, 1)   # puntos 1,2,3
    return abs(d2_01 - d2_12) < tol


# ─────────────────────────────────────────────
# 7. MISMA PARÁBOLA (pinza doble)
#    Compara segunda diferencia en descenso y ascenso.
#    Si coinciden dentro de tolerancia → misma curva.
# ─────────────────────────────────────────────
def misma_parabola(desc: list, asc: list) -> bool:
    if len(desc) < 3 or len(asc) < 3:
        return False
    d2_desc = segunda_diferencia(desc)
    d2_asc  = segunda_diferencia(asc)
    return abs(d2_desc - d2_asc) < TOL_D2


# ─────────────────────────────────────────────
# 8. VELOCIDAD MEDIA (discreta)
#    Cambio de valor decimal entre primero y último.
# ─────────────────────────────────────────────
def velocidad_media(lista: list) -> float:
    if len(lista) < 2:
        return 0.0
    dy = val(lista[-1]) - val(lista[0])
    dm = lista[-1][0] - lista[0][0]
    if dm == 0:
        return 0.0
    return dy / dm


# ─────────────────────────────────────────────
# 9. SALTO DIFERENCIAL
#    Estimación de cuántos pasos saltar hacia el cero.
#    Basado en velocidad media discreta.
#    Clampado para evitar saltos fuera de rango.
# ─────────────────────────────────────────────
def calcular_salto(desc: list, m_max: int) -> int:
    v = velocidad_media(desc)
    if v == 0.0:
        return 1
    y0 = val(desc[0])
    # cuántos pasos hasta que y ≈ 0 extrapolando linealmente
    estimado = int(abs(y0 / v))
    estimado = max(1, min(estimado, m_max // 4))
    return estimado


# ─────────────────────────────────────────────
# 10. VÉRTICE PARABÓLICO (tres puntos exactos)
#     Resuelve sistema 3x3 entero, devuelve vértice en float.
# ─────────────────────────────────────────────
def vertice_parabola(lista: list):
    if len(lista) < 3:
        return None
    x1, y1 = lista[0][0], val(lista[0])
    x2, y2 = lista[1][0], val(lista[1])
    x3, y3 = lista[2][0], val(lista[2])

    denom = (x1-x2)*(x1-x3)*(x2-x3)
    if denom == 0:
        return None

    a = (y1*(x2-x3) + y2*(x3-x1) + y3*(x1-x2)) / denom
    b = (y1*(x3**2-x2**2) + y2*(x1**2-x3**2) + y3*(x2**2-x1**2)) / denom

    if a == 0:
        return None

    xv = -b / (2*a)
    return xv


# ─────────────────────────────────────────────
# 11. ALGORITMO PRINCIPAL
# ─────────────────────────────────────────────
def factorizar(N: int, verbose: bool = True) -> int | None:
    """
    Factoriza N usando el Método Diofántico Cinemático (MDC).
    Retorna un factor no trivial, o None si no encontró.

    Estrategia de inicio:
      - Empieza desde sqrt(N) hacia abajo (factores equilibrados primero).
      - Si no encuentra, sube desde m=1 (factores desbalanceados).
    """
    if N % 2 == 0:
        return 2
    if N % 3 == 0:
        return 3

    m_max = (N // 6) + 1
    log_init()

    # Calcular m de inicio desde sqrt(N)
    isqrt_N = int(N**0.5)
    m_start = max(1, (isqrt_N - 3) // 2)

    # Dos pasadas: desde sqrt bajando, luego desde 1 subiendo
    rangos = [
        range(m_start, 0, -1),               # desde sqrt hacia 1
        range(m_start + 1, min(MAX_ITER, m_max)),  # desde sqrt hacia arriba
    ]

    for rango in rangos:
        m_lista = list(rango)
        idx = 0
        while idx < len(m_lista):
            m = m_lista[idx]
            if m < 1:
                idx += 1
                continue

            # ── COMPROBACIÓN DIRECTA ──────────────────────
            ok, denom = es_factor_exacto(N, m)
            if ok:
                if verbose:
                    print(f"[EXACTO] Factor en m={m} → {denom}")
                return denom

            # ── MEDICIÓN DESCENSO ─────────────────────────
            paso = 1 if rango.step > 0 else -1
            desc = medir_lista(N, m, PINZA_PASOS, paso)

            # ── TERCERA MEDICIÓN: control de escala ───────
            if not tercera_medicion_ok(desc):
                idx += 1
                continue

            # ── MEDICIÓN ASCENSO (dirección opuesta) ──────
            asc = medir_lista(N, m, PINZA_PASOS, -paso)
            asc.reverse()

            # ── COMPROBACIÓN EXACTA EN LAS LISTAS ────────
            for p in desc + asc:
                mi, num, denom = p
                if num == 0 and denom > 1 and denom != N:
                    if verbose:
                        print(f"[LISTA] Factor en m={mi} → {denom}")
                    return denom

            # ── MISMA PARÁBOLA? ───────────────────────────
            if misma_parabola(desc, asc):
                xv = vertice_parabola(desc)
                if xv is not None:
                    mv = int(round(xv))
                    if 1 <= mv < m_max:
                        ok, denom = es_factor_exacto(N, mv)
                        if ok:
                            if verbose:
                                print(f"[VÉRTICE] Factor en m={mv} → {denom}")
                            return denom
                dm = calcular_salto(desc, m_max)
                idx += dm
            else:
                idx += 1

    if verbose:
        print(f"[FIN] No encontrado.")
    return None


# ─────────────────────────────────────────────
# 12. FACTORIZACIÓN COMPLETA
#     Recupera todos los factores primos.
# ─────────────────────────────────────────────
def factorizacion_completa(N: int) -> list:
    """Factorización completa recursiva en primos."""
    if N <= 1:
        return []
    if es_primo_simple(N):
        return [N]
    f = factorizar(N, verbose=False)
    if f is None or f == N or f == 1:
        return [N]
    return sorted(factorizacion_completa(f) + factorizacion_completa(N // f))


def es_primo_simple(n: int) -> bool:
    """Criba básica hasta sqrt(n). Solo para números pequeños."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i*i <= n:
        if n % i == 0 or n % (i+2) == 0:
            return False
        i += 6
    return True


# ─────────────────────────────────────────────
# 13. BATERÍA DE TESTS
# ─────────────────────────────────────────────
CASOS = [
    (143,             [11, 13]),
    (221,             [13, 17]),
    (323,             [17, 19]),
    (10403,           [101, 103]),
    (15251,           [101, 151]),
    (46189,           [11, 13, 17, 19]),
    (999961,          [999961]),           # primo
    (1_000_003,       [1_000_003]),        # primo
    (1_000_033,       [1000033]),          # verificar
    (104_729 * 104_723, [104723, 104729]), # semiprimo equilibrado
    (9999991,         [9999991]),          # primo
]

def ejecutar_tests():
    import time
    print("=" * 70)
    print("MDC - Batería de tests")
    print("=" * 70)
    for N, esperado in CASOS:
        t0 = time.perf_counter()
        factores = factorizacion_completa(N)
        t1 = time.perf_counter()
        ms = (t1-t0)*1000
        ok = "✓" if sorted(factores) == sorted(esperado) else "✗"
        print(f"  {ok}  N={N:>16,}  →  {factores}  ({ms:.2f} ms)")
    print("=" * 70)


# ─────────────────────────────────────────────
# 14. MODO INTERACTIVO
# ─────────────────────────────────────────────
def modo_interactivo():
    print("\nMDC factorizador — modo interactivo")
    print("Escribe un entero (o 'q' para salir)\n")
    while True:
        entrada = input("N = ").strip()
        if entrada.lower() in ("q", "exit", "salir"):
            break
        try:
            N = int(entrada)
            if N < 4:
                print("Introduce un entero >= 4")
                continue
            import time
            t0 = time.perf_counter()
            factores = factorizacion_completa(N)
            t1 = time.perf_counter()
            print(f"  Factores: {factores}  ({(t1-t0)*1000:.3f} ms)\n")
        except ValueError:
            print("  No es un entero válido\n")


# ─────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == "test":
            ejecutar_tests()
        else:
            try:
                N = int(arg)
                factores = factorizacion_completa(N)
                print(f"N={N}  →  {factores}")
            except ValueError:
                print("Uso: python mdc_libro5.py [N | test]")
    else:
        ejecutar_tests()
        print()
        modo_interactivo()
