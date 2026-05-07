import math
import time
from decimal import Decimal, getcontext
import matplotlib.pyplot as plt

# Configuración de precisión (ajustable)
PREC = 80
getcontext().prec = PREC

# ------------------------------------------------------------
# Funciones base
# ------------------------------------------------------------
def compute_v(N, m):
    """v = N/(2m+3) - 1.5"""
    return Decimal(N) / Decimal(2*m + 3) - Decimal('1.5')

def frac(x):
    return x - x.to_integral_value()

def phi(N, m):
    """Fase: distancia al entero más cercano (diente de sierra)"""
    v = compute_v(N, m)
    f = frac(v)
    return float(min(f, 1 - f))

def medir_tramo(N, m0, paso, puntos=16):
    """
    Toma 'puntos' mediciones en una dirección.
    Devuelve listas de m, fase, y derivadas discretas.
    """
    ms = [m0 + i*paso for i in range(puntos)]
    fs = [phi(N, m) for m in ms]
    # derivadas
    vel = [fs[i+1] - fs[i] for i in range(puntos-1)]
    acc = [vel[i+1] - vel[i] for i in range(puntos-2)]
    jerk = [acc[i+1] - acc[i] for i in range(puntos-3)] if len(acc) >= 2 else []
    return ms, fs, vel, acc, jerk

def ajuste_parabolico(ms, fs):
    """
    Usa los primeros 3 puntos para estimar el mínimo de la parábola local.
    """
    m1, m2, m3 = ms[0], ms[1], ms[2]
    f1, f2, f3 = fs[0], fs[1], fs[2]
    denom = (f3 - 2*f2 + f1)
    if abs(denom) < 1e-30:
        return m2
    delta = (f3 - f1) / (2 * denom)
    return m2 - delta

def es_mismo_diente(v_vals):
    """Comprueba si todos los v caen dentro del mismo intervalo de 1 unidad"""
    v_min = min(v_vals)
    v_max = max(v_vals)
    return (v_max - v_min) < 1.0

# ------------------------------------------------------------
# Algoritmo principal con 32 vectores (16+16)
# ------------------------------------------------------------
def factorizar_32v(N, max_iter=2000, verbose=False):
    """
    Factoriza N usando 16 mediciones en cada dirección (32 vectores total),
    con detección de cambio de diente y simetría espejo.
    Devuelve (p, q, iteraciones, historial).
    """
    N_dec = Decimal(N)
    # punto inicial cerca de √N
    m0 = int((N_dec.sqrt() - 3) // 2)
    if m0 < 1:
        m0 = 1

    paso = 1
    hist = []       # (m, fase_min, salto, direccion)

    for it in range(max_iter):
        # --- Medir adelante y atrás ---
        ms_f, fs_f, vel_f, acc_f, jerk_f = medir_tramo(N, m0, paso, puntos=16)
        ms_b, fs_b, vel_b, acc_b, jerk_b = medir_tramo(N, m0, -paso, puntos=16)

        # --- Detectar cambio de diente ---
        v_vals = [compute_v(N, m) for m in ms_f + ms_b]
        if not es_mismo_diente(v_vals):
            # Salto de rama: reiniciamos con paso más fino
            paso = max(1, paso // 2)
            if verbose:
                print(f"  [Cambio de diente] paso -> {paso}")
            continue

        # --- Estimar mínimo por parábola (usamos los 3 primeros puntos de cada dirección) ---
        m_min_f = ajuste_parabolico(ms_f, fs_f)
        m_min_b = ajuste_parabolico(ms_b, fs_b)

        # Elegir la dirección con la fase más baja
        min_f = min(fs_f)
        min_b = min(fs_b)
        if min_f < min_b:
            m_next = int(round(m_min_f))
            fase_min = min_f
            direccion = 1
        else:
            m_next = int(round(m_min_b))
            fase_min = min_b
            direccion = -1

        # --- Validar si es factor ---
        v = compute_v(N, m_next)
        if abs(frac(v)) < Decimal('1e-25'):
            v_int = int(v)
            p = 2*m_next + 3
            q = 2*v_int + 3
            return p, q, it+1, hist

        # --- Simetría espejo: tachar espacio opuesto (aniquilación) ---
        # En la práctica, esto se traduce en un salto más agresivo si la fase es estable
        # Aquí lo simulamos aumentando el paso si la curvatura es pequeña y la fase central está lejos de 0
        curv = sum(acc_f)/len(acc_f) if acc_f else 0
        if abs(curv) < 1e-10 and 0.2 < fase_min < 0.8:
            # Zona estable → salto grande (espejo activo)
            paso = min(paso * 2, 1000)
        else:
            # Ajuste fino
            paso = max(1, paso // 2)

        # --- Guardar historial ---
        hist.append((m0, fase_min, paso, direccion))

        # --- Avanzar ---
        m0 = m_next

        if verbose and it % 100 == 0:
            print(f"Iter {it}: m={m0}, fase={fase_min:.6f}, paso={paso}")

    return None, None, max_iter, hist

# ------------------------------------------------------------
# Función de prueba y comparativa
# ------------------------------------------------------------
def probar_metodo(N, metodo='32v', verbose=False):
    start = time.time()
    if metodo == '32v':
        p, q, it, hist = factorizar_32v(N, max_iter=5000, verbose=verbose)
    else:
        raise ValueError("Método no implementado")
    elapsed = time.time() - start
    if p:
        print(f"[{metodo}] N={N} → {p} × {q}  en {it} iteraciones, {elapsed:.3f} s")
        return p, q, it, elapsed, hist
    else:
        print(f"[{metodo}] No encontrado en {it} iteraciones, {elapsed:.3f} s")
        return None, None, it, elapsed, hist

def comparar(N):
    print("\n" + "="*60)
    print(f"Comparativa para N = {N}")
    print("="*60)
    # Versión con 32 vectores
    probar_metodo(N, '32v', verbose=False)
    # Aquí podrías añadir la versión con 8 vectores si la tienes
    # probar_metodo(N, '8v', verbose=False)

# ------------------------------------------------------------
# Ejecución de ejemplo
# ------------------------------------------------------------
if __name__ == "__main__":
    # Números de prueba (hasta 14-15 cifras)
    tests = [
        10403,              # 101 × 103
        12001,              # 11 × 1091? (ajustar)
        251953,             # 493 × 511
        10000000000000003,  # semiprimo de 17 cifras (ejemplo)
        9999999999999997,   # otro
    ]
    for N in tests:
        comparar(N)
        # Puedes activar verbose para ver detalle
        # probar_metodo(N, '32v', verbose=True)