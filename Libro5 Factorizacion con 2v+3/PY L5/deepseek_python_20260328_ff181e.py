from decimal import Decimal, getcontext
import math
import time

# Configurar precisión alta (ajustar según necesidad)
PRECISION = 100
getcontext().prec = PRECISION

def compute_v(N, m):
    """v = (N/(2m+3)) - 1.5, con Decimal"""
    return Decimal(N) / Decimal(2*m + 3) - Decimal('1.5')

def frac(x):
    """Parte decimal de x (Decimal)"""
    return x - x.to_integral_value()

def phi(N, m):
    """Señal de fase: distancia al entero más cercano"""
    v = compute_v(N, m)
    f = frac(v)
    return min(f, 1 - f)

def measure_block(N, m, step, puntos=16):
    """
    Toma 'puntos' mediciones en una dirección.
    Devuelve listas de m y fase.
    """
    ms = []
    ys = []
    for i in range(puntos):
        mi = m + i * step
        ms.append(mi)
        ys.append(float(phi(N, mi)))
    return ms, ys

def ajuste_parabolico(ms, ys):
    """
    Ajusta una parábola a los primeros 3 puntos y estima el mínimo.
    Se usa como predictor local.
    """
    m1, m2, m3 = ms[:3]
    y1, y2, y3 = ys[:3]
    denom = (y3 - 2*y2 + y1)
    if abs(denom) < 1e-30:
        return m2
    delta = (y3 - y1) / (2 * denom)
    return m2 - delta

def derivadas_locales(ys):
    """
    Calcula primeras derivadas (velocidad, aceleración, jerk, snap)
    a partir de una lista de valores de fase.
    """
    vel = [ys[i+1] - ys[i] for i in range(len(ys)-1)]
    acc = [vel[i+1] - vel[i] for i in range(len(vel)-1)]
    jerk = [acc[i+1] - acc[i] for i in range(len(acc)-1)]
    snap = [jerk[i+1] - jerk[i] for i in range(len(jerk)-1)] if len(jerk) > 1 else []
    return vel, acc, jerk, snap

def misma_parabola(vel_f, acc_f, vel_b, acc_b, tol=1e-8):
    """
    Compara curvatura (aceleración) en ambas direcciones para decidir si
    estamos en la misma parábola.
    """
    # Usamos la media de las aceleraciones como medida
    curv_f = sum(acc_f) / len(acc_f) if acc_f else 0
    curv_b = sum(acc_b) / len(acc_b) if acc_b else 0
    return abs(curv_f - curv_b) < tol

def salto_adaptativo(N, m, paso, puntos=16):
    """
    Realiza una medición en ambas direcciones con 'puntos' vectores,
    calcula derivadas, decide si misma parábola, y devuelve el siguiente m
    candidato y la fase mínima observada.
    """
    # Medir en ambas direcciones
    ms_f, ys_f = measure_block(N, m, paso, puntos)
    ms_b, ys_b = measure_block(N, m, -paso, puntos)

    # Derivadas
    vel_f, acc_f, jerk_f, _ = derivadas_locales(ys_f)
    vel_b, acc_b, jerk_b, _ = derivadas_locales(ys_b)

    # Comprobar si estamos en la misma parábola
    misma = misma_parabola(vel_f, acc_f, vel_b, acc_b)

    # Elegir dirección basada en el mínimo observado
    min_f = min(ys_f)
    min_b = min(ys_b)
    if min_f < min_b:
        ms_sel = ms_f
        ys_sel = ys_f
        vel_sel = vel_f
        acc_sel = acc_f
        direccion = 1
    else:
        ms_sel = ms_b
        ys_sel = ys_b
        vel_sel = vel_b
        acc_sel = acc_b
        direccion = -1

    # Estimar mínimo local con ajuste parabólico (mejor que tomar el mínimo muestreado)
    m_est = ajuste_parabolico(ms_sel, ys_sel)
    m_next = int(round(m_est))

    # Si estamos en la misma parábola y la fase es estable, podemos saltar más
    # Usamos la aceleración media como indicador de curvatura
    curv = sum(acc_sel) / len(acc_sel) if acc_sel else 0

    # Calcular el nuevo paso (heuristico)
    if misma and abs(curv) < 1e-12 and min(ys_sel) > 0.01:
        # Zona estable: aumentar paso (pero limitado)
        nuevo_paso = min(paso * 2, 1000)
    else:
        # Zona de cambio o cerca de un mínimo: paso fino
        nuevo_paso = max(1, paso // 2)

    return m_next, min(ys_sel), nuevo_paso, misma, curv

def factorizar_con_32_vectores(N, max_iter=10000, precision_inicial=PRECISION):
    """
    Intenta factorizar N usando 32 vectores por iteración (16+16).
    Devuelve (factor1, factor2, iteraciones) o None.
    """
    getcontext().prec = precision_inicial
    N_dec = Decimal(N)

    # Punto inicial cerca de la raíz
    m = int((N_dec.sqrt() - 3) // 2)
    if m < 0:
        m = 1

    paso = 1
    for it in range(max_iter):
        m_next, fase_min, paso, misma, curv = salto_adaptativo(N, m, paso, puntos=16)

        # Verificar si es factor exacto
        v = compute_v(N, m_next)
        if abs(frac(v)) < Decimal('1e-25'):
            v_int = int(v)
            A = 2*m_next + 3
            B = 2*v_int + 3
            if A * B == N:
                return A, B, it+1

        # Actualizar posición
        m = m_next

        # Mostrar progreso cada 100 iteraciones
        if it % 100 == 0:
            print(f"Iter {it}: m={m}, fase={fase_min:.6f}, paso={paso}, misma={misma}, curv={curv:.2e}")

    return None

# ---------- Ejemplo de uso ----------
if __name__ == "__main__":
    # Números de prueba (semiprimos)
    # 17 cifras: 10000000000000003 * 10000000000000009
    N = 10000000000000003 * 10000000000000009
    print(f"Intentando factorizar N={N} (17 cifras)")
    inicio = time.time()
    res = factorizar_con_32_vectores(N)
    fin = time.time()
    if res:
        p, q, it = res
        print(f"Factores encontrados: {p} × {q} = {p*q}")
        print(f"Iteraciones: {it}, tiempo: {fin-inicio:.3f}s")
    else:
        print("No se encontraron factores.")