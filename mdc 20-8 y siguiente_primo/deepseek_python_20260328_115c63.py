import math
import time
from decimal import Decimal, getcontext
import matplotlib.pyplot as plt

# Configuración de precisión
PREC = 80
getcontext().prec = PREC

# ------------------------------------------------------------
# Funciones base
# ------------------------------------------------------------
def compute_v(N, m):
    return Decimal(N) / Decimal(2*m + 3) - Decimal('1.5')

def frac(x):
    return x - x.to_integral_value()

def phi(N, m):
    v = compute_v(N, m)
    f = frac(v)
    return float(min(f, 1 - f))

def medir_tramo(N, m0, paso, puntos=32):
    """
    Toma 'puntos' mediciones en una dirección.
    Devuelve ms, fs, y derivadas discretas (vel, acc, jerk).
    """
    ms = [m0 + i*paso for i in range(puntos)]
    fs = [phi(N, m) for m in ms]
    # derivadas
    vel = [fs[i+1] - fs[i] for i in range(puntos-1)]
    acc = [vel[i+1] - vel[i] for i in range(puntos-2)]
    jerk = [acc[i+1] - acc[i] for i in range(puntos-3)] if len(acc) >= 2 else []
    return ms, fs, vel, acc, jerk

def firma_rama(vel, acc):
    """Descriptor de la rama: media de velocidad y aceleración"""
    if not vel or not acc:
        return (0.0, 0.0)
    vel_mean = sum(vel) / len(vel)
    acc_mean = sum(acc) / len(acc)
    return (vel_mean, acc_mean)

def misma_parabola(firma1, firma2, tol_vel=0.05, tol_acc=0.01):
    """Compara si dos ramas pertenecen a la misma parábola"""
    vel1, acc1 = firma1
    vel2, acc2 = firma2
    return abs(vel1 - vel2) < tol_vel and abs(acc1 - acc2) < tol_acc

def ajuste_parabolico(ms, fs):
    m1, m2, m3 = ms[0], ms[1], ms[2]
    f1, f2, f3 = fs[0], fs[1], fs[2]
    denom = (f3 - 2*f2 + f1)
    if abs(denom) < 1e-30:
        return m2
    delta = (f3 - f1) / (2 * denom)
    return m2 - delta

def es_mismo_diente(v_vals):
    v_min = min(v_vals)
    v_max = max(v_vals)
    return (v_max - v_min) < 1.0

# ------------------------------------------------------------
# Algoritmo con 64 vectores y espejo explícito
# ------------------------------------------------------------
def factorizar_64v(N, max_iter=5000, verbose=False):
    N_dec = Decimal(N)
    m0 = int((N_dec.sqrt() - 3) // 2)
    if m0 < 1:
        m0 = 1

    paso = 1
    hist = []  # (m, fase_min, paso, firma_coherencia)
    coherencia_espejo = 0  # contador de coherencia para salto

    for it in range(max_iter):
        # Medir ambos lados
        ms_f, fs_f, vel_f, acc_f, _ = medir_tramo(N, m0, paso, puntos=32)
        ms_b, fs_b, vel_b, acc_b, _ = medir_tramo(N, m0, -paso, puntos=32)

        # Detectar cambio de diente
        v_vals = [compute_v(N, m) for m in ms_f + ms_b]
        if not es_mismo_diente(v_vals):
            paso = max(1, paso // 2)
            if verbose:
                print(f"  [Cambio de diente] paso -> {paso}")
            continue

        # Comparar ramas para determinar si estamos en la misma parábola
        firma_f = firma_rama(vel_f, acc_f)
        firma_b = firma_rama(vel_b, acc_b)
        misma = misma_parabola(firma_f, firma_b)

        # Estimar mínimos
        m_min_f = ajuste_parabolico(ms_f, fs_f)
        m_min_b = ajuste_parabolico(ms_b, fs_b)

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

        # Validar factor
        v = compute_v(N, m_next)
        if abs(frac(v)) < Decimal('1e-25'):
            v_int = int(v)
            p = 2*m_next + 3
            q = 2*v_int + 3
            return p, q, it+1, hist

        # --- Lógica del espejo ---
        # Si ambas ramas pertenecen a la misma parábola y la fase está lejos de 0,
        # podemos saltar más porque el espejo garantiza que no hay factor en la zona opuesta.
        if misma and 0.2 < fase_min < 0.8:
            coherencia_espejo += 1
            if coherencia_espejo >= 2:   # confirmación de estabilidad
                # Salto grande: usamos la pendiente de la fase para estimar la distancia al siguiente mínimo
                # Usamos velocidad media como estimador de frecuencia
                vel_media = (abs(sum(vel_f)/len(vel_f)) + abs(sum(vel_b)/len(vel_b))) / 2
                if vel_media > 1e-6:
                    # La fase cambia lentamente -> podemos saltar más
                    salto_estimado = max(2, int(0.5 / vel_media))
                else:
                    salto_estimado = 2 * paso
                paso = min(salto_estimado, 10000)
                coherencia_espejo = 0
            else:
                # Esperamos a confirmar
                paso = max(1, paso)
        else:
            coherencia_espejo = 0
            paso = max(1, paso // 2)

        # Guardar historial
        hist.append((m0, fase_min, paso, misma))

        # Avanzar
        m0 = m_next

        if verbose and it % 200 == 0:
            print(f"Iter {it}: m={m0}, fase={fase_min:.6f}, paso={paso}, misma={misma}")

    return None, None, max_iter, hist

# ------------------------------------------------------------
# Función de prueba
# ------------------------------------------------------------
def probar(N, verbose=False):
    print(f"\n--- Factorizando N = {N} ---")
    start = time.time()
    p, q, it, hist = factorizar_64v(N, max_iter=5000, verbose=verbose)
    elapsed = time.time() - start
    if p:
        print(f"✓ Factores: {p} × {q}  | iteraciones: {it} | tiempo: {elapsed:.3f} s")
    else:
        print(f"✗ No encontrado tras {it} iteraciones, {elapsed:.3f} s")
    return p, q, it, elapsed, hist

# ------------------------------------------------------------
# Visualización del historial
# ------------------------------------------------------------
def graficar_historial(hist, N):
    if not hist:
        return
    iters = list(range(len(hist)))
    ms = [h[0] for h in hist]
    fases = [h[1] for h in hist]
    coherencias = [h[3] for h in hist]   # True/False

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    ax1.plot(iters, ms, 'b-', marker='o', markersize=3)
    ax1.set_ylabel('m')
    ax1.set_title(f'Evolución de m para N={N}')
    ax1.grid(True)

    ax2.plot(iters, fases, 'r-', marker='o', markersize=3)
    ax2.set_ylabel('Fase mínima')
    ax2.grid(True)

    ax3.fill_between(iters, 0, 1, where=coherencias, color='green', alpha=0.3, label='Misma parábola')
    ax3.set_xlabel('Iteración')
    ax3.set_ylabel('Coherencia')
    ax3.set_yticks([])
    ax3.legend()
    ax3.grid(True)

    plt.tight_layout()
    plt.show()

# ------------------------------------------------------------
# Ejemplo
# ------------------------------------------------------------
if __name__ == "__main__":
    # Números de prueba (ajusta según quieras)
    tests = [
        10403,          # 101×103
        12001,          # 11×1091? (en realidad 12001 = 11×1091)
        251953,         # 493×511
        10000000000000003,  # 17 cifras (ejemplo)
    ]
    for N in tests:
        p, q, it, t, hist = probar(N, verbose=False)
        if p:
            graficar_historial(hist, N)