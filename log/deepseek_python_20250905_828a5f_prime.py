import math
import time

def log_doble_ajuste_simultaneo(a, b, factor1=1.3, factor2=3.5, factor3=-0.85, factor4=1.8, factor5=-0.85, factor6=1.8):
    """
    🏆 VERSIÓN CON DOBLE AJUSTE SIMULTÁNEO Y PARÁMETROS AJUSTABLES
    """
    if a <= 0 or b <= 0 or b == 1:
        raise ValueError("Base debe ser > 0 y ≠ 1, número debe ser > 0")
    
    j = 1.0
    d2, d1 = 0.0, 0.0
    precision = 1e-15
    historial = []
    
    for iteracion in range(15):
        historial.append(j)
        
        # Prevención de sobrepaso optimizada
        if a / (b ** j) >= b ** j:
            j *= factor1
            continue
            
        # Cálculo base
        j = (a + (a / (b ** j)) - 1) / a
        j1 = j
        d3, d2, d1 = d2, d1, j - j1
        
        # ¡AJUSTE DOBLE SIMULTÁNEO!
        # 1. Ajuste multiplicativo
        j *= (1 + factor2 * abs(d1))
        
        # 2. Ajuste exponencial
        if abs(d3 - d1) > 1e-12:
            factor = (factor3 + factor4 * (d3 - d2) / (d3 - d1))
            j = j ** abs(factor)
            
            # Cálculo base duplicado (intencional)
            j = (a + (a / (b ** j)) - 1) / a
            j1 = j
            d3, d2, d1 = d2, d1, j - j1
            factor = (factor5 + factor6 * (d3 - d2) / (d3 - d1))
            j = j ** abs(factor)
        
        if abs(j - j1) < precision:
            return j, iteracion + 1, historial
    
    return j, 15, historial

# Función para probar con diferentes parámetros
def probar_parametros(parametros, casos_prueba):
    resultados = []
    for params in parametros:
        tiempos = []
        errores = []
        iteraciones = []
        for a, b in casos_prueba:
            inicio = time.perf_counter_ns()
            res, it, _ = log_doble_ajuste_simultaneo(a, b, *params)
            fin = time.perf_counter_ns()
            tiempo = (fin - inicio) / 1000  # μs
            real = math.log(a, b)
            error = abs(res - real)
            
            tiempos.append(tiempo)
            errores.append(error)
            iteraciones.append(it)
        
        # Estadísticas promedio
        tiempo_prom = sum(tiempos) / len(tiempos)
        error_prom = sum(errores) / len(errores)
        iter_prom = sum(iteraciones) / len(iteraciones)
        
        resultados.append({
            'parametros': params,
            'tiempo_prom': tiempo_prom,
            'error_prom': error_prom,
            'iter_prom': iter_prom
        })
    
    return resultados

# Casos de prueba
casos_prueba = [(100, 10), (8, 2), (27, 3), (1024, 2), (10000, 10)]

# Generar combinaciones de parámetros a probar
parametros_a_probar = []

# Variar factor1 (1.3, 1.4, 1.5)
# Variar factor2 (3.5, 3.6, 3.7)
# Variar factor3 (-0.85, -0.80, -0.75)
for f1 in [1.3, 1.4, 1.5]:
    for f2 in [3.5, 3.6, 3.7]:
        for f3 in [-0.85, -0.80, -0.75]:
            parametros_a_probar.append((f1, f2, f3, 1.8, -0.85, 1.8))

# Ejecutar pruebas
resultados = probar_parametros(parametros_a_probar, casos_prueba)

# Mostrar mejores resultados ordenados por tiempo promedio
resultados_ordenados = sorted(resultados, key=lambda x: x['tiempo_prom'])
print("Mejores combinaciones de parámetros (por tiempo):")
for i, res in enumerate(resultados_ordenados[:5]):
    print(f"{i+1}. Parámetros: {res['parametros']}")
    print(f"   Tiempo: {res['tiempo_prom']:.2f}μs, Error: {res['error_prom']:.2e}, Iter: {res['iter_prom']:.1f}")
    print()