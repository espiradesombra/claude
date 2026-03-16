def evaluar_parametros(f1, f2, f3, f4, f5, f6):
    """Evalúa un conjunto de parámetros y devuelve su puntuación"""
    
    resultados = []
    for a, b in casos_prueba:
        try:
            resultado, iteraciones, _ = log_doble_ajuste_simultaneo(a, b, f1, f2, f3, f4, f5, f6)
            real = math.log(a, b)
            error = abs(resultado - real)
            resultados.append({
                'error': error,
                'iteraciones': iteraciones,
                'convergio': error < 1e-10
            })
        except:
            resultados.append({'error': float('inf'), 'iteraciones': 50, 'convergio': False})
    
    # Calcular métricas combinadas
    errores = [r['error'] for r in resultados if r['convergio']]
    iteraciones = [r['iteraciones'] for r in resultados if r['convergio']]
    tasa_exito = sum(1 for r in resultados if r['convergio']) / len(resultados)
    
    if not errores:  # Si no convergió ningún caso
        return float('inf'), float('inf'), 0.0
    
    # Puntuación combinada (error promedio + iteraciones promedio - tasa de éxito)
    puntuacion = np.mean(errores) + np.mean(iteraciones)/10 - tasa_exito
    return puntuacion, np.mean(errores), tasa_exito