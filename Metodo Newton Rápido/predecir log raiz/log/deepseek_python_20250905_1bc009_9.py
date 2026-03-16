def evaluacion_rapida(f1, f2, f4):
    """Evaluación rápida de solo 3 parámetros"""
    metricas = []
    for a, b in casos_rapidos:
        try:
            resultado, iteraciones, _ = log_doble_ajuste_simultaneo(a, b, f1, f2, F3_FIJO, f4, F5_FIJO, F6_FIJO)
            real = math.log(a, b)
            error = abs(resultado - real)
            metricas.append({
                'error': error,
                'iteraciones': iteraciones,
                'valido': error < 1e-10
            })
        except:
            metricas.append({'error': 100, 'iteraciones': 50, 'valido': False})
    
    # Calcular puntuación combinada
    errores = [m['error'] for m in metricas if m['valido']]
    if not errores:
        return float('inf')
    
    return np.mean(errores) + np.mean([m['iteraciones'] for m in metricas])/20