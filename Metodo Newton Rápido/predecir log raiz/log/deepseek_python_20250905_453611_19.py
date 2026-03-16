# USANDO LOS MEJORES PARÁMETROS DE LA FASE 1
MEJOR_F1 = 1.2
MEJOR_F2 = 3.4  
MEJOR_F4 = 1.9

# OPTIMIZAR PARÁMETROS SECUNDARIOS (3, 5, 6)
parametros_secundarios = {
    'factor3': [-0.95, -0.90, -0.85, -0.80, -0.75],
    'factor5': [-0.95, -0.90, -0.85, -0.80, -0.75],
    'factor6': [1.7, 1.8, 1.9, 2.0, 2.1]
}

# Función de evaluación para secundarios
def evaluacion_rapida_secundaria(f1, f2, f4, f3, f5, f6):
    metricas = []
    for a, b in casos_rapidos:
        try:
            j = 1.0
            d2, d1 = 0.0, 0.0
            precision = 1e-15
            
            for iteracion in range(15):
                if a / (b ** j) >= b ** j:
                    j *= f1
                    continue
                    
                j = (a + (a / (b ** j)) - 1) / a
                j1 = j
                d3, d2, d1 = d2, d1, j - j1
                
                # AJUSTE DOBLE CON NUEVOS PARÁMETROS
                j *= (1 + f2 * abs(d1))
                
                if abs(d3 - d1) > 1e-12:
                    factor = (f3 + f4 * (d3 - d2) / (d3 - d1))
                    j = j ** abs(factor)
                    
                if abs(j - j1) < precision:
                    break
                    
            real = math.log(a, b)
            error = abs(j - real)
            metricas.append(error + iteracion/10)
            
        except:
            metricas.append(1000)
            
    return np.mean(metricas)

print("🔍 Optimizando parámetros secundarios (3, 5, 6)...")
mejores_secundarios = []

total_combinaciones = len(parametros_secundarios['factor3']) * len(parametros_secundarios['factor5']) * len(parametros_secundarios['factor6'])
completadas = 0

for f3 in parametros_secundarios['factor3']:
    for f5 in parametros_secundarios['factor5']:
        for f6 in parametros_secundarios['factor6']:
            puntuacion = evaluacion_rapida_secundaria(MEJOR_F1, MEJOR_F2, MEJOR_F4, f3, f5, f6)
            mejores_secundarios.append((puntuacion, f3, f5, f6))
            completadas += 1
            print(f"✓ {completadas}/{total_combinaciones}: f3={f3}, f5={f5}, f6={f6} → Score: {puntuacion:.6f}")

# RESULTADOS FINALES
mejores_secundarios.sort(key=lambda x: x[0])
mejor_punt_sec, mejor_f3, mejor_f5, mejor_f6 = mejores_secundarios[0]

print("\n" + "="*60)
print("🎯 OPTIMIZACIÓN COMPLETA - RESULTADOS FINALES")
print("="*60)

print(f"\n🏆 PARÁMETROS ÓPTIMOS ENCONTRADOS:")
print(f"   • factor1 (sobrepaso): {MEJOR_F1}")
print(f"   • factor2 (multiplicativo): {MEJOR_F2}") 
print(f"   • factor3 (exponencial): {mejor_f3}")
print(f"   • factor4 (exponencial): {MEJOR_F4}")
print(f"   • factor5 (secundario): {mejor_f5}")
print(f"   • factor6 (secundario): {mejor_f6}")
print(f"   • Puntuación final: {mejor_punt_sec:.6f}")

# COMPARATIVA FINAL
punt_original = evaluacion_rapida_secundaria(1.3, 3.5, 1.8, -0.85, -0.85, 1.8)
mejora_porcentual = ((punt_original - mejor_punt_sec) / punt_original) * 100

print(f"\n📊 MEJORA TOTAL vs PARÁMETROS ORIGINALES:")
print(f"   Original: Score = {punt_original:.6f}")
print(f"   Optimizado: Score = {mejor_punt_sec:.6f}")
print(f"   ✅ MEJORA: {mejora_porcentual:.1f}%")

print(f"\n🚀 VALIDACIÓN EN CASOS DE PRUEBA:")
test_cases = [(100, 10), (8, 2), (1024, 2), (math.pi, 10), (1000000, 10), (0.001, 10)]
for a, b in test_cases:
    resultado, iteraciones, _ = log_doble_ajuste_simultaneo(a, b, MEJOR_F1, MEJOR_F2, mejor_f3, MEJOR_F4, mejor_f5, mejor_f6)
    real = math.log(a, b)
    error = abs(resultado - real)
    print(f"   log_{b}({a}) = {resultado:.8f} (error: {error:.2e}, iter: {iteraciones})")

print("\n" + "="*60)
print("¡OPTIMIZACIÓN FINALIZADA! 🎉")
print("="*60)