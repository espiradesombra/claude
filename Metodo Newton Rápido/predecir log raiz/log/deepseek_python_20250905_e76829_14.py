# PARÁMETROS CRÍTICOS A OPTIMIZAR
parametros_criticos = {
    'factor1': [1.1, 1.2, 1.3, 1.4, 1.5],      # Sobrepaso
    'factor2': [3.0, 3.2, 3.4, 3.6, 3.8, 4.0],  # Multiplicativo  
    'factor4': [1.5, 1.6, 1.7, 1.8, 1.9, 2.0]   # Exponencial
}

# Parámetros fijos (menos críticos)
F3_FIJO = -0.85
F5_FIJO = -0.85  
F6_FIJO = 1.8

# Casos de prueba rápidos pero representativos
casos_rapidos = [(100, 10), (8, 2), (1024, 2), (math.pi, 10)]

# Función de evaluación optimizada
def evaluacion_rapida(f1, f2, f4):
    metricas = []
    for a, b in casos_rapidos:
        try:
            # TU ALGORITMO con parámetros actuales
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
                
                # AJUSTE DOBLE
                j *= (1 + f2 * abs(d1))  # Ajuste multiplicativo
                
                if abs(d3 - d1) > 1e-12:
                    factor = (F3_FIJO + f4 * (d3 - d2) / (d3 - d1))
                    j = j ** abs(factor)
                    
                if abs(j - j1) < precision:
                    break
                    
            real = math.log(a, b)
            error = abs(j - real)
            metricas.append(error + iteracion/10)  # Ponderación error + iteraciones
            
        except:
            metricas.append(1000)  # Penalización por error
            
    return np.mean(metricas)

# EJECUTANDO OPTIMIZACIÓN RÁPIDA
print("🔍 Optimizando parámetros críticos...")
mejores_resultados = []

for f1 in parametros_criticos['factor1']:
    for f2 in parametros_criticos['factor2']:
        for f4 in parametros_criticos['factor4']:
            puntuacion = evaluacion_rapida(f1, f2, f4)
            mejores_resultados.append((puntuacion, f1, f2, f4))
            print(f"✓ Probado: f1={f1}, f2={f2}, f4={f4} → Puntuación: {puntuacion:.4f}")

# RESULTADOS
mejores_resultados.sort(key=lambda x: x[0])
top_5 = mejores_resultados[:5]

print("\n🎯 TOP 5 COMBINACIONES ÓPTIMAS:")
print("=" * 50)
for i, (punt, f1, f2, f4) in enumerate(top_5):
    print(f"{i+1}. f1={f1}, f2={f2}, f4={f4} → Score: {punt:.6f}")

# MEJOR COMBINACIÓN
MEJOR_F1, MEJOR_F2, MEJOR_F4 = top_5[0][1], top_5[0][2], top_5[0][3]
print(f"\n🏆 MEJOR COMBINACIÓN: f1={MEJOR_F1}, f2={MEJOR_F2}, f4={MEJOR_F4}")