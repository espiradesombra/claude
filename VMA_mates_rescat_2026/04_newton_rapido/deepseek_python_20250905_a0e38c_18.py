# USANDO LOS MEJORES PARÁMETROS ENCONTRADOS
MEJOR_F1 = 1.2
MEJOR_F2 = 3.4  
MEJOR_F4 = 1.9

# OPTIMIZAR PARÁMETROS SECUNDARIOS
parametros_secundarios = {
    'factor3': [-0.95, -0.90, -0.85, -0.80, -0.75],
    'factor5': [-0.95, -0.90, -0.85, -0.80, -0.75],
    'factor6': [1.7, 1.8, 1.9, 2.0, 2.1]
}

print("🔍 Optimizando parámetros secundarios...")
mejores_secundarios = []

for f3 in parametros_secundarios['factor3']:
    for f5 in parametros_secundarios['factor5']:
        for f6 in parametros_secundarios['factor6']:
            puntuacion = evaluacion_rapida_secundaria(MEJOR_F1, MEJOR_F2, MEJOR_F4, f3, f5, f6)
            mejores_secundarios.append((puntuacion, f3, f5, f6))
            print(f"✓ Probado: f3={f3}, f5={f5}, f6={f6} → Score: {puntuacion:.4f}")

# FUNCIÓN DE EVALUACIÓN PARA SECUNDARIOS
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

# RESULTADOS FINALES
mejores_secundarios.sort(key=lambda x: x[0])
mejor_f3, mejor_f5, mejor_f6 = mejores_secundarios[0][1], mejores_secundarios[0][2], mejores_secundarios[0][3]

print(f"\n🏆 OPTIMIZACIÓN COMPLETA!")
print(f"   factor1: {MEJOR_F1} (sobrepaso)")
print(f"   factor2: {MEJOR_F2} (multiplicativo)") 
print(f"   factor3: {mejor_f3} (exponencial)")
print(f"   factor4: {MEJOR_F4} (exponencial)")
print(f"   factor5: {mejor_f5} (secundario)")
print(f"   factor6: {mejor_f6} (secundario)")