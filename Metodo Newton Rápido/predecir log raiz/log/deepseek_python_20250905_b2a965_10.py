# Solo 5x6x6 = 180 combinaciones (¡solo segundos!)
mejores_resultados = []

for f1 in parametros_criticos['factor1']:
    for f2 in parametros_criticos['factor2']:
        for f4 in parametros_criticos['factor4']:
            puntuacion = evaluacion_rapida(f1, f2, f4)
            mejores_resultados.append((puntuacion, f1, f2, f4))

# Ordenar y tomar las mejores
mejores_resultados.sort(key=lambda x: x[0])
top_5 = mejores_resultados[:5]

print("🔝 TOP 5 COMBINACIONES RÁPIDAS:")
for i, (punt, f1, f2, f4) in enumerate(top_5):
    print(f"{i+1}. f1={f1}, f2={f2}, f4={f4} → Puntuación: {punt:.6f}")