# Primero: Ajustar factores 1-3 manteniendo 4-6 fijos
mejores_combinaciones = []

# Probaremos combinaciones de factores 1-3
for f1, f2, f3 in product(rangos['factor1'], rangos['factor2'], rangos['factor3']):
    puntuacion, error, exito = evaluar_parametros(f1, f2, f3, 1.8, -0.85, 1.8)
    mejores_combinaciones.append((puntuacion, error, exito, (f1, f2, f3, 1.8, -0.85, 1.8)))

# Ordenar por mejor puntuación
mejores_combinaciones.sort(key=lambda x: x[0])

print("TOP 5 COMBINACIONES (factores 1-3):")
for i, (punt, err, ext, params) in enumerate(mejores_combinaciones[:5]):
    print(f"{i+1}. Parámetros: {params}")
    print(f"   Puntuación: {punt:.4f}, Error: {err:.2e}, Éxito: {ext:.1%}")