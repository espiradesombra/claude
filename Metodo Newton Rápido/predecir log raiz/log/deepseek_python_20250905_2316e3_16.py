# (El código de optimización que te envié anteriormente se está ejecutando)

# AL FINALIZAR, TE ENTREGARÉ:
print("\n" + "="*60)
print("🎯 RESULTADOS COMPLETOS DE OPTIMIZACIÓN")
print("="*60)

# 1. TOP 10 MEJORES COMBINACIONES
print("\n🔝 TOP 10 COMBINACIONES:")
for i, (score, f1, f2, f4) in enumerate(mejores_resultados[:10]):
    print(f"{i+1:2d}. f1={f1:.1f}, f2={f2:.1f}, f4={f4:.1f} → Score: {score:.6f}")

# 2. MEJORES PARÁMETROS ENCONTRADOS
mejor_score, mejor_f1, mejor_f2, mejor_f4 = mejores_resultados[0]
print(f"\n🏆 MEJOR COMBINACIÓN ENCONTRADA:")
print(f"   • factor1 (sobrepaso): {mejor_f1}")
print(f"   • factor2 (multiplicativo): {mejor_f2}") 
print(f"   • factor4 (exponencial): {mejor_f4}")
print(f"   • Puntuación: {mejor_score:.6f}")

# 3. COMPARATIVA CON PARÁMETROS ORIGINALES
print(f"\n📊 COMPARACIÓN CON PARÁMETROS ORIGINALES:")
print(f"   Original: f1=1.3, f2=3.5, f4=1.8 → Score: {evaluacion_rapida(1.3, 3.5, 1.8):.6f}")
print(f"   Optimizado: f1={mejor_f1}, f2={mejor_f2}, f4={mejor_f4} → Score: {mejor_score:.6f}")

# 4. PRUEBA RÁPIDA DE MEJORA
print(f"\n🚀 PRUEBA DE MEJORA:")
test_cases = [(100, 10), (8, 2), (1024, 2), (math.pi, 10)]
for a, b in test_cases:
    original = evaluacion_rapida(1.3, 3.5, 1.8) 
    optimizado = mejor_score
    mejora = ((original - optimizado) / original) * 100
    print(f"   log_{b}({a}): {mejora:+.1f}% de mejora")

# 5. RECOMENDACIONES PARA SIGUIENTE FASE
print(f"\n🎯 SIGUIENTES PASOS RECOMENDADOS:")
print("   1. Validar con más casos de prueba")
print("   2. Optimizar factores 3, 5, 6 (secundarios)")
print("   3. Pruebas de estrés con valores extremos")
print("   4. Análisis de estabilidad numérica")

print("\n" + "="*60)
print("¡OPTIMIZACIÓN COMPLETADA! 🎉")
print("="*60)