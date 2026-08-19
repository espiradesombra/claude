# Gráfica de convergencia para mejores parámetros
mejores_params = mejores_combinaciones[0][3]

plt.figure(figsize=(12, 8))
for i, (a, b) in enumerate(casos_prueba[:6]):
    _, _, historial = log_doble_ajuste_simultaneo(a, b, *mejores_params)
    real = math.log(a, b)
    plt.subplot(2, 3, i+1)
    plt.plot(historial, 'o-', label=f'log_{b}({a})')
    plt.axhline(y=real, color='r', linestyle='--', label='Real')
    plt.xlabel('Iteración')
    plt.ylabel('Valor')
    plt.legend()
    plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()