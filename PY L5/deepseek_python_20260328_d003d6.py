def graficar_historial(hist):
    if not hist:
        return
    iters = list(range(len(hist)))
    ms = [h[0] for h in hist]
    fases = [h[1] for h in hist]
    plt.figure(figsize=(12,4))
    plt.subplot(1,2,1)
    plt.plot(iters, ms, 'o-')
    plt.xlabel('Iteración')
    plt.ylabel('m')
    plt.title('Evolución de m')
    plt.grid(True)
    plt.subplot(1,2,2)
    plt.plot(iters, fases, 'o-', color='orange')
    plt.xlabel('Iteración')
    plt.ylabel('Fase mínima')
    plt.title('Evolución de la fase')
    plt.grid(True)
    plt.tight_layout()
    plt.show()