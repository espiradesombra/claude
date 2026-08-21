# 🌐 GEMelo 1: Red de 5 Molinos Autoregulable (Kuramoto Anidado + Tensor 3x3)
# ========================================================================
# Descripción: Simulación de sincronización de fase en una red de 5 molinos
#             (2 centrales + 3 anillo) con tensor de conexión antisimétrico.
# Autor: Lee Arbusto & Víctor M. A. (espiradesombra)
# Fecha: 21-08-2026
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# ========================================================================
# 📌 PARÁMETROS DEL SISTEMA (Óptimos según tus repos)
# ========================================================================

# Frecuencias naturales de los molinos (rad/s)
omega = np.array([2.0, 2.0, 2.0, 2.0, 2.0])  # M1, M2 (centrales), M3, M4, M5 (anillo)

# Acoplamientos
K = 0.5       # Acoplamiento entre molinos (Kuramoto)
K_bus = 0.8   # Acoplamiento al bus común

# Tensiones (V)
V_bus = 1.0   # Tensión del bus de alta tensión
V_ref = 1.0   # Tensión de referencia (centrales)

# Tensor de conexión 3x3 (antisimétrico)
# Representa cómo se acoplan las fases entre centrales (C1, C2) y anillos (A1, A2, A3)
# u, v, w = fases desfasadas 120° (u=0°, v=120°, w=240°)
T = np.array([
    [0, 1, 1],   # A1 = v(C1) + w(C2)
    [1, 0, 1],   # A2 = w(C1) + u(C2)
    [1, 1, 0]    # A3 = u(C1) + v(C2)
])

# ========================================================================
# 🔬 MODELO MATEMÁTICO (Ecuaciones de Kuramoto)
# ========================================================================

def kuramoto(theta, t, omega, K, K_bus, V_bus, V_ref, T):
    """
    Ecuaciones de Kuramoto para 5 molinos con tensor de conexión.
    
    Args:
        theta: Array de fases de los molinos [θ1, θ2, θ3, θ4, θ5]
        t: Tiempo
        omega: Frecuencias naturales
        K: Acoplamiento entre molinos
        K_bus: Acoplamiento al bus
        V_bus: Tensión del bus
        V_ref: Tensión de referencia
        T: Tensor de conexión (3x3)
    
    Returns:
        dtheta/dt: Derivadas de las fases
    """
    N = len(theta)
    dtheta = np.zeros(N)
    
    # Término de frecuencia natural
    dtheta += omega
    
    # Término de acoplamiento entre molinos (Kuramoto clásico)
    for i in range(N):
        for j in range(N):
            if i != j:
                dtheta[i] += K * np.sin(theta[j] - theta[i])
    
    # Término de acoplamiento al bus (con tensor)
    # Los centrales (M1, M2) marcan V_ref; los anillos (M3-M5) se acoplan a V_bus
    for i in range(2, N):  # Solo para anillos (M3, M4, M5)
        # Aplicar tensor: A1 = v(C1) + w(C2), etc.
        # Simplificación: Usamos V_bus - V_ref como diferencia de potencial
        dtheta[i] += K_bus * (V_bus - V_ref) * np.sum(T[i-2] * np.sin(theta[:2] - theta[i]))
    
    return dtheta

# ========================================================================
# 📊 SIMULACIÓN
# ========================================================================

# Condiciones iniciales: fases desincronizadas
theta0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5])  # rad

# Tiempo de simulación (0 a 10 segundos)
t = np.linspace(0, 10, 1000)

# Resolver ODE
sol = odeint(kuramoto, theta0, t, args=(omega, K, K_bus, V_bus, V_ref, T))

# ========================================================================
# 📈 VISUALIZACIÓN
# ========================================================================

plt.figure(figsize=(15, 10))

# --- Gráfico 1: Evolución de las fases ---
plt.subplot(3, 1, 1)
for i in range(5):
    plt.plot(t, sol[:, i], label=f'Molino {i+1}')
plt.title('Evolución de las fases de los molinos (Kuramoto + Tensor)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Fase (rad)')
plt.legend()
plt.grid(True)

# --- Gráfico 2: Diferencias de fase (sincronización) ---
plt.subplot(3, 1, 2)
for i in range(1, 5):
    plt.plot(t, sol[:, i] - sol[:, 0], label=f'Δθ M{i+1}-M1')
plt.title('Diferencias de fase respecto a M1 (Central)')
plt.xlabel('Tiempo (s)')
plt.ylabel('ΔFase (rad)')
plt.legend()
plt.grid(True)

# --- Gráfico 3: Tensor de conexión (representación 3D) ---
plt.subplot(3, 1, 3, projection='3d')
X, Y = np.meshgrid([0, 1, 2], [0, 1, 2])
Z = T
fig = plt.gca()
fig.plot_surface(X, Y, Z, cmap='viridis')
plt.title('Tensor de Conexión 3x3 (Antisimétrico)')
plt.xticks([0, 1, 2], ['C1', 'C2', ' '])
plt.yticks([0, 1, 2], ['C1', 'C2', ' '])
plt.zticks([0, 1], ['0', '1'])

plt.tight_layout()
plt.savefig('kuramoto_5molinos_tensor.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 🔍 ANÁLISIS DE RESULTADOS
# ========================================================================

# Tiempo de sincronización (cuando todas las fases están a < 0.1 rad de diferencia)
sync_time = None
for i in range(len(t)):
    if np.all(np.abs(sol[i, :] - sol[i, 0]) < 0.1):
        sync_time = t[i]
        break

print(f"🔹 Tiempo de sincronización: {sync_time:.2f} segundos")
print(f"🔹 Frecuencia final (M1): {sol[-1, 0]:.2f} rad")
print(f"🔹 Diferencias finales de fase: {sol[-1, :] - sol[-1, 0]}")

# ========================================================================
# 📌 PRÓXIMOS PASOS (Para Víctor)
# ========================================================================
# 1. Ajustar K y K_bus para optimizar el tiempo de sincronización.
# 2. Validar el tensor T con tus datos reales de molinos.
# 3. Integrar con el Gemelo 2 (Quijote) para inercia variable.
# 4. Conectar a un bus real de alta tensión (simulación en PSS®E o DIgSILENT).

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código simula la sincronización de 5 molinos con un tensor 3x3.
# - El tensor T define cómo se acoplan las fases entre centrales y anillos.
# - Los parámetros K y K_bus son críticos para la estabilidad.
# - Para validar experimentalmente, necesitarás datos reales de molinos.
# ========================================================================