# 🔄 MÓDULO AUTOREGULANTE: Integración Kuramoto + Quijote (3vs7)
# ========================================================================
# Descripción: Simulación del sistema híbrido que integra:
#             - Gemelo 1: Red de 5 molinos (Kuramoto + Tensor 3x3)
#             - Gemelo 2: Quijote con baile de pesos (3vs7)
# Autor: Lee Arbusto & Víctor M. A. (espiradesombra)
# Fecha: 21-08-2026
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# ========================================================================
# 📌 PARÁMETROS GLOBALES (Combinación de Gemelo 1 y Gemelo 2)
# ========================================================================

# --- Parámetros de Kuramoto (Gemelo 1) ---
N_molinos = 5  # 2 centrales (M1, M2) + 3 anillo (M3, M4, M5)
omega_natural = np.array([2.0, 2.0, 2.0, 2.0, 2.0])  # Frecuencia natural (rad/s)
K = 0.5  # Acoplamiento entre molinos
K_bus = 0.8  # Acoplamiento al bus
V_bus = 1.0  # Tensión del bus
V_ref = 1.0  # Tensión de referencia

# Tensor de conexión 3x3 (antisimétrico)
T = np.array([
    [0, 1, 1],  # A1 = v(C1) + w(C2)
    [1, 0, 1],  # A2 = w(C1) + u(C2)
    [1, 1, 0]   # A3 = u(C1) + v(C2)
])

# --- Parámetros de Quijote (Gemelo 2) ---
N_blades = 3  # Número de aspas por molino
M_Q = 4.0  # Masa desplazable (kg)
J_G = 10.0  # Inercia del generador (kg·m²)
omega_rated = 2.0  # Velocidad angular nominal (rad/s)
v_wind_rated = 12.0  # Velocidad del viento nominal (m/s)
K_Q_OM = 0.1  # Ganancia de velocidad angular
K_Q_V = 0.05  # Ganancia de viento

# --- Parámetros de simulación ---
t_sim = 20.0  # Tiempo total (s)
dt = 0.01  # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)  # Vector de tiempo

# ========================================================================
# 🔬 FUNCIONES DE LOS GEMELOS INDIVIDUALES
# ========================================================================

# --- Gemelo 1: Kuramoto (sincronización de fase) ---
def kuramoto(theta, t, omega, K, K_bus, V_bus, V_ref, T):
    """Ecuaciones de Kuramoto para 5 molinos con tensor de conexión."""
    N = len(theta)
    dtheta = np.zeros(N)
    dtheta += omega  # Frecuencia natural
    
    # Acoplamiento entre molinos
    for i in range(N):
        for j in range(N):
            if i != j:
                dtheta[i] += K * np.sin(theta[j] - theta[i])
    
    # Acoplamiento al bus (con tensor para anillos)
    for i in range(2, N):  # Solo anillos (M3, M4, M5)
        dtheta[i] += K_bus * (V_bus - V_ref) * np.sum(T[i-2] * np.sin(theta[:2] - theta[i]))
    
    return dtheta

# --- Gemelo 2: Quijote (baile de pesos) ---
def quijote_dynamics(r_q, t, M_Q, J_G, N_blades, omega, v_wind, K_Q_OM, K_Q_V):
    """Dinámica del baile de pesos en las aspas."""
    v_slide = np.zeros(N_blades)
    for i in range(N_blades):
        v_slide[i] = K_Q_OM * (omega - omega_rated) + K_Q_V * (v_wind - v_wind_rated)
    return v_slide

# ========================================================================
# 🔄 SISTEMA INTEGRADO: Kuramoto + Quijote
# ========================================================================

def sistema_autoregulante(state, t, M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V, K, K_bus, V_bus, V_ref, T):
    """
    Sistema híbrido: 5 molinos (Kuramoto) + Quijote (3vs7) en cada molino.
    
    Estado: [theta1, theta2, ..., theta5, r1_1, r1_2, r1_3, ..., r5_1, r5_2, r5_3, omega, v_wind]
    - theta1..theta5: Fases de los 5 molinos (Kuramoto).
    - r1_1..r5_3: Posiciones radiales de las masas en cada aspa de cada molino (Quijote).
    - omega: Velocidad angular del rotor (común para todos los molinos).
    - v_wind: Velocidad del viento.
    """
    # Extraer estados
    theta = state[:N_molinos]  # Fases de los molinos
    r_q = state[N_molinos:N_molinos + N_molinos * N_blades]  # Posiciones radiales (5 molinos × 3 aspas)
    omega = state[N_molinos + N_molinos * N_blades]  # Velocidad angular
    v_wind = state[N_molinos + N_molinos * N_blades + 1]  # Velocidad del viento
    
    # Reshape r_q a una matriz (5 molinos × 3 aspas)
    r_q_matrix = r_q.reshape((N_molinos, N_blades))
    
    # --- Dinámica de Kuramoto (fases) ---
    dtheta_dt = kuramoto(theta, t, omega_natural, K, K_bus, V_bus, V_ref, T)
    
    # --- Dinámica de Quijote (baile de pesos) ---
    dr_q_dt = np.zeros(N_molinos * N_blades)
    for i in range(N_molinos):
        # Para cada molino, calcular el baile de pesos en sus 3 aspas
        start_idx = i * N_blades
        dr_q_dt[start_idx:start_idx + N_blades] = quijote_dynamics(
            r_q_matrix[i], t, M_Q, J_G, N_blades, omega, v_wind, K_Q_OM, K_Q_V
        )
    
    # --- Dinámica del rotor (omega) ---
    # Inercia total del sistema (suma de inercias de todos los molinos)
    J_total = 0.0
    for i in range(N_molinos):
        J_total += J_G + N_blades * M_Q * np.sum(r_q_matrix[i]**2)
    
    # Dinámica de omega (simplificada: depende de v_wind y J_total)
    K_wind = 0.1
    K_inercia = 0.01
    domega_dt = K_wind * v_wind - K_inercia * J_total * omega
    
    # --- Dinámica del viento (oscilante para pruebas) ---
    dv_wind_dt = 0.5 * np.sin(0.5 * t)
    
    # Devolver derivadas en el mismo orden que el estado
    return np.concatenate([dtheta_dt, dr_q_dt, [domega_dt], [dv_wind_dt]])

# ========================================================================
# 📊 SIMULACIÓN DEL SISTEMA INTEGRADO
# ========================================================================

# --- Condiciones iniciales ---
# Fases de los molinos (desincronizadas)
theta0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5])

# Posiciones radiales de las masas (todas en r_min inicialmente)
r_q0 = np.full(N_molinos * N_blades, 5.0)  # 5 molinos × 3 aspas = 15 masas

# Velocidad angular y viento iniciales
omega0 = 1.5
v_wind0 = 10.0

# Estado inicial: [theta1..theta5, r1_1..r5_3, omega, v_wind]
state0 = np.concatenate([theta0, r_q0, [omega0], [v_wind0]])

# --- Resolver ODE ---
sol = odeint(sistema_autoregulante, state0, t, 
            args=(M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V, K, K_bus, V_bus, V_ref, T))

# Extraer resultados
theta_sim = sol[:, :N_molinos]  # Fases de los molinos
r_q_sim = sol[:, N_molinos:N_molinos + N_molinos * N_blades]  # Posiciones radiales
omega_sim = sol[:, N_molinos + N_molinos * N_blades]  # Velocidad angular
v_wind_sim = sol[:, N_molinos + N_molinos * N_blades + 1]  # Velocidad del viento

# Reshape r_q_sim a (tiempo × 5 molinos × 3 aspas)
r_q_sim_matrix = r_q_sim.reshape((len(t), N_molinos, N_blades))

# ========================================================================
# 📈 VISUALIZACIÓN
# ========================================================================

plt.figure(figsize=(18, 12))

# --- Gráfico 1: Fases de los molinos (Kuramoto) ---
plt.subplot(4, 1, 1)
for i in range(N_molinos):
    plt.plot(t, theta_sim[:, i], label=f'Molino {i+1}')
plt.title('Sincronización de Fases (Kuramoto)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Fase (rad)')
plt.legend()
plt.grid(True)

# --- Gráfico 2: Posiciones radiales (Quijote) ---
plt.subplot(4, 1, 2)
for i in range(N_molinos):
    for j in range(N_blades):
        plt.plot(t, r_q_sim_matrix[:, i, j], label=f'Molino {i+1}, Aspa {j+1}')
plt.title('Baile de Pesos (Quijote 3vs7)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Posición radial (m)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)

# --- Gráfico 3: Velocidad angular del rotor ---
plt.subplot(4, 1, 3)
plt.plot(t, omega_sim, label='ω (rad/s)', color='red')
plt.axhline(y=omega_rated, color='green', linestyle='--', label='ω_rated')
plt.title('Velocidad Angular del Rotor')
plt.xlabel('Tiempo (s)')
plt.ylabel('ω (rad/s)')
plt.legend()
plt.grid(True)

# --- Gráfico 4: Velocidad del viento ---
plt.subplot(4, 1, 4)
plt.plot(t, v_wind_sim, label='v_wind (m/s)', color='blue')
plt.axhline(y=v_wind_rated, color='orange', linestyle='--', label='v_wind_rated')
plt.title('Velocidad del Viento (Oscilante)')
plt.xlabel('Tiempo (s)')
plt.ylabel('v_wind (m/s)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('modulo_autoregulante_kuramoto_quijote.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 🔍 ANÁLISIS DE RESULTADOS
# ========================================================================

# --- Tiempo de sincronización (Kuramoto) ---
sync_time = None
for i in range(len(t)):
    if np.all(np.abs(theta_sim[i, :] - theta_sim[i, 0]) < 0.1):
        sync_time = t[i]
        break

print(f"🔹 Tiempo de sincronización (Kuramoto): {sync_time:.2f} segundos")

# --- Energía cinética total ---
J_total = np.zeros(len(t))
for i in range(len(t)):
    for j in range(N_molinos):
        J_total[i] += J_G + N_blades * M_Q * np.sum(r_q_sim_matrix[i, j]**2)
E_cin = 0.5 * J_total * omega_sim**2

# --- Energía potencial total ---
g = 9.81
E_pot = np.zeros(len(t))
for i in range(len(t)):
    for j in range(N_molinos):
        E_pot[i] += N_blades * M_Q * g * np.mean(r_q_sim_matrix[i, j])

# --- Energía total ---
E_total = E_cin + E_pot

print(f"🔹 Energía cinética final: {E_cin[-1]:.2f} J")
print(f"🔹 Energía potencial final: {E_pot[-1]:.2f} J")
print(f"🔹 Energía total final: {E_total[-1]:.2f} J")
print(f"🔹 Velocidad angular final: {omega_sim[-1]:.2f} rad/s")

# --- Posiciones radiales finales ---
print(f"🔹 Posiciones radiales finales (Molino 1): {r_q_sim_matrix[-1, 0]}")

# ========================================================================
# 🎯 VISUALIZACIÓN 3D: Baile de Pesos en un Molino
# ========================================================================

# Crear figura 3D para el Molino 1
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Datos para el Molino 1 (3 aspas)
molino_1_r = r_q_sim_matrix[:, 0, :]  # Posiciones radiales de las 3 aspas del Molino 1

# Coordenadas polares a cartesianas (para visualizar el baile)
theta_aspas = np.linspace(0, 2*np.pi, 100)  # Ángulo de las aspas (fijo)
for i in range(N_blades):
    x = molino_1_r[:, i] * np.cos(theta_aspas[i])
    y = molino_1_r[:, i] * np.sin(theta_aspas[i])
    z = np.zeros_like(t)  # Altura constante (simplificación)
    ax.plot(x, y, z, label=f'Aspa {i+1}')

ax.set_title('Baile de Pesos en Molino 1 (3D)')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.legend()
plt.savefig('baile_pesos_3d_molino1.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 📌 PRÓXIMOS PASOS (Para Víctor)
# ========================================================================
# 1. Validar la sincronización de fases (Kuramoto) + baile de pesos (Quijote).
# 2. Ajustar parámetros (K, K_bus, K_Q_OM, K_Q_V) para optimizar el sistema.
# 3. Prototipar en un entorno real (ej: 1 molino con Quijote + sincronización Kuramoto).
# 4. Integrar con el Gemelo 4 (Kilómetro) para almacenamiento de energía.

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código integra los Gemelos 1 (Kuramoto) y 2 (Quijote) en un módulo autoregulante.
# - El sistema es estable si:
#   - Las fases se sincronizan en < 2 segundos (Kuramoto).
#   - Las masas se desplazan correctamente para resistir ausencias de viento (Quijote).
# - Para validar experimentalmente, necesitarás datos reales de molinos.
# ========================================================================