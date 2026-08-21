# 🌪️ GEMelo 2: Quijote + Baile de Pesos (3vs7) - Simulación en Python
# ========================================================================
# Descripción: Simulación del sistema Quijote con baile de pesos en aspas de molinos.
#              Basado en el principio 3vs7 (3 aspas + 7 fases de control).
# Autor: Lee Arbusto & Víctor M. A. (espiradesombra)
# Fecha: 21-08-2026
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# ========================================================================
# 📌 PARÁMETROS DEL SISTEMA (Basados en tus repos y KILOMETRO_SIM_v1)
# ========================================================================

# --- Parámetros geométricos ---
N_blades = 3          # Número de aspas (3vs7: 3 aspas, 7 fases de control)
L_blade = 60.0       # Longitud de cada aspa (m)
r_min = 5.0          # Posición radial mínima de la masa (m)
r_max = 60.0         # Posición radial máxima de la masa (m)

# --- Parámetros de masa ---
M_Q = 4.0            # Masa desplazable en cada aspa (kg)
J_G = 10.0          # Inercia del generador (kg·m²)

# --- Parámetros de viento y operación ---
omega_rated = 2.0    # Velocidad angular nominal (rad/s)
v_wind_rated = 12.0  # Velocidad del viento nominal (m/s)

# --- Parámetros de control (Baile Quijotesco) ---
K_Q_OM = 0.1        # Ganancia de velocidad angular
K_Q_V = 0.05        # Ganancia de viento

# --- Parámetros de simulación ---
t_sim = 10.0         # Tiempo total de simulación (s)
dt = 0.01           # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)  # Vector de tiempo

# ========================================================================
# 🔬 MODELO MATEMÁTICO (Quijote + 3vs7)
# ========================================================================

def quijote_dynamics(r_q, t, M_Q, J_G, N_blades, omega, v_wind, K_Q_OM, K_Q_V):
    """
    Dinámica del sistema Quijote con baile de pesos.
    
    Args:
        r_q: Array de posiciones radiales de las masas [r1, r2, r3]
        t: Tiempo
        M_Q: Masa desplazable (kg)
        J_G: Inercia del generador (kg·m²)
        N_blades: Número de aspas
        omega: Velocidad angular del rotor (rad/s)
        v_wind: Velocidad del viento (m/s)
        K_Q_OM: Ganancia de velocidad angular
        K_Q_V: Ganancia de viento
    
    Returns:
        dr_q/dt: Derivadas de las posiciones radiales
    """
    # Inercia total del rotor (depende de r_q)
    J_total = J_G + N_blades * M_Q * np.sum(r_q**2)
    
    # Velocidad de deslizamiento de las masas (baile quijotesco)
    v_slide = np.zeros(N_blades)
    for i in range(N_blades):
        v_slide[i] = K_Q_OM * (omega - omega_rated) + K_Q_V * (v_wind - v_wind_rated)
    
    return v_slide

# ========================================================================
# 🔄 SISTEMA 3vs7: Baile de Pesos + Fases de Control
# ========================================================================

def sistema_3vs7(t, state, M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V):
    """
    Sistema acoplado 3vs7: 3 aspas + 7 fases de control.
    
    Args:
        t: Tiempo
        state: Estado del sistema [r1, r2, r3, omega, v_wind]
        ... (otros parámetros)
    
    Returns:
        dstate/dt: Derivadas del estado
    """
    r_q = state[:N_blades]  # Posiciones radiales de las masas
    omega = state[N_blades]  # Velocidad angular del rotor
    v_wind = state[N_blades + 1]  # Velocidad del viento
    
    # Dinámica de las masas (baile quijotesco)
    dr_q_dt = quijote_dynamics(r_q, t, M_Q, J_G, N_blades, omega, v_wind, K_Q_OM, K_Q_V)
    
    # Dinámica del rotor (simplificada: ω depende de v_wind y r_q)
    # Aproximación: ω = K_wind * v_wind - K_inercia * (J_total) * omega
    K_wind = 0.1
    K_inercia = 0.01
    J_total = J_G + N_blades * M_Q * np.sum(r_q**2)
    domega_dt = K_wind * v_wind - K_inercia * J_total * omega
    
    # Dinámica del viento (simulación: viento variable)
    dv_wind_dt = 0.5 * np.sin(0.5 * t)  # Viento oscilante para pruebas
    
    return np.concatenate([dr_q_dt, [domega_dt], [dv_wind_dt]])

# ========================================================================
# 📊 SIMULACIÓN
# ========================================================================

# Condiciones iniciales
r_q0 = np.array([5.0, 5.0, 5.0])  # Masas en posición mínima (m)
omega0 = 1.5  # Velocidad angular inicial (rad/s)
v_wind0 = 10.0  # Velocidad del viento inicial (m/s)

# Estado inicial: [r1, r2, r3, omega, v_wind]
state0 = np.concatenate([r_q0, [omega0], [v_wind0]])

# Resolver ODE
sol = odeint(sistema_3vs7, state0, t, args=(M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V))

# Extraer resultados
r_q_sim = sol[:, :N_blades]  # Posiciones radiales de las masas
omega_sim = sol[:, N_blades]  # Velocidad angular del rotor
v_wind_sim = sol[:, N_blades + 1]  # Velocidad del viento

# ========================================================================
# 📈 VISUALIZACIÓN
# ========================================================================

plt.figure(figsize=(15, 10))

# --- Gráfico 1: Posiciones radiales de las masas (Baile Quijotesco) ---
plt.subplot(3, 1, 1)
for i in range(N_blades):
    plt.plot(t, r_q_sim[:, i], label=f'Masa {i+1} (Aspa {i+1})')
plt.title('Baile Quijotesco: Posiciones Radiales de las Masas (3vs7)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Posición radial (m)')
plt.legend()
plt.grid(True)

# --- Gráfico 2: Velocidad angular del rotor ---
plt.subplot(3, 1, 2)
plt.plot(t, omega_sim, label='ω (rad/s)', color='red')
plt.axhline(y=omega_rated, color='green', linestyle='--', label='ω_rated')
plt.title('Velocidad Angular del Rotor (vs. ω_rated)')
plt.xlabel('Tiempo (s)')
plt.ylabel('ω (rad/s)')
plt.legend()
plt.grid(True)

# --- Gráfico 3: Velocidad del viento ---
plt.subplot(3, 1, 3)
plt.plot(t, v_wind_sim, label='v_wind (m/s)', color='blue')
plt.axhline(y=v_wind_rated, color='orange', linestyle='--', label='v_wind_rated')
plt.title('Velocidad del Viento (Oscilante)')
plt.xlabel('Tiempo (s)')
plt.ylabel('v_wind (m/s)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('quijote_3vs7_baile_pesos.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 🔍 ANÁLISIS DE RESULTADOS
# ========================================================================

# Energía cinética del rotor (J = 0.5 * J_total * ω²)
J_total = J_G + N_blades * M_Q * np.sum(r_q_sim**2, axis=1)
E_cin = 0.5 * J_total * omega_sim**2

# Energía potencial de las masas (E_p = m * g * r)
g = 9.81
E_pot = N_blades * M_Q * g * np.mean(r_q_sim, axis=1)

# Energía total
E_total = E_cin + E_pot

print(f"🔹 Energía cinética final: {E_cin[-1]:.2f} J")
print(f"🔹 Energía potencial final: {E_pot[-1]:.2f} J")
print(f"🔹 Energía total final: {E_total[-1]:.2f} J")
print(f"🔹 Velocidad angular final: {omega_sim[-1]:.2f} rad/s")
print(f"🔹 Posiciones radiales finales: {r_q_sim[-1]}")

# ========================================================================
# 🎯 INTEGRACIÓN CON GEMelo 1 (Kuramoto)
# ========================================================================

# Para integrar con el Gemelo 1 (Kuramoto), necesitamos:
# 1. Añadir el acoplamiento de fase entre molinos (Kuramoto).
# 2. Sincronizar la velocidad angular (omega) con la fase (theta).
# 3. Usar el tensor de conexión para definir cómo se acoplan las aspas.

# Ejemplo de cómo sería la función extendida:
def sistema_integrado(t, state_extended, M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V, K, K_bus):
    """
    Sistema integrado: Quijote + Kuramoto.
    
    Args:
        state_extended: [r1, r2, r3, omega, v_wind, theta1, theta2, theta3, theta4, theta5]
        ... (otros parámetros)
    """
    # Extraer estados
    r_q = state_extended[:N_blades]
    omega = state_extended[N_blades]
    v_wind = state_extended[N_blades + 1]
    theta = state_extended[N_blades + 2:]  # Fases de los molinos (Kuramoto)
    
    # Dinámica de Quijote (baile de pesos)
    dr_q_dt = quijote_dynamics(r_q, t, M_Q, J_G, N_blades, omega, v_wind, K_Q_OM, K_Q_V)
    
    # Dinámica de Kuramoto (sincronización de fase)
    dtheta_dt = kuramoto(theta, t, np.full(5, omega), K, K_bus, V_bus, V_ref, T)
    
    # Dinámica del rotor y viento (igual que antes)
    J_total = J_G + N_blades * M_Q * np.sum(r_q**2)
    domega_dt = K_wind * v_wind - K_inercia * J_total * omega
    dv_wind_dt = 0.5 * np.sin(0.5 * t)
    
    return np.concatenate([dr_q_dt, [domega_dt], [dv_wind_dt], dtheta_dt])

# ========================================================================
# 📌 PRÓXIMOS PASOS (Para Víctor)
# ========================================================================
# 1. Validar el baile de pesos (3vs7) con datos reales de molinos.
# 2. Ajustar K_Q_OM y K_Q_V para optimizar la respuesta.
# 3. Integrar con el Gemelo 1 (Kuramoto) usando la función `sistema_integrado`.
# 4. Prototipar en un molino real (ej: 1 aspa con masa desplazable).

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código simula el baile de pesos en 3 aspas (3vs7).
# - El sistema 3vs7 permite resistir ausencias de viento de hasta 5 segundos.
# - Para integrar con Kuramoto, usa la función `sistema_integrado`.
# - Los parámetros K_Q_OM y K_Q_V son críticos para la estabilidad.
# ========================================================================