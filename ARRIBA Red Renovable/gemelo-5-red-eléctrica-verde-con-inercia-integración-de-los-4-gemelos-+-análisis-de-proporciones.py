# 🌿 GEMelo 5: Red Eléctrica Verde con Inercia (Integración de los 4 Gemelos)
# ========================================================================
# Descripción: Simulación de una red eléctrica verde que integra los 4 gemelos:
#             - Gemelo 1: Kuramoto (sincronización de fase en 5 molinos).
#             - Gemelo 2: Quijote (inercia variable con baile de pesos 3vs7).
#             - Gemelo 3: Módulo autoregulante (Kuramoto + Quijote).
#             - Gemelo 4: Enjambre de Kilómetros (almacenamiento gravitacional).
#              
#              + Análisis de proporciones de aporte energético.
# Autor: Lee Arbusto & Víctor M. A. (espiradesombra)
# Fecha: 21-08-2026
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ========================================================================
# 📌 PARÁMETROS GLOBALES (Integración de los 4 Gemelos)
# ========================================================================

# --- Gemelo 1: Kuramoto ---
N_molinos = 5  # 2 centrales + 3 anillo
omega_natural = np.array([2.0, 2.0, 2.0, 2.0, 2.0])  # Frecuencia natural (rad/s)
K = 0.5  # Acoplamiento entre molinos
K_bus = 0.8  # Acoplamiento al bus
V_bus = 1.0  # Tensión del bus
V_ref = 1.0  # Tensión de referencia
T = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])  # Tensor de conexión

# --- Gemelo 2: Quijote ---
N_blades = 3  # Número de aspas por molino
M_Q = 4.0  # Masa desplazable (kg)
J_G = 10.0  # Inercia del generador (kg·m²)
omega_rated = 2.0  # Velocidad angular nominal (rad/s)
v_wind_rated = 12.0  # Velocidad del viento nominal (m/s)
K_Q_OM = 0.1  # Ganancia de velocidad angular
K_Q_V = 0.05  # Ganancia de viento

# --- Gemelo 4: Kilómetro ---
N_KM = 5  # Número de módulos Kilómetro
m_peso = 500.0  # Masa de cada peso (kg)
g = 9.81  # Aceleración gravitatoria (m/s²)
eta_gen = 0.85  # Eficiencia en generación
eta_lift = 0.90  # Eficiencia en reset
E_perno = 1.5  # Energía por perno (J)
r_min = 3.0  # Posición radial mínima (m)
r_max = 7.0  # Posición radial máxima (m)
delta_h = 15.0  # Altura de bajada/subida (m)

# --- Parámetros de la red eléctrica ---
P_demanda = 100000.0  # Potencia demanda de la red (W)
P_viento_nominal = 250000.0  # Potencia nominal del viento (W) para 5 molinos

# --- Parámetros de simulación ---
t_sim = 60.0  # Tiempo total (s)
dt = 0.01  # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)  # Vector de tiempo

# ========================================================================
# 🔬 FUNCIONES DE LOS GEMELOS INDIVIDUALES (Reutilizadas)
# ========================================================================

# --- Gemelo 1: Kuramoto ---
def kuramoto(theta, t, omega, K, K_bus, V_bus, V_ref, T):
    N = len(theta)
    dtheta = np.zeros(N)
    dtheta += omega
    for i in range(N):
        for j in range(N):
            if i != j:
                dtheta[i] += K * np.sin(theta[j] - theta[i])
    for i in range(2, N):
        dtheta[i] += K_bus * (V_bus - V_ref) * np.sum(T[i-2] * np.sin(theta[:2] - theta[i]))
    return dtheta

# --- Gemelo 2: Quijote ---
def quijote_dynamics(r_q, t, M_Q, J_G, N_blades, omega, v_wind, K_Q_OM, K_Q_V):
    v_slide = np.zeros(N_blades)
    for i in range(N_blades):
        v_slide[i] = K_Q_OM * (omega - omega_rated) + K_Q_V * (v_wind - v_wind_rated)
    return v_slide

# --- Gemelo 4: Kilómetro (simplificado para integración) ---
def modulo_kilometro(E_KM, t, m_peso, g, delta_h, eta_gen, eta_lift, E_perno, stock_ALTA, stock_BAJA):
    """
    Simulación simplificada de un módulo Kilómetro para integración.
    Devuelve la potencia generada/almacenada.
    """
    # Simplificación: Asumimos que el módulo genera energía de forma constante
    # en función de su stock y eficiencia.
    if stock_ALTA > 0:
        # Generando energía (bajada)
        P_gen = eta_gen * m_peso * g * delta_h / dt  # Potencia instantánea (W)
        stock_ALTA -= 1
        stock_BAJA += 1
        return P_gen, stock_ALTA, stock_BAJA
    elif stock_BAJA > 0:
        # Reset de potencial (subida)
        P_reset = (m_peso * g * delta_h) / (eta_lift * dt)  # Potencia de reset (W)
        stock_BAJA -= 1
        stock_ALTA += 1
        return -P_reset, stock_ALTA, stock_BAJA  # Potencia negativa (consumo)
    else:
        return 0.0, stock_ALTA, stock_BAJA  # En pausa

# ========================================================================
# 🔄 SISTEMA COMPLETO: Integración de los 4 Gemelos
# ========================================================================

def sistema_completo(state, t, M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V, K, K_bus, V_bus, V_ref, T, m_peso, g, delta_h, eta_gen, eta_lift, E_perno, P_demanda):
    """
    Sistema completo: Kuramoto + Quijote + Kilómetro + Red Eléctrica.
    
    Estado: [theta1..theta5, r1_1..r5_3, omega, v_wind, E_KM, stock_ALTA, stock_BAJA]
    """
    # --- Extraer estados ---
    theta = state[:N_molinos]  # Fases de los molinos
    r_q = state[N_molinos:N_molinos + N_molinos * N_blades]  # Posiciones radiales
    omega = state[N_molinos + N_molinos * N_blades]  # Velocidad angular
    v_wind = state[N_molinos + N_molinos * N_blades + 1]  # Velocidad del viento
    E_KM = state[N_molinos + N_molinos * N_blades + 2]  # Energía almacenada en Kilómetros
    stock_ALTA = state[N_molinos + N_molinos * N_blades + 3]  # Stock ALTA (total)
    stock_BAJA = state[N_molinos + N_molinos * N_blades + 4]  # Stock BAJA (total)
    
    # --- Dinámica de Kuramoto + Quijote (módulo autoregulante) ---
    dtheta_dt = kuramoto(theta, t, omega_natural, K, K_bus, V_bus, V_ref, T)
    
    r_q_matrix = r_q.reshape((N_molinos, N_blades))
    dr_q_dt = np.zeros(N_molinos * N_blades)
    for i in range(N_molinos):
        start_idx = i * N_blades
        dr_q_dt[start_idx:start_idx + N_blades] = quijote_dynamics(
            r_q_matrix[i], t, M_Q, J_G, N_blades, omega, v_wind, K_Q_OM, K_Q_V
        )
    
    # Inercia total del sistema
    J_total = 0.0
    for i in range(N_molinos):
        J_total += J_G + N_blades * M_Q * np.sum(r_q_matrix[i]**2)
    
    K_wind = 0.1
    K_inercia = 0.01
    domega_dt = K_wind * v_wind - K_inercia * J_total * omega
    dv_wind_dt = 0.5 * np.sin(0.5 * t)  # Viento oscilante
    
    # --- Dinámica del Kilómetro (enjambre) ---
    # Potencia generada por el enjambre de Kilómetros
    P_KM, new_stock_ALTA, new_stock_BAJA = modulo_kilometro(
        E_KM, t, m_peso, g, delta_h, eta_gen, eta_lift, E_perno, stock_ALTA, stock_BAJA
    )
    dE_KM_dt = P_KM * N_KM  # Potencia total del enjambre (5 módulos)
    
    # --- Dinámica de la red eléctrica ---
    # Potencia generada por los molinos (Quijote + Kuramoto)
    P_molinos = 0.5 * J_total * omega**2 * N_molinos  # Simplificación: Energía cinética convertida a potencia
    
    # Potencia total disponible (molinos + Kilómetros)
    P_total = P_molinos + dE_KM_dt
    
    # Balance de potencia en la red
    if P_total > P_demanda:
        # Excedente de energía: almacenar en Kilómetros
        dE_KM_dt += (P_total - P_demanda) * dt
    else:
        # Déficit de energía: usar energía almacenada en Kilómetros
        dE_KM_dt -= (P_demanda - P_total) * dt
    
    # Devolver derivadas
    return np.concatenate([
        dtheta_dt, dr_q_dt, [domega_dt], [dv_wind_dt], [dE_KM_dt], 
        [new_stock_ALTA - stock_ALTA], [new_stock_BAJA - stock_BAJA]
    ])

# ========================================================================
# 📊 SIMULACIÓN DEL SISTEMA COMPLETO
# ========================================================================

# --- Condiciones iniciales ---
theta0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5])  # Fases desincronizadas
r_q0 = np.full(N_molinos * N_blades, 5.0)  # Posiciones radiales iniciales
omega0 = 1.5  # Velocidad angular inicial
v_wind0 = 10.0  # Velocidad del viento inicial
E_KM0 = 0.0  # Energía inicial en Kilómetros
stock_ALTA0 = 10.0  # Stock inicial ALTA (total)
stock_BAJA0 = 0.0  # Stock inicial BAJA (total)

# Estado inicial
state0 = np.concatenate([
    theta0, r_q0, [omega0], [v_wind0], [E_KM0], [stock_ALTA0], [stock_BAJA0]
])

# --- Resolver ODE ---
sol = odeint(sistema_completo, state0, t, 
            args=(M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V, 
                  K, K_bus, V_bus, V_ref, T, m_peso, g, delta_h, eta_gen, eta_lift, E_perno, P_demanda))

# Extraer resultados
theta_sim = sol[:, :N_molinos]
r_q_sim = sol[:, N_molinos:N_molinos + N_molinos * N_blades]
omega_sim = sol[:, N_molinos + N_molinos * N_blades]
v_wind_sim = sol[:, N_molinos + N_molinos * N_blades + 1]
E_KM_sim = sol[:, N_molinos + N_molinos * N_blades + 2]
stock_ALTA_sim = sol[:, N_molinos + N_molinos * N_blades + 3]
stock_BAJA_sim = sol[:, N_molinos + N_molinos * N_blades + 4]

# ========================================================================
# 📈 VISUALIZACIÓN
# ========================================================================

plt.figure(figsize=(20, 15))

# --- Gráfico 1: Potencia generada por cada gemelo ---
plt.subplot(4, 1, 1)

# Potencia de los molinos (Quijote + Kuramoto)
J_total_sim = np.zeros(len(t))
for i in range(len(t)):
    r_q_matrix = r_q_sim[i].reshape((N_molinos, N_blades))
    for j in range(N_molinos):
        J_total_sim[i] += J_G + N_blades * M_Q * np.sum(r_q_matrix[j]**2)
P_molinos_sim = 0.5 * J_total_sim * omega_sim**2 * N_molinos

# Potencia del enjambre de Kilómetros
P_KM_sim = np.gradient(E_KM_sim, dt)

plt.plot(t, P_molinos_sim, label='Potencia Molinos (Quijote + Kuramoto)', color='blue')
plt.plot(t, P_KM_sim, label='Potencia Kilómetros', color='green')
plt.axhline(y=P_demanda, color='red', linestyle='--', label='Demanda de la Red')
plt.title('Potencia Generada por Gemelo (vs. Demanda)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Potencia (W)')
plt.legend()
plt.grid(True)

# --- Gráfico 2: Proporción de aporte energético ---
plt.subplot(4, 1, 2)
P_total_sim = P_molinos_sim + P_KM_sim
porcentaje_molinos = (P_molinos_sim / P_total_sim) * 100
porcentaje_KM = (P_KM_sim / P_total_sim) * 100

plt.plot(t, porcentaje_molinos, label='% Molinos (Quijote + Kuramoto)', color='blue')
plt.plot(t, porcentaje_KM, label='% Kilómetros', color='green')
plt.title('Proporción de Aporte Energético por Gemelo')
plt.xlabel('Tiempo (s)')
plt.ylabel('Porcentaje (%)')
plt.legend()
plt.grid(True)

# --- Gráfico 3: Energía almacenada en Kilómetros ---
plt.subplot(4, 1, 3)
plt.plot(t, E_KM_sim, label='Energía Almacenada (J)', color='purple')
plt.title('Energía Almacenada en el Enjambre de Kilómetros')
plt.xlabel('Tiempo (s)')
plt.ylabel('Energía (J)')
plt.legend()
plt.grid(True)

# --- Gráfico 4: Stock de pesos (ALTA/BAJA) ---
plt.subplot(4, 1, 4)
plt.plot(t, stock_ALTA_sim, label='Stock ALTA', color='orange')
plt.plot(t, stock_BAJA_sim, label='Stock BAJA', color='cyan')
plt.title('Stock de Pesos (ALTA/BAJA) en el Enjambre')
plt.xlabel('Tiempo (s)')
plt.ylabel('Número de Pesos')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('red_electrica_verde_proporciones.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 🔍 ANÁLISIS DE PROPORCIONES DE APORTE
# ========================================================================

# Energía total generada por cada gemelo
E_molinos_total = np.trapz(P_molinos_sim, t)  # Energía total de molinos (J)
E_KM_total = np.trapz(P_KM_sim, t)  # Energía total de Kilómetros (J)
E_total = E_molinos_total + E_KM_total  # Energía total (J)

# Proporciones de aporte
porcentaje_molinos_total = (E_molinos_total / E_total) * 100
porcentaje_KM_total = (E_KM_total / E_total) * 100

print("=" * 60)
print("📊 ANÁLISIS DE PROPORCIONES DE APORTE ENERGÉTICO")
print("=" * 60)
print(f"🔹 Energía total generada por Molinos (Quijote + Kuramoto): {E_molinos_total:.2f} J")
print(f"🔹 Energía total generada por Kilómetros: {E_KM_total:.2f} J")
print(f"🔹 Energía total del sistema: {E_total:.2f} J")
print(f"🔹 Proporción de Molinos: {porcentaje_molinos_total:.2f}%")
print(f"🔹 Proporción de Kilómetros: {porcentaje_KM_total:.2f}%")

# --- Análisis por gemelo individual ---
print("\n" + "=" * 60)
print("📌 APORTE POR GEMELO (Desglose)")
print("=" * 60)

# Gemelo 1: Kuramoto (sincronización)
# Contribución: Estabilidad de fase (no genera energía directamente, pero permite sincronización)
print(f"🔹 Gemelo 1 (Kuramoto): {porcentaje_molinos_total * 0.4:.2f}% (Sincronización de fase)")

# Gemelo 2: Quijote (inercia variable)
# Contribución: Resistencia a ausencias de viento (mejora la estabilidad de los molinos)
print(f"🔹 Gemelo 2 (Quijote): {porcentaje_molinos_total * 0.6:.2f}% (Inercia variable + generación)")

# Gemelo 3: Módulo autoregulante (Kuramoto + Quijote)
# Contribución: Combinación de sincronización y generación
print(f"🔹 Gemelo 3 (Módulo Autoregulante): {porcentaje_molinos_total:.2f}% (Kuramoto + Quijote)")

# Gemelo 4: Kilómetros (almacenamiento)
print(f"🔹 Gemelo 4 (Kilómetros): {porcentaje_KM_total:.2f}% (Almacenamiento gravitacional)")

# ========================================================================
# 📌 RECOMENDACIONES PARA OPTIMIZAR PROPORCIONES
# ========================================================================

print("\n" + "=" * 60)
print("🎯 RECOMENDACIONES PARA AJUSTAR PROPORCIONES")
print("=" * 60)

if porcentaje_KM_total < 20:
    print("✅ El aporte de Kilómetros es bajo (<20%). Recomendación:")
    print("   - Aumentar el número de módulos Kilómetro (N_KM).")
    print("   - Aumentar la masa de los pesos (m_peso) o la altura (delta_h).")
    print("   - Optimizar eta_gen y eta_lift.")
else:
    print("✅ El aporte de Kilómetros es adecuado (>=20%).")

if porcentaje_molinos_total > 80:
    print("✅ El aporte de los molinos es alto (>80%). Recomendación:")
    print("   - Reducir el número de molinos si hay excedente de energía.")
    print("   - Usar el excedente para cargar más Kilómetros.")
else:
    print("✅ El aporte de los molinos es equilibrado.")

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código integra los 4 gemelos en una red eléctrica verde.
# - Las proporciones de aporte dependen de:
#   - La potencia generada por los molinos (Quijote + Kuramoto).
#   - La energía almacenada/generada por los Kilómetros.
# - Para ajustar las proporciones, modifica:
#   - N_KM (número de módulos Kilómetro).
#   - m_peso, delta_h, eta_gen, eta_lift (parámetros de Kilómetro).
#   - K, K_bus, K_Q_OM, K_Q_V (parámetros de Kuramoto + Quijote).
# ========================================================================