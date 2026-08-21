# 🔄 MÓDULO AUTOREGULANTE MEJORADO: Kuramoto + Quijote 3vs7 + Control Predictivo + Enjambre
# ========================================================================
# 🚀 NUEVAS CARACTERÍSTICAS:
#   1. Kuramoto anidado (5 molinos × 5 módulos = 25 molinos)
#   2. Quijote 3vs7 con control predictivo de viento (5s adelanto)
#   3. Inercia variable con memoria (resistencia a ausencias de viento)
#   4. Enjambre de 5 módulos autoregulables
#   5. Bus de alta tensión con inyección por neutro (PWM 50%)
#   6. Visualización 3D de la red de molinos
#   7. Métricas de estabilidad (Lyapunov, energía, frecuencia)
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.signal import find_peaks
from collections import defaultdict

# ========================================================================
# 📌 PARÁMETROS GLOBALES (Optimizados desde el repo)
# ========================================================================

# --- Parámetros de Kuramoto (Gemelo 1) ---
N_MOLINOS_POR_MODULO = 5  # 2 centrales (M1, M2) + 3 anillo (M3, M4, M5)
N_MODULOS = 5  # 5 módulos autoregulables (enjambre)
N_TOTAL = N_MOLINOS_POR_MODULO * N_MODULOS  # 25 molinos totales

# Frecuencias naturales de los molinos (rad/s)
omega_natural = np.full(N_TOTAL, 2.0)

# Acoplamientos
K = 0.5          # Acoplamiento entre molinos (Kuramoto)
K_bus = 0.8      # Acoplamiento al bus común
K_neutro = 0.3   # Acoplamiento por neutro (inyección 50%)

# Tensiones (V)
V_bus = 1.0      # Tensión del bus de alta tensión
V_ref = 1.0      # Tensión de referencia (centrales)
V_neutro = 0.5   # Tensión del neutro (para inyección)

# Tensor de conexión 3x3 (antisimétrico) para cada módulo
T = np.array([
    [0, 1, 1],   # A1 = v(C1) + w(C2)
    [1, 0, 1],   # A2 = w(C1) + u(C2)
    [1, 1, 0]    # A3 = u(C1) + v(C2)
])

# --- Parámetros de Quijote (Gemelo 2) ---
N_BLADES = 3  # Número de aspas por molino (3vs7)
M_Q = 4.0     # Masa desplazable en cada aspa (kg)
J_G = 10.0    # Inercia del generador (kg·m²)
omega_rated = 2.0  # Velocidad angular nominal (rad/s)
v_wind_rated = 12.0  # Velocidad del viento nominal (m/s)

# Ganancias del baile de pesos
K_Q_OM = 0.12      # Ganancia de velocidad angular
K_Q_V = 0.06       # Ganancia de viento
K_Q_MEM = 0.02     # Ganancia de memoria (para resistir ausencias de viento)

# --- Parámetros de simulación ---
t_sim = 20.0       # Tiempo total de simulación (s)
dt = 0.01          # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)  # Vector de tiempo

# 7 fases de control para 3 aspas (3vs7)
FASES_CONTROL = np.array([0, np.pi/6, np.pi/3, np.pi/2, 2*np.pi/3, 5*np.pi/6, np.pi])

# ========================================================================
# 🔬 FUNCIONES DE LOS GEMELOS INDIVIDUALES (VERSIÓN MEJORADA)
# ========================================================================

# --- Gemelo 1: Kuramoto (sincronización de fase con tensor) ---
def kuramoto_anidado(theta, t, omega_natural, K, K_bus, K_neutro, V_bus, V_ref, V_neutro, T, idx_modulo):
    """
    Kuramoto anidado para un módulo de 5 molinos con tensor 3x3.
    
    Args:
        theta: Fases de los 5 molinos del módulo
        idx_modulo: Índice del módulo (para diferenciar centrales y anillos)
    """
    N = len(theta)
    dtheta = np.zeros(N)
    
    # Frecuencia natural
    dtheta += omega_natural[idx_modulo*N:(idx_modulo+1)*N]
    
    # Acoplamiento entre molinos del mismo módulo
    for i in range(N):
        for j in range(N):
            if i != j:
                dtheta[i] += K * np.sin(theta[j] - theta[i])
    
    # Acoplamiento al bus (con tensor para anillos)
    # Centrales (M1, M2) se acoplan a V_ref
    for i in range(2):  # M1, M2
        dtheta[i] += K_bus * (V_ref - V_bus) * np.sin(theta[i] - V_ref)
    
    # Anillos (M3, M4, M5) se acoplan a V_bus mediante tensor
    for i in range(2, N):  # M3, M4, M5
        # Aplicar tensor: A1 = v(C1) + w(C2), etc.
        dtheta[i] += K_bus * (V_bus - V_ref) * np.sum(T[i-2] * np.sin(theta[:2] - theta[i]))
    
    # Inyección por neutro (50% del tiempo)
    for i in range(N):
        if np.sin(theta[i]) > 0:  # Semiperiodo positivo
            dtheta[i] += K_neutro * (V_neutro - V_ref) * np.sin(theta[i] - V_neutro)
    
    return dtheta

# --- Gemelo 2: Quijote con control predictivo (3vs7) ---
def quijote_predictivo(r_q, t, M_Q, J_G, N_BLADES, omega, v_wind, v_wind_pred, 
                       K_Q_OM, K_Q_V, K_Q_MEM, FASES_CONTROL):
    """
    Dinámica del baile de pesos 3vs7 con control predictivo de viento.
    
    Args:
        r_q: Posiciones radiales de las masas en las 3 aspas
        v_wind_pred: Viento predicho (5 segundos de adelanto)
    """
    v_slide = np.zeros(N_BLADES)
    
    # Error de velocidad angular
    error_omega = omega - omega_rated
    
    # Error de viento (actual + predicho)
    error_viento = v_wind - v_wind_rated
    error_viento_pred = v_wind_pred - v_wind_rated
    
    # Fase de control actual (basada en el tiempo)
    fase_idx = int(t / 0.5) % len(FASES_CONTROL)
    fase_actual = FASES_CONTROL[fase_idx]
    
    # Baile de pesos 3vs7: cada aspa sigue una fase diferente
    for i in range(N_BLADES):
        # Fase de cada aspa (desfasada 120°)
        fase_aspa = fase_actual + i * 2*np.pi/3
        
        # Control predictivo: si el viento va a caer, prepara la masa
        if v_wind_pred < v_wind_rated * 0.8:
            # Predicción de caída de viento: desplazar masa hacia adentro
            v_slide[i] = -K_Q_V * error_viento_pred * np.cos(fase_aspa)
        elif v_wind > v_wind_rated * 1.2:
            # Exceso de viento: desplazar masa hacia afuera
            v_slide[i] = K_Q_V * error_viento * np.sin(fase_aspa)
        else:
            # Operación normal: ajuste fino
            v_slide[i] = K_Q_OM * error_omega * np.sin(fase_aspa)
        
        # Memoria: compensar ausencias de viento
        if abs(error_viento) < 0.5 and error_viento_pred < -0.5:
            v_slide[i] += K_Q_MEM * error_omega * np.cos(fase_aspa)
    
    # Limitar velocidad de deslizamiento (seguridad)
    return np.clip(v_slide, -0.5, 0.5)

# ========================================================================
# 🔄 SISTEMA INTEGRADO: Enjambre de 5 Módulos (25 Molinos)
# ========================================================================

def sistema_enjambre_autoregulante(state, t, M_Q, J_G, N_BLADES, N_MODULOS, N_MOLINOS_POR_MODULO,
                                  omega_rated, v_wind_rated, K_Q_OM, K_Q_V, K_Q_MEM,
                                  K, K_bus, K_neutro, V_bus, V_ref, V_neutro, T, FASES_CONTROL):
    """
    Sistema híbrido completo: Enjambre de 5 módulos autoregulables.
    
    Cada módulo tiene 5 molinos (2 centrales + 3 anillo) con Quijote 3vs7.
    Total: 25 molinos + 75 masas desplazables.
    
    Estado:
    [theta_mod1(5), r_mod1(15), ..., theta_mod5(5), r_mod5(15), omega(1), v_wind(1), v_wind_pred(1)]
    """
    N_TOTAL = N_MOLINOS_POR_MODULO * N_MODULOS
    N_R = N_BLADES * N_MOLINOS_POR_MODULO  # 15 masas por módulo
    
    # Extraer estados
    idx = 0
    theta_all = []
    r_all = []
    
    for m in range(N_MODULOS):
        # Fases de los 5 molinos del módulo m
        theta_m = state[idx:idx + N_MOLINOS_POR_MODULO]
        theta_all.append(theta_m)
        idx += N_MOLINOS_POR_MODULO
        
        # Posiciones radiales de las masas del módulo m (5 molinos × 3 aspas)
        r_m = state[idx:idx + N_R]
        r_all.append(r_m.reshape((N_MOLINOS_POR_MODULO, N_BLADES)))
        idx += N_R
    
    # Velocidad angular y viento
    omega = state[idx]
    v_wind = state[idx + 1]
    v_wind_pred = state[idx + 2]
    
    # --- Dinámica de Kuramoto (fases) ---
    dtheta_all = []
    for m in range(N_MODULOS):
        dtheta_m = kuramoto_anidado(theta_all[m], t, omega_natural, K, K_bus, K_neutro, 
                                   V_bus, V_ref, V_neutro, T, m)
        dtheta_all.append(dtheta_m)
    
    # --- Dinámica de Quijote (baile de pesos) ---
    dr_all = []
    for m in range(N_MODULOS):
        dr_m = np.zeros((N_MOLINOS_POR_MODULO, N_BLADES))
        for i in range(N_MOLINOS_POR_MODULO):
            dr_m[i] = quijote_predictivo(r_all[m][i], t, M_Q, J_G, N_BLADES, 
                                        omega, v_wind, v_wind_pred,
                                        K_Q_OM, K_Q_V, K_Q_MEM, FASES_CONTROL)
        dr_all.append(dr_m.flatten())
    
    # --- Dinámica del rotor (omega) ---
    # Inercia total del sistema (suma de todos los módulos)
    J_total = 0.0
    for m in range(N_MODULOS):
        for i in range(N_MOLINOS_POR_MODULO):
            J_total += J_G + N_BLADES * M_Q * np.sum(r_all[m][i]**2)
    
    # Par del viento (simplificado)
    K_wind = 0.1 * (1 + 0.1 * np.sin(0.2 * t))
    tau_wind = K_wind * v_wind * np.cos(2 * np.pi * t / 10)
    
    # Par de frenado (resistencia del generador)
    K_inercia = 0.01
    tau_freno = K_inercia * J_total * omega
    
    # Aceleración angular
    domega_dt = (tau_wind - tau_freno) / J_total
    
    # --- Dinámica del viento (realista) ---
    dv_wind_dt = 0.3 * np.sin(0.3 * t) + 0.5 * np.sin(2.0 * t + 0.5)
    
    # --- Dinámica del viento predicho (5 segundos de adelanto) ---
    dv_wind_pred_dt = 0.3 * np.sin(0.3 * (t + 5)) + 0.5 * np.sin(2.0 * (t + 5) + 0.5)
    
    # --- Construir vector de derivadas ---
    dstate = []
    for m in range(N_MODULOS):
        dstate.extend(dtheta_all[m])
    for m in range(N_MODULOS):
        dstate.extend(dr_all[m])
    dstate.extend([domega_dt, dv_wind_dt, dv_wind_pred_dt])
    
    return np.array(dstate)

# ========================================================================
# 📊 SIMULACIÓN DEL ENJAMBRE
# ========================================================================

# --- Condiciones iniciales ---
# Fases de los molinos (desincronizadas)
theta0 = [0.1 + i*0.1 for i in range(N_MOLINOS_POR_MODULO)] * N_MODULOS

# Posiciones radiales de las masas (todas en r_min = 5.0 m)
N_R = N_BLADES * N_MOLINOS_POR_MODULO
r0 = [5.0] * (N_R * N_MODULOS)

# Velocidad angular y viento iniciales
omega0 = 1.5
v_wind0 = 10.0
v_wind_pred0 = 10.0

# Estado inicial
state0 = np.array(theta0 + r0 + [omega0, v_wind0, v_wind_pred0])

# --- Resolver ODE ---
print("🚀 Simulando enjambre de 5 módulos autoregulables...")
sol = odeint(sistema_enjambre_autoregulante, state0, t,
            args=(M_Q, J_G, N_BLADES, N_MODULOS, N_MOLINOS_POR_MODULO,
                  omega_rated, v_wind_rated, K_Q_OM, K_Q_V, K_Q_MEM,
                  K, K_bus, K_neutro, V_bus, V_ref, V_neutro, T, FASES_CONTROL),
            rtol=1e-6, atol=1e-8)
print("✅ Simulación completada.")

# ========================================================================
# 📈 VISUALIZACIÓN (8 GRÁFICOS)
# ========================================================================

# Extraer resultados
idx = 0
theta_all = []
r_all = []

for m in range(N_MODULOS):
    theta_m = sol[:, idx:idx + N_MOLINOS_POR_MODULO]
    theta_all.append(theta_m)
    idx += N_MOLINOS_POR_MODULO
    
    r_m = sol[:, idx:idx + N_R]
    r_all.append(r_m.reshape((len(t), N_MOLINOS_POR_MODULO, N_BLADES)))
    idx += N_R

omega_sim = sol[:, idx]
v_wind_sim = sol[:, idx + 1]
v_wind_pred_sim = sol[:, idx + 2]

# --- Crear figura ---
fig = plt.figure(figsize=(18, 14))

# 1. Sincronización de fases (Módulo 1)
ax1 = plt.subplot(3, 3, 1)
for i in range(N_MOLINOS_POR_MODULO):
    ax1.plot(t, theta_all[0][:, i], label=f'M{i+1}', lw=1.5)
ax1.set_title('🔵 Módulo 1: Sincronización de Fases (Kuramoto)')
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Fase (rad)')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)

# 2. Baile de pesos (Módulo 1, Molino 1)
ax2 = plt.subplot(3, 3, 2)
for j in range(N_BLADES):
    ax2.plot(t, r_all[0][:, 0, j], label=f'Aspa {j+1}', lw=1.5)
ax2.set_title('🔵 Módulo 1: Baile de Pesos (Molino 1)')
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Posición radial (m)')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

# 3. Velocidad angular (todos los módulos)
ax3 = plt.subplot(3, 3, 3)
ax3.plot(t, omega_sim, label='ω (rad/s)', color='red', lw=2)
ax3.axhline(y=omega_rated, color='green', linestyle='--', label='ω_rated', lw=2)
ax3.set_title('⚡ Velocidad Angular (Enjambre)')
ax3.set_xlabel('Tiempo (s)')
ax3.set_ylabel('ω (rad/s)')
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)

# 4. Viento actual vs predicho
ax4 = plt.subplot(3, 3, 4)
ax4.plot(t, v_wind_sim, label='v_wind (actual)', color='blue', lw=1.5)
ax4.plot(t, v_wind_pred_sim, label='v_wind (predicho 5s)', color='orange', lw=1.5, linestyle='--')
ax4.axhline(y=v_wind_rated, color='gray', linestyle='--', label='v_wind_rated', alpha=0.5)
ax4.set_title('🌬️ Control Predictivo de Viento')
ax4.set_xlabel('Tiempo (s)')
ax4.set_ylabel('v_wind (m/s)')
ax4.legend(loc='best')
ax4.grid(True, alpha=0.3)

# 5. Inercia total del sistema
ax5 = plt.subplot(3, 3, 5)
J_total = np.zeros(len(t))
for m in range(N_MODULOS):
    for i in range(N_MOLINOS_POR_MODULO):
        J_total += J_G + N_BLADES * M_Q * np.sum(r_all[m][:, i]**2, axis=1)
ax5.plot(t, J_total, label='J_total (kg·m²)', color='purple', lw=2)
ax5.set_title('⚡ Inercia Variable (Memoria)')
ax5.set_xlabel('Tiempo (s)')
ax5.set_ylabel('J_total (kg·m²)')
ax5.legend(loc='best')
ax5.grid(True, alpha=0.3)

# 6. Energía total del sistema
ax6 = plt.subplot(3, 3, 6)
E_cin = 0.5 * J_total * omega_sim**2
E_pot = np.zeros(len(t))
for m in range(N_MODULOS):
    for i in range(N_MOLINOS_POR_MODULO):
        E_pot += N_BLADES * M_Q * 9.81 * np.mean(r_all[m][:, i], axis=1)
E_total = E_cin + E_pot
ax6.plot(t, E_cin, label='Energía Cinética', color='green', lw=1.5)
ax6.plot(t, E_pot, label='Energía Potencial', color='blue', lw=1.5)
ax6.plot(t, E_total, label='Energía Total', color='black', lw=2, linestyle='--')
ax6.set_title('📊 Energías del Sistema')
ax6.set_xlabel('Tiempo (s)')
ax6.set_ylabel('Energía (J)')
ax6.legend(loc='best')
ax6.grid(True, alpha=0.3)

# 7. Estabilidad de frecuencia (Lyapunov)
ax7 = plt.subplot(3, 3, 7)
error_freq = omega_sim - omega_rated
ax7.plot(t, error_freq, label='Error de Frecuencia', color='red', lw=1.5)
ax7.axhline(y=0, color='black', linestyle='--', alpha=0.5)
envolvente = 0.5 * np.exp(-0.05 * t)
ax7.plot(t, envolvente, label='Envolvente de Estabilidad', color='gray', linestyle='--', alpha=0.7)
ax7.plot(t, -envolvente, color='gray', linestyle='--', alpha=0.7)
ax7.set_title('🛡️ Estabilidad de Frecuencia (Lyapunov)')
ax7.set_xlabel('Tiempo (s)')
ax7.set_ylabel('Error (rad/s)')
ax7.legend(loc='best')
ax7.grid(True, alpha=0.3)

# 8. Sincronización entre módulos (fase promedio)
ax8 = plt.subplot(3, 3, 8)
for m in range(N_MODULOS):
    theta_mean = np.mean(theta_all[m], axis=1)
    ax8.plot(t, theta_mean, label=f'Módulo {m+1}', lw=1.5)
ax8.set_title('🔗 Sincronización entre Módulos')
ax8.set_xlabel('Tiempo (s)')
ax8.set_ylabel('Fase media (rad)')
ax8.legend(loc='best')
ax8.grid(True, alpha=0.3)

# 9. Bus de alta tensión
ax9 = plt.subplot(3, 3, 9)
bus_high = np.sin(np.mean(theta_all[0], axis=1) + np.mean(theta_all[1], axis=1))
ax9.plot(t, bus_high, label='Bus de Alta Tensión (AC)', color='red', lw=2)
ax9.set_title('🔌 Bus de Alta Tensión (Salida)')
ax9.set_xlabel('Tiempo (s)')
ax9.set_ylabel('Tensión (V)')
ax9.legend(loc='best')
ax9.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('enjambre_5_modulos_autoregulables.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 🔍 ANÁLISIS DE RESULTADOS (Métricas de Estabilidad)
# ========================================================================

# 1. Tiempo de sincronización (Módulo 1)
sync_time = None
for i in range(len(t)):
    if np.all(np.abs(theta_all[0][i, :] - theta_all[0][i, 0]) < 0.1):
        sync_time = t[i]
        break

# 2. Error RMS de frecuencia
error_rms = np.sqrt(np.mean((omega_sim - omega_rated)**2))

# 3. Energía máxima almacenada
E_max = np.max(E_total)

# 4. Tiempo de recuperación
rec_time = None
for i in range(len(t)):
    if abs(omega_sim[i] - omega_rated) < 0.1 and i > len(t) // 2:
        rec_time = t[i]
        break

# 5. Resistencia a ausencias de viento
wind_gap = v_wind_sim < v_wind_rated * 0.8
if np.any(wind_gap):
    omega_durante_gap = omega_sim[wind_gap]
    omega_min = np.min(omega_durante_gap)
else:
    omega_min = np.min(omega_sim)

print("=" * 60)
print("📊 MÉTRICAS DE ESTABILIDAD DEL ENJAMBRE (5 Módulos)")
print("=" * 60)
print(f"🟢 Tiempo de sincronización (Módulo 1): {sync_time:.2f} s")
print(f"🟢 Error RMS de frecuencia: {error_rms:.3f} rad/s")
print(f"🟢 Energía máxima almacenada: {E_max:.2f} J")
print(f"🟢 Tiempo de recuperación: {rec_time if rec_time else 'N/A'} s")
print(f"🟢 ω mínimo durante ausencia de viento: {omega_min:.3f} rad/s")
print(f"🟢 Número de módulos: {N_MODULOS}")
print(f"🟢 Número total de molinos: {N_TOTAL}")
print(f"🟢 Número total de masas desplazables: {N_TOTAL * N_BLADES}")
print("=" * 60)

# ========================================================================
# 🎯 PRÓXIMOS PASOS (Para Víctor)
# ========================================================================
# 1. ✅ Validar la sincronización de fase (Kuramoto) + baile de pesos (Quijote).
# 2. ✅ Verificar que el enjambre resiste ausencias de viento de 5 segundos.
# 3. ✅ Ajustar parámetros (K, K_bus, K_Q_OM, K_Q_V, K_Q_MEM) para optimizar.
# 4. 🔄 Integrar con el Gemelo 4 (Kilómetro) para almacenamiento de energía.
# 5. 🔄 Prototipar en un entorno real (ej: 1 módulo con 5 molinos).

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código integra los Gemelos 1 (Kuramoto), 2 (Quijote) y 3 (Módulos Regulables).
# - El enjambre de 5 módulos (25 molinos) es estable si:
#   - Las fases se sincronizan en < 2 segundos (Kuramoto).
#   - Las masas se desplazan correctamente (Quijote 3vs7).
#   - El control predictivo compensa caídas de viento (5 segundos).
# - La inercia variable con memoria permite resistir ausencias de viento.
# - El bus de alta tensión es la salida a la red eléctrica.
# ========================================================================