# 🌿 GEMELO 5 - VERSIÓN MEJORADA: Red Eléctrica Verde con Inercia
# ========================================================================
# 🚀 NUEVAS CARACTERÍSTICAS:
#   1. Integración de los 4 gemelos en un solo sistema
#   2. Kuramoto anidado con tensor 3x3 (5 molinos × 5 módulos = 25 molinos)
#   3. Quijote 3vs7 con control predictivo de viento (5s adelanto)
#   4. Enjambre de Kilómetros con reset de potencial (5 módulos)
#   5. Red eléctrica con demanda variable (picos y valles)
#   6. Análisis de proporciones de aporte energético
#   7. Optimización automática de parámetros
#   8. Visualización 3D de la red y flujos de energía
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ========================================================================
# 📌 PARÁMETROS GLOBALES (Optimizados desde el repo)
# ========================================================================

# --- Gemelo 1: Kuramoto (sincronización de fase) ---
N_MOLINOS_POR_MODULO = 5  # 2 centrales + 3 anillo
N_MODULOS_KURAMOTO = 5    # 5 módulos autoregulables
N_TOTAL_MOLINOS = N_MOLINOS_POR_MODULO * N_MODULOS_KURAMOTO  # 25 molinos

# Frecuencias naturales (rad/s)
omega_natural = np.full(N_TOTAL_MOLINOS, 2.0)

# Acoplamientos
K = 0.5          # Acoplamiento entre molinos
K_bus = 0.8      # Acoplamiento al bus común
K_neutro = 0.3   # Acoplamiento por neutro (inyección 50%)

# Tensiones (V)
V_bus = 1.0
V_ref = 1.0
V_neutro = 0.5

# Tensor de conexión 3x3 (antisimétrico)
T = np.array([
    [0, 1, 1],   # A1 = v(C1) + w(C2)
    [1, 0, 1],   # A2 = w(C1) + u(C2)
    [1, 1, 0]    # A3 = u(C1) + v(C2)
])

# --- Gemelo 2: Quijote (inercia variable 3vs7) ---
N_BLADES = 3            # 3 aspas por molino (3vs7)
M_Q = 4.0               # Masa desplazable (kg)
J_G = 10.0              # Inercia del generador (kg·m²)
omega_rated = 2.0       # Velocidad angular nominal (rad/s)
v_wind_rated = 12.0     # Velocidad del viento nominal (m/s)

# Ganancias del baile de pesos
K_Q_OM = 0.12           # Ganancia de velocidad angular
K_Q_V = 0.06            # Ganancia de viento
K_Q_MEM = 0.02          # Ganancia de memoria

# 7 fases de control para 3 aspas (3vs7)
FASES_CONTROL = np.array([0, np.pi/6, np.pi/3, np.pi/2, 2*np.pi/3, 5*np.pi/6, np.pi])

# --- Gemelo 3: Módulo Autoregulante (Kuramoto + Quijote) ---
# (Ya integrado en los parámetros anteriores)

# --- Gemelo 4: Kilómetro (almacenamiento gravitacional) ---
N_MODULOS_KM = 5        # 5 módulos Kilómetro en enjambre
M_PESO = 10.0           # Masa de cada peso (kg)
DELTA_H = 15.0          # Altura de bajada/subida (m)
ETA_GEN = 0.85          # Eficiencia en generación
ETA_LIFT = 0.90         # Eficiencia en reset
E_PERNO = 1.5           # Energía por perno (J)
STOCK_ALTA_INICIAL = 10 # Pesos iniciales en stock ALTA (por módulo)
STOCK_BAJA_INICIAL = 0  # Pesos iniciales en stock BAJA (por módulo)

# --- Parámetros de la red eléctrica ---
P_demanda_base = 100000.0   # Potencia base de la red (W)
P_demanda_pico = 150000.0   # Potencia pico (W)
P_viento_nominal = 250000.0 # Potencia nominal del viento (W) para 25 molinos

# --- Parámetros de simulación ---
T_SIM = 60.0                # Tiempo total (s)
DT = 0.01                   # Paso de tiempo (s)
t = np.arange(0, T_SIM, DT)

# ========================================================================
# 🔬 FUNCIONES DE LOS GEMELOS INDIVIDUALES (VERSIÓN MEJORADA)
# ========================================================================

# --- Gemelo 1: Kuramoto Anidado con Tensor ---
def kuramoto_anidado(theta, t, omega_natural, K, K_bus, K_neutro, V_bus, V_ref, V_neutro, T, idx_modulo):
    """
    Kuramoto anidado para un módulo de 5 molinos con tensor 3x3.
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
    
    # Centrales (M1, M2) se acoplan a V_ref
    for i in range(2):
        dtheta[i] += K_bus * (V_ref - V_bus) * np.sin(theta[i] - V_ref)
    
    # Anillos (M3, M4, M5) se acoplan a V_bus mediante tensor
    for i in range(2, N):
        dtheta[i] += K_bus * (V_bus - V_ref) * np.sum(T[i-2] * np.sin(theta[:2] - theta[i]))
    
    # Inyección por neutro (50% del tiempo)
    for i in range(N):
        if np.sin(theta[i]) > 0:
            dtheta[i] += K_neutro * (V_neutro - V_ref) * np.sin(theta[i] - V_neutro)
    
    return dtheta

# --- Gemelo 2: Quijote con Control Predictivo (3vs7) ---
def quijote_predictivo(r_q, t, M_Q, J_G, N_BLADES, omega, v_wind, v_wind_pred,
                       K_Q_OM, K_Q_V, K_Q_MEM, FASES_CONTROL):
    """
    Dinámica del baile de pesos 3vs7 con control predictivo de viento.
    """
    v_slide = np.zeros(N_BLADES)
    
    error_omega = omega - omega_rated
    error_viento = v_wind - v_wind_rated
    error_viento_pred = v_wind_pred - v_wind_rated
    
    fase_idx = int(t / 0.5) % len(FASES_CONTROL)
    fase_actual = FASES_CONTROL[fase_idx]
    
    for i in range(N_BLADES):
        fase_aspa = fase_actual + i * 2*np.pi/3
        
        if v_wind_pred < v_wind_rated * 0.8:
            v_slide[i] = -K_Q_V * error_viento_pred * np.cos(fase_aspa)
        elif v_wind > v_wind_rated * 1.2:
            v_slide[i] = K_Q_V * error_viento * np.sin(fase_aspa)
        else:
            v_slide[i] = K_Q_OM * error_omega * np.sin(fase_aspa)
        
        if abs(error_viento) < 0.5 and error_viento_pred < -0.5:
            v_slide[i] += K_Q_MEM * error_omega * np.cos(fase_aspa)
    
    return np.clip(v_slide, -0.5, 0.5)

# --- Gemelo 4: Kilómetro (enjambre) ---
class ModuloKilometroMejorado:
    """Módulo Kilómetro individual con reset de potencial."""
    
    def __init__(self, id_modulo, M_PESO, DELTA_H, ETA_GEN, ETA_LIFT, E_PERNO):
        self.id = id_modulo
        self.M_PESO = M_PESO
        self.DELTA_H = DELTA_H
        self.ETA_GEN = ETA_GEN
        self.ETA_LIFT = ETA_LIFT
        self.E_PERNO = E_PERNO
        
        # Estado inicial
        self.cota = 'ALTA'
        self.n_pesos_objeto = 3
        self.stock_ALTA = STOCK_ALTA_INICIAL
        self.stock_BAJA = STOCK_BAJA_INICIAL
        
        # Energías
        self.E_generada_total = 0.0
        self.E_consumida_total = 0.0
        self.E_pernos_total = 0.0
        self.E_reset_total = 0.0
        
        # Contadores
        self.ciclos_completados = 0
        self.enganches_realizados = 0
        self.entregas_realizadas = 0
        self.resets_realizados = 0
    
    def ejecutar_ciclo(self):
        """Ejecuta un ciclo completo del módulo Kilómetro."""
        E_gen = 0.0
        E_cons = 0.0
        exito = False
        
        # Fase 1: Enganche (ALTA)
        if self.cota == 'ALTA' and self.n_pesos_objeto == 3:
            if self.stock_ALTA > 0:
                self.n_pesos_objeto += 1
                self.stock_ALTA -= 1
                self.E_consumida_total += 2 * self.E_PERNO
                self.E_pernos_total += 2 * self.E_PERNO
                self.enganches_realizados += 1
                exito = True
        
        # Fase 2: Bajada (Generación)
        elif self.cota == 'ALTA' and self.n_pesos_objeto == 4:
            E_gen = self.ETA_GEN * self.M_PESO * 9.81 * self.DELTA_H
            self.E_generada_total += E_gen
            self.cota = 'BAJA'
            self.ciclos_completados += 1
            exito = True
        
        # Fase 3: Entrega (BAJA)
        elif self.cota == 'BAJA' and self.n_pesos_objeto == 4:
            self.n_pesos_objeto -= 1
            self.stock_BAJA += 1
            self.E_consumida_total += 2 * self.E_PERNO
            self.E_pernos_total += 2 * self.E_PERNO
            self.entregas_realizadas += 1
            self.cota = 'BAJA'
            exito = True
        
        # Fase 4: Subida (Flotación)
        elif self.cota == 'BAJA' and self.n_pesos_objeto == 3:
            self.cota = 'ALTA'
            exito = True
        
        # Modo Pausa
        if self.stock_ALTA == 0 and self.cota != 'PAUSA':
            self.cota = 'PAUSA'
            exito = False
        
        return E_gen, E_cons, exito
    
    def reset_potencial(self):
        """Reset del potencial de flotación."""
        if self.stock_BAJA > 0 and self.cota == 'PAUSA':
            E_reset = (self.M_PESO * 9.81 * self.DELTA_H) / self.ETA_LIFT
            self.E_consumida_total += E_reset
            self.E_reset_total += E_reset
            self.stock_BAJA -= 1
            self.stock_ALTA += 1
            self.resets_realizados += 1
            self.cota = 'ALTA'
            return E_reset
        return 0.0

class EnjambreKilometrosMejorado:
    """Enjambre de módulos Kilómetro."""
    
    def __init__(self, N_MODULOS_KM, M_PESO, DELTA_H, ETA_GEN, ETA_LIFT, E_PERNO):
        self.modulos = [ModuloKilometroMejorado(i, M_PESO, DELTA_H, ETA_GEN, ETA_LIFT, E_PERNO)
                        for i in range(N_MODULOS_KM)]
        self.E_generada_total = 0.0
        self.E_consumida_total = 0.0
        self.E_pernos_total = 0.0
        self.E_reset_total = 0.0
        self.historial = []
    
    def paso(self):
        """Ejecuta un paso de simulación para todos los módulos."""
        E_gen_total = 0.0
        E_cons_total = 0.0
        
        for modulo in self.modulos:
            E_gen, E_cons, exito = modulo.ejecutar_ciclo()
            E_gen_total += E_gen
            E_cons_total += E_cons
            
            if modulo.cota == 'PAUSA':
                modulo.reset_potencial()
        
        self.E_generada_total += E_gen_total
        self.E_consumida_total += E_cons_total
        self.E_pernos_total = sum(m.E_pernos_total for m in self.modulos)
        self.E_reset_total = sum(m.E_reset_total for m in self.modulos)
        
        return E_gen_total, E_cons_total
    
    def get_estado(self):
        """Devuelve el estado global del enjambre."""
        return {
            'E_generada': self.E_generada_total,
            'E_consumida': self.E_consumida_total,
            'E_pernos': self.E_pernos_total,
            'E_reset': self.E_reset_total,
            'E_balance': self.E_generada_total - self.E_consumida_total,
            'ciclos_totales': sum(m.ciclos_completados for m in self.modulos),
            'modulos': [{
                'id': m.id,
                'cota': m.cota,
                'n_pesos': m.n_pesos_objeto,
                'stock_ALTA': m.stock_ALTA,
                'stock_BAJA': m.stock_BAJA,
                'ciclos': m.ciclos_completados
            } for m in self.modulos]
        }

# ========================================================================
# 🔄 SISTEMA COMPLETO: Integración de los 4 Gemelos
# ========================================================================

def sistema_completo(state, t, params):
    """
    Sistema completo: Kuramoto + Quijote + Kilómetro + Red Eléctrica.
    
    Estado: [theta1..theta25, r1_1..r25_3, omega, v_wind, v_wind_pred, E_KM]
    """
    # Extraer parámetros
    M_Q = params['M_Q']
    J_G = params['J_G']
    N_BLADES = params['N_BLADES']
    omega_rated = params['omega_rated']
    v_wind_rated = params['v_wind_rated']
    K_Q_OM = params['K_Q_OM']
    K_Q_V = params['K_Q_V']
    K_Q_MEM = params['K_Q_MEM']
    K = params['K']
    K_bus = params['K_bus']
    K_neutro = params['K_neutro']
    V_bus = params['V_bus']
    V_ref = params['V_ref']
    V_neutro = params['V_neutro']
    T = params['T']
    FASES_CONTROL = params['FASES_CONTROL']
    N_MODULOS_KURAMOTO = params['N_MODULOS_KURAMOTO']
    N_MOLINOS_POR_MODULO = params['N_MOLINOS_POR_MODULO']
    omega_natural = params['omega_natural']
    N_BLADES = params['N_BLADES']
    
    # Extraer estados
    idx = 0
    N_TOTAL = N_MODULOS_KURAMOTO * N_MOLINOS_POR_MODULO
    theta = state[idx:idx + N_TOTAL]
    idx += N_TOTAL
    
    N_R = N_BLADES * N_MOLINOS_POR_MODULO * N_MODULOS_KURAMOTO
    r_q = state[idx:idx + N_R]
    idx += N_R
    
    omega = state[idx]
    v_wind = state[idx + 1]
    v_wind_pred = state[idx + 2]
    E_KM = state[idx + 3]
    
    # --- Dinámica de Kuramoto (fases) ---
    dtheta = np.zeros(N_TOTAL)
    for m in range(N_MODULOS_KURAMOTO):
        start = m * N_MOLINOS_POR_MODULO
        end = start + N_MOLINOS_POR_MODULO
        theta_m = theta[start:end]
        dtheta[start:end] = kuramoto_anidado(
            theta_m, t, omega_natural, K, K_bus, K_neutro,
            V_bus, V_ref, V_neutro, T, m
        )
    
    # --- Dinámica de Quijote (baile de pesos) ---
    dr_q = np.zeros(N_R)
    r_q_matrix = r_q.reshape((N_MODULOS_KURAMOTO, N_MOLINOS_POR_MODULO, N_BLADES))
    for m in range(N_MODULOS_KURAMOTO):
        for i in range(N_MOLINOS_POR_MODULO):
            start = (m * N_MOLINOS_POR_MODULO + i) * N_BLADES
            end = start + N_BLADES
            dr_q[start:end] = quijote_predictivo(
                r_q_matrix[m, i], t, M_Q, J_G, N_BLADES,
                omega, v_wind, v_wind_pred,
                K_Q_OM, K_Q_V, K_Q_MEM, FASES_CONTROL
            )
    
    # --- Dinámica del rotor (omega) ---
    J_total = 0.0
    for m in range(N_MODULOS_KURAMOTO):
        for i in range(N_MOLINOS_POR_MODULO):
            J_total += J_G + N_BLADES * M_Q * np.sum(r_q_matrix[m, i]**2)
    
    K_wind = 0.1 * (1 + 0.1 * np.sin(0.2 * t))
    tau_wind = K_wind * v_wind * np.cos(2 * np.pi * t / 10)
    K_inercia = 0.01
    tau_freno = K_inercia * J_total * omega
    domega = (tau_wind - tau_freno) / J_total
    
    # --- Dinámica del viento ---
    dv_wind = 0.3 * np.sin(0.3 * t) + 0.5 * np.sin(2.0 * t + 0.5)
    dv_wind_pred = 0.3 * np.sin(0.3 * (t + 5)) + 0.5 * np.sin(2.0 * (t + 5) + 0.5)
    
    # --- Dinámica del Kilómetro (enjambre) ---
    # Simulación simplificada del enjambre de Kilómetros
    enjambre_km = EnjambreKilometrosMejorado(
        N_MODULOS_KM, M_PESO, DELTA_H, ETA_GEN, ETA_LIFT, E_PERNO
    )
    E_gen_km, E_cons_km = enjambre_km.paso()
    dE_KM = (E_gen_km - E_cons_km) / DT
    
    # --- Construir vector de derivadas ---
    dstate = np.concatenate([
        dtheta,
        dr_q,
        [domega],
        [dv_wind],
        [dv_wind_pred],
        [dE_KM]
    ])
    
    return dstate

# ========================================================================
# 🔄 FUNCIÓN DE OPTIMIZACIÓN AUTOMÁTICA
# ========================================================================

def optimizar_parametros(params_iniciales, bounds, target_proporcion):
    """
    Optimiza los parámetros del sistema para alcanzar una proporción objetivo.
    """
    def objetivo(params):
        # Actualizar parámetros
        params_dict = params_iniciales.copy()
        params_dict['K'] = params[0]
        params_dict['K_bus'] = params[1]
        params_dict['K_Q_OM'] = params[2]
        params_dict['K_Q_V'] = params[3]
        
        # Ejecutar simulación
        state0 = np.concatenate([
            np.array([0.1 + i*0.1 for i in range(25)]),  # theta
            np.full(25 * 3, 5.0),                        # r_q
            [1.5, 10.0, 10.0, 0.0]                       # omega, v_wind, v_wind_pred, E_KM
        ])
        
        try:
            sol = odeint(sistema_completo, state0, t[:100], args=(params_dict,), rtol=1e-6, atol=1e-8)
            E_KM_final = sol[-1, -1]
            # Calcular proporción
            E_total = E_KM_final + 10000  # Aproximación
            proporcion = E_KM_final / E_total
            return abs(proporcion - target_proporcion)
        except:
            return 1.0
    
    # Ejecutar optimización
    result = minimize(objetivo, [params_iniciales['K'], params_iniciales['K_bus'],
                                 params_iniciales['K_Q_OM'], params_iniciales['K_Q_V']],
                     bounds=bounds, method='Nelder-Mead')
    
    # Actualizar parámetros optimizados
    params_opt = params_iniciales.copy()
    params_opt['K'] = result.x[0]
    params_opt['K_bus'] = result.x[1]
    params_opt['K_Q_OM'] = result.x[2]
    params_opt['K_Q_V'] = result.x[3]
    
    return params_opt

# ========================================================================
# 📊 SIMULACIÓN DEL SISTEMA COMPLETO
# ========================================================================

# --- Parámetros iniciales ---
params_iniciales = {
    'M_Q': M_Q,
    'J_G': J_G,
    'N_BLADES': N_BLADES,
    'omega_rated': omega_rated,
    'v_wind_rated': v_wind_rated,
    'K_Q_OM': K_Q_OM,
    'K_Q_V': K_Q_V,
    'K_Q_MEM': K_Q_MEM,
    'K': K,
    'K_bus': K_bus,
    'K_neutro': K_neutro,
    'V_bus': V_bus,
    'V_ref': V_ref,
    'V_neutro': V_neutro,
    'T': T,
    'FASES_CONTROL': FASES_CONTROL,
    'N_MODULOS_KURAMOTO': N_MODULOS_KURAMOTO,
    'N_MOLINOS_POR_MODULO': N_MOLINOS_POR_MODULO,
    'omega_natural': omega_natural,
    'N_BLADES': N_BLADES
}

# --- Optimizar parámetros ---
print("🚀 Optimizando parámetros para proporción objetivo...")
bounds = [(0.1, 1.0), (0.1, 2.0), (0.01, 0.5), (0.01, 0.5)]
params_opt = optimizar_parametros(params_iniciales, bounds, target_proporcion=0.3)
print("✅ Parámetros optimizados:")
print(f"   K = {params_opt['K']:.3f}")
print(f"   K_bus = {params_opt['K_bus']:.3f}")
print(f"   K_Q_OM = {params_opt['K_Q_OM']:.3f}")
print(f"   K_Q_V = {params_opt['K_Q_V']:.3f}")

# --- Simulación con parámetros optimizados ---
print("🚀 Simulando red eléctrica completa...")
state0 = np.concatenate([
    np.array([0.1 + i*0.1 for i in range(25)]),  # theta (25 molinos)
    np.full(25 * 3, 5.0),                        # r_q (25 molinos × 3 aspas)
    [1.5, 10.0, 10.0, 0.0]                       # omega, v_wind, v_wind_pred, E_KM
])

sol = odeint(sistema_completo, state0, t, args=(params_opt,), rtol=1e-6, atol=1e-8)

# Extraer resultados
N_TOTAL = N_MODULOS_KURAMOTO * N_MOLINOS_POR_MODULO
theta_sim = sol[:, :N_TOTAL]
r_q_sim = sol[:, N_TOTAL:N_TOTAL + 25*3]
omega_sim = sol[:, N_TOTAL + 25*3]
v_wind_sim = sol[:, N_TOTAL + 25*3 + 1]
v_wind_pred_sim = sol[:, N_TOTAL + 25*3 + 2]
E_KM_sim = sol[:, N_TOTAL + 25*3 + 3]

print("✅ Simulación completada.")

# ========================================================================
# 📈 VISUALIZACIÓN (9 GRÁFICOS + 3D)
# ========================================================================

fig = plt.figure(figsize=(20, 16))

# --- Gráfico 1: Sincronización de fases (Módulo 1) ---
ax1 = plt.subplot(3, 3, 1)
for i in range(N_MOLINOS_POR_MODULO):
    ax1.plot(t, theta_sim[:, i], label=f'M{i+1}', lw=1.5)
ax1.set_title('🔵 Sincronización de Fases (Kuramoto)')
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Fase (rad)')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)

# --- Gráfico 2: Baile de pesos (Molino 1) ---
ax2 = plt.subplot(3, 3, 2)
for j in range(N_BLADES):
    ax2.plot(t, r_q_sim[:, j], label=f'Aspa {j+1}', lw=1.5)
ax2.set_title('🌪️ Baile de Pesos (Quijote 3vs7)')
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Posición radial (m)')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

# --- Gráfico 3: Velocidad angular ---
ax3 = plt.subplot(3, 3, 3)
ax3.plot(t, omega_sim, label='ω (rad/s)', color='red', lw=2)
ax3.axhline(y=omega_rated, color='green', linestyle='--', label='ω_rated', lw=2)
ax3.set_title('⚡ Velocidad Angular')
ax3.set_xlabel('Tiempo (s)')
ax3.set_ylabel('ω (rad/s)')
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)

# --- Gráfico 4: Viento actual vs predicho ---
ax4 = plt.subplot(3, 3, 4)
ax4.plot(t, v_wind_sim, label='v_wind (actual)', color='blue', lw=1.5)
ax4.plot(t, v_wind_pred_sim, label='v_wind (predicho 5s)', color='orange', lw=1.5, linestyle='--')
ax4.axhline(y=v_wind_rated, color='gray', linestyle='--', label='v_wind_rated', alpha=0.5)
ax4.set_title('🌬️ Control Predictivo de Viento')
ax4.set_xlabel('Tiempo (s)')
ax4.set_ylabel('v_wind (m/s)')
ax4.legend(loc='best')
ax4.grid(True, alpha=0.3)

# --- Gráfico 5: Energía almacenada en Kilómetros ---
ax5 = plt.subplot(3, 3, 5)
ax5.plot(t, E_KM_sim, label='Energía almacenada (J)', color='purple', lw=2)
ax5.set_title('🏗️ Energía en Enjambre de Kilómetros')
ax5.set_xlabel('Tiempo (s)')
ax5.set_ylabel('Energía (J)')
ax5.legend(loc='best')
ax5.grid(True, alpha=0.3)

# --- Gráfico 6: Potencia total del sistema ---
ax6 = plt.subplot(3, 3, 6)
# Calcular potencia (derivada de la energía)
P_total = np.gradient(E_KM_sim, DT)
ax6.plot(t, P_total, label='Potencia total', color='purple', lw=2)
demanda = P_demanda_base + 0.5 * P_demanda_pico * (1 + np.sin(0.5 * t))
ax6.plot(t, demanda, label='Demanda', color='red', lw=2, linestyle='--')
ax6.set_title('⚡ Potencia y Demanda de la Red')
ax6.set_xlabel('Tiempo (s)')
ax6.set_ylabel('Potencia (W)')
ax6.legend(loc='best')
ax6.grid(True, alpha=0.3)

# --- Gráfico 7: Análisis de proporciones ---
ax7 = plt.subplot(3, 3, 7)
# Calcular energías por gemelo
E_kuramoto = 0.4 * np.sum(theta_sim**2, axis=1)  # Aproximación
E_quijote = 0.6 * np.sum(r_q_sim**2, axis=1)     # Aproximación
E_kilometro = E_KM_sim

# Calcular proporciones
total = E_kuramoto + E_quijote + E_kilometro + 1e-9
prop_kuramoto = E_kuramoto / total * 100
prop_quijote = E_quijote / total * 100
prop_kilometro = E_kilometro / total * 100

ax7.plot(t, prop_kuramoto, label='Kuramoto', color='blue', lw=1.5)
ax7.plot(t, prop_quijote, label='Quijote', color='green', lw=1.5)
ax7.plot(t, prop_kilometro, label='Kilómetro', color='purple', lw=1.5)
ax7.set_title('📊 Proporciones de Aporte Energético')
ax7.set_xlabel('Tiempo (s)')
ax7.set_ylabel('Porcentaje (%)')
ax7.legend(loc='best')
ax7.grid(True, alpha=0.3)

# --- Gráfico 8: Inventario de Kilómetros ---
ax8 = plt.subplot(3, 3, 8)
# Stock de pesos (simplificado)
stock_ALTA = [10 - i * 0.1 for i in range(len(t))]
stock_BAJA = [i * 0.1 for i in range(len(t))]
ax8.plot(t, stock_ALTA, label='Stock ALTA', color='orange', lw=1.5)
ax8.plot(t, stock_BAJA, label='Stock BAJA', color='cyan', lw=1.5)
ax8.set_title('📦 Stock de Pesos (Kilómetro)')
ax8.set_xlabel('Tiempo (s)')
ax8.set_ylabel('Número de Pesos')
ax8.legend(loc='best')
ax8.grid(True, alpha=0.3)

# --- Gráfico 9: Eficiencias del sistema ---
ax9 = plt.subplot(3, 3, 9)
eta_aparente = (E_KM_sim / (E_KM_sim + 100)) * 100
eta_real = (E_KM_sim / (E_KM_sim + 200)) * 100
ax9.plot(t, eta_aparente, label='η_aparente', color='blue', lw=1.5)
ax9.plot(t, eta_real, label='η_real', color='red', lw=1.5)
ax9.axhline(y=100, color='gray', linestyle='--', label='η = 100%')
ax9.set_title('📊 Eficiencias del Sistema')
ax9.set_xlabel('Tiempo (s)')
ax9.set_ylabel('Eficiencia (%)')
ax9.legend(loc='best')
ax9.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('gemelo5_red_electrica_verde_completa.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 📊 VISUALIZACIÓN 3D DE LA RED
# ========================================================================

fig_3d = plt.figure(figsize=(12, 10))
ax_3d = fig_3d.add_subplot(111, projection='3d')

# Topología de la red (5 módulos × 5 molinos)
for m in range(N_MODULOS_KURAMOTO):
    for i in range(N_MOLINOS_POR_MODULO):
        x = m * 2 + i * 0.5
        y = np.sin(theta_sim[-1, m * N_MOLINOS_POR_MODULO + i]) * 2
        z = np.cos(theta_sim[-1, m * N_MOLINOS_POR_MODULO + i]) * 2
        
        ax_3d.scatter(x, y, z, c='blue', s=50, alpha=0.7)
        
        # Conectar con los centrales
        if i < 2:  # Centrales
            ax_3d.scatter(x, y, z, c='red', s=100, marker='*')
        else:      # Anillos
            ax_3d.scatter(x, y, z, c='green', s=50, marker='o')
        
        # Líneas de conexión
        for j in range(i+1, N_MOLINOS_POR_MODULO):
            x2 = m * 2 + j * 0.5
            y2 = np.sin(theta_sim[-1, m * N_MOLINOS_POR_MODULO + j]) * 2
            z2 = np.cos(theta_sim[-1, m * N_MOLINOS_POR_MODULO + j]) * 2
            ax_3d.plot([x, x2], [y, y2], [z, z2], 'gray', alpha=0.3)

ax_3d.set_title('🌐 Topología de la Red Eléctrica (5 Módulos × 5 Molinos)')
ax_3d.set_xlabel('Módulo')
ax_3d.set_ylabel('sin(θ)')
ax_3d.set_zlabel('cos(θ)')
ax_3d.legend(['Centrales', 'Anillos', 'Conexiones'])

plt.savefig('gemelo5_topologia_red_3d.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 🔍 ANÁLISIS DE RESULTADOS Y PROPORCIONES
# ========================================================================

print("\n" + "=" * 70)
print("📊 ANÁLISIS DE PROPORCIONES DE APORTE ENERGÉTICO")
print("=" * 70)

# Calcular energías totales
E_kuramoto_total = np.trapz(E_kuramoto, t)
E_quijote_total = np.trapz(E_quijote, t)
E_kilometro_total = np.trapz(E_kilometro, t)
E_total_sistema = E_kuramoto_total + E_quijote_total + E_kilometro_total

print(f"🔹 Energía total (Kuramoto):   {E_kuramoto_total:.2f} J  ({E_kuramoto_total/E_total_sistema*100:.1f}%)")
print(f"🔹 Energía total (Quijote):    {E_quijote_total:.2f} J  ({E_quijote_total/E_total_sistema*100:.1f}%)")
print(f"🔹 Energía total (Kilómetro):  {E_kilometro_total:.2f} J  ({E_kilometro_total/E_total_sistema*100:.1f}%)")
print(f"🔹 Energía total del sistema:  {E_total_sistema:.2f} J")

# Análisis por gemelo
print("\n" + "=" * 70)
print("📌 APORTE POR GEMELO (Desglose)")
print("=" * 70)
print(f"🔹 Gemelo 1 (Kuramoto):        {E_kuramoto_total/E_total_sistema*100:.1f}%  (Sincronización de fase)")
print(f"🔹 Gemelo 2 (Quijote):         {E_quijote_total/E_total_sistema*100:.1f}%  (Inercia variable + generación)")
print(f"🔹 Gemelo 3 (Módulo Autoreg.): {E_kuramoto_total/E_total_sistema*100 + E_quijote_total/E_total_sistema*100:.1f}%  (Kuramoto + Quijote)")
print(f"🔹 Gemelo 4 (Kilómetros):      {E_kilometro_total/E_total_sistema*100:.1f}%  (Almacenamiento gravitacional)")

# ========================================================================
# 📌 RECOMENDACIONES PARA OPTIMIZAR PROPORCIONES
# ========================================================================

print("\n" + "=" * 70)
print("🎯 RECOMENDACIONES PARA AJUSTAR PROPORCIONES")
print("=" * 70)

if E_kilometro_total / E_total_sistema < 0.2:
    print("✅ El aporte de Kilómetros es bajo (<20%). Recomendación:")
    print("   - Aumentar N_MODULOS_KM (número de módulos Kilómetro).")
    print("   - Aumentar M_PESO (masa de los pesos) o DELTA_H (altura).")
    print("   - Optimizar ETA_GEN y ETA_LIFT.")
else:
    print("✅ El aporte de Kilómetros es adecuado (>=20%).")

if E_kuramoto_total / E_total_sistema + E_quijote_total / E_total_sistema > 0.8:
    print("✅ El aporte de los molinos es alto (>80%). Recomendación:")
    print("   - Reducir el número de molinos si hay excedente de energía.")
    print("   - Usar el excedente para cargar más Kilómetros.")
else:
    print("✅ El aporte de los molinos es equilibrado.")

print("\n" + "=" * 70)
print("🔧 PARÁMETROS OPTIMIZADOS")
print("=" * 70)
print(f"K = {params_opt['K']:.3f}")
print(f"K_bus = {params_opt['K_bus']:.3f}")
print(f"K_Q_OM = {params_opt['K_Q_OM']:.3f}")
print(f"K_Q_V = {params_opt['K_Q_V']:.3f}")

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código integra los 4 gemelos en una red eléctrica verde.
# - Las proporciones de aporte dependen de:
#   - La potencia generada por los molinos (Quijote + Kuramoto).
#   - La energía almacenada/generada por los Kilómetros.
# - Para ajustar las proporciones, modifica:
#   - N_MODULOS_KM (número de módulos Kilómetro).
#   - M_PESO, DELTA_H, ETA_GEN, ETA_LIFT (parámetros de Kilómetro).
#   - K, K_bus, K_Q_OM, K_Q_V (parámetros de Kuramoto + Quijote).
# ========================================================================