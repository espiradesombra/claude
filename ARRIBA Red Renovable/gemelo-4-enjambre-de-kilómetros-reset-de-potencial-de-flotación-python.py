# 🏗️ GEMelo 4: Enjambre de Kilómetros (Reset de Potencial de Flotación por Perneado)
# ========================================================================
# Descripción: Simulación de un enjambre de módulos Kilómetro para almacenamiento
#              gravitacional con reset de potencial de flotación.
# Autor: Lee Arbusto & Víctor M. A. (espiradesombra)
# Fecha: 21-08-2026
# Basado en: KILOMETRO_SIM_v1_PARAMETROS.md (repo espiradesombra/claude)
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ========================================================================
# 📌 PARÁMETROS DEL SISTEMA (Basados en KILOMETRO_SIM_v1_PARAMETROS.md)
# ========================================================================

# --- Parámetros geométricos ---
R_rotacion = 5.0       # Radio de rotación (m)
L_carril = 2.0         # Longitud del carril (m)
N_masas = 3           # Número de masas por módulo Kilómetro
m_peso = 500.0        # Masa de cada peso (kg)

# --- Parámetros de operación ---
r_min = 3.0           # Posición radial mínima (m)
r_max = 7.0           # Posición radial máxima (m)
omega_rot = 1.0       # Velocidad angular (rad/s)
T_ciclo = 2 * np.pi / omega_rot  # Período de ciclo (s)

# --- Parámetros de energía ---
g = 9.81              # Aceleración gravitatoria (m/s²)
rho_fluido = 1000.0   # Densidad del fluido (agua, kg/m³)
mu_k = 0.05           # Coeficiente de fricción
c = 100.0             # Factor de amortiguamiento (N·s/m)

# --- Parámetros de control ---
E_perno = 1.5          # Energía por perno (J)
eta_gen = 0.85        # Eficiencia en generación
eta_lift = 0.90       # Eficiencia en reset

# --- Parámetros del enjambre ---
N_KM = 5              # Número de módulos Kilómetro en el enjambre
stock_ALTA = 10      # Pesos iniciales en stock ALTA (por módulo)
stock_BAJA = 0       # Pesos iniciales en stock BAJA (por módulo)

# --- Parámetros de simulación ---
t_sim = 60.0          # Tiempo total de simulación (s)
dt = 0.01             # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)  # Vector de tiempo

# ========================================================================
# 🔬 MODELO MATEMÁTICO: Módulo Kilómetro Individual
# ========================================================================

class ModuloKilometro:
    """Clase para simular un módulo Kilómetro individual."""
    
    def __init__(self, id_modulo, m_peso, g, eta_gen, eta_lift, E_perno, r_min, r_max):
        self.id = id_modulo
        self.m_peso = m_peso
        self.g = g
        self.eta_gen = eta_gen
        self.eta_lift = eta_lift
        self.E_perno = E_perno
        self.r_min = r_min
        self.r_max = r_max
        
        # Estado inicial
        self.cota = 'ALTA'  # Cota actual (ALTA o BAJA)
        self.n_pesos = 3    # Pesos actuales en el objeto (3 o 4)
        self.stock_ALTA = stock_ALTA  # Pesos en stock ALTA
        self.stock_BAJA = stock_BAJA  # Pesos en stock BAJA
        self.energia_generada = 0.0  # Energía generada acumulada (J)
        self.energia_consumida = 0.0  # Energía consumida en reset (J)
        self.ciclos_completados = 0   # Contador de ciclos
        
    def ciclo(self, delta_h):
        """
        Ejecuta un ciclo completo del módulo Kilómetro:
        1. Enganche (ALTA)
        2. Bajada (Generación)
        3. Entrega (BAJA)
        4. Subida (Flotación)
        """
        # --- Fase 1: Enganche (ALTA) ---
        if self.cota == 'ALTA' and self.n_pesos == 3:
            # Enganchar 1 peso del stock ALTA
            if self.stock_ALTA > 0:
                self.n_pesos = 4
                self.stock_ALTA -= 1
                self.energia_consumida += 4 * self.E_perno  # Coste: 4 pernos
                self.cota = 'ALTA'  # Sigue en ALTA hasta bajar
                
        # --- Fase 2: Bajada (Generación) ---
        elif self.cota == 'ALTA' and self.n_pesos == 4:
            # Bajada con 4 pesos (se hunde)
            E_gen = self.eta_gen * self.m_peso * self.g * delta_h
            self.energia_generada += E_gen
            self.cota = 'BAJA'
            self.ciclos_completados += 1
            
        # --- Fase 3: Entrega (BAJA) ---
        elif self.cota == 'BAJA' and self.n_pesos == 4:
            # Soltar el peso extra
            self.n_pesos = 3
            self.stock_BAJA += 1
            self.energia_consumida += 4 * self.E_perno  # Coste: 4 pernos
            
        # --- Fase 4: Subida (Flotación) ---
        elif self.cota == 'BAJA' and self.n_pesos == 3:
            # Subida con 3 pesos (flota)
            self.cota = 'ALTA'
            
        # --- Modo Pausa (Reset de flotación) ---
        if self.stock_ALTA == 0:
            # Si no hay pesos en ALTA, entrar en modo pausa
            self.cota = 'PAUSA'
            
    def reset_potencial(self, delta_h):
        """Reset del potencial de flotación (mover pesos de BAJA a ALTA)."""
        if self.stock_BAJA > 0 and self.cota == 'PAUSA':
            # Subir un peso de BAJA a ALTA (coste: m_peso * g * delta_h / eta_lift)
            E_reset = (self.m_peso * self.g * delta_h) / self.eta_lift
            self.energia_consumida += E_reset
            self.stock_BAJA -= 1
            self.stock_ALTA += 1
            self.cota = 'ALTA'
            
    def get_estado(self):
        """Devuelve el estado actual del módulo."""
        return {
            'id': self.id,
            'cota': self.cota,
            'n_pesos': self.n_pesos,
            'stock_ALTA': self.stock_ALTA,
            'stock_BAJA': self.stock_BAJA,
            'E_generada': self.energia_generada,
            'E_consumida': self.energia_consumida,
            'ciclos': self.ciclos_completados
        }

# ========================================================================
# 🔄 SISTEMA: Enjambre de Kilómetros
# ========================================================================

class EnjambreKilometros:
    """Clase para simular un enjambre de módulos Kilómetro."""
    
    def __init__(self, N_KM, m_peso, g, eta_gen, eta_lift, E_perno, r_min, r_max, delta_h):
        self.modulos = [ModuloKilometro(i, m_peso, g, eta_gen, eta_lift, E_perno, r_min, r_max) 
                        for i in range(N_KM)]
        self.delta_h = delta_h
        self.energia_total_generada = 0.0
        self.energia_total_consumida = 0.0
        
    def paso(self):
        """Ejecuta un paso de simulación para todos los módulos."""
        for modulo in self.modulos:
            modulo.ciclo(self.delta_h)
            modulo.reset_potencial(self.delta_h)
            
        # Actualizar energía total del enjambre
        self.energia_total_generada = sum(m.energia_generada for m in self.modulos)
        self.energia_total_consumida = sum(m.energia_consumida for m in self.modulos)
        
    def get_estado_global(self):
        """Devuelve el estado global del enjambre."""
        return {
            'E_total_generada': self.energia_total_generada,
            'E_total_consumida': self.energia_total_consumida,
            'eta_aparente': self.energia_total_generada / (4 * self.modulos[0].E_perno * sum(m.ciclos_completados for m in self.modulos)),
            'eta_real': self.energia_total_generada / self.energia_total_consumida,
            'modulos': [m.get_estado() for m in self.modulos]
        }

# ========================================================================
# 📊 SIMULACIÓN DEL ENJAMBRE
# ========================================================================

# --- Inicializar enjambre ---
delta_h = 15.0  # Altura de bajada/subida (m)
enjambre = EnjambreKilometros(N_KM, m_peso, g, eta_gen, eta_lift, E_perno, r_min, r_max, delta_h)

# --- Simulación paso a paso ---
estados = []
for _ in t:
    enjambre.paso()
    estados.append(enjambre.get_estado_global())

# ========================================================================
# 📈 VISUALIZACIÓN
# ========================================================================

plt.figure(figsize=(18, 12))

# --- Gráfico 1: Energía generada y consumida por el enjambre ---
plt.subplot(3, 1, 1)
E_gen = [estado['E_total_generada'] for estado in estados]
E_con = [estado['E_total_consumida'] for estado in estados]
plt.plot(t, E_gen, label='Energía Generada (J)', color='green')
plt.plot(t, E_con, label='Energía Consumida (J)', color='red')
plt.title('Energía Generada vs. Consumida (Enjambre de Kilómetros)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Energía (J)')
plt.legend()
plt.grid(True)

# --- Gráfico 2: Eficiencias (Aparente y Real) ---
plt.subplot(3, 1, 2)
eta_aparente = [estado['eta_aparente'] for estado in estados]
eta_real = [estado['eta_real'] for estado in estados]
plt.plot(t, eta_aparente, label='Eficiencia Aparente (η > 1)', color='blue')
plt.plot(t, eta_real, label='Eficiencia Real (η < 1)', color='orange')
plt.axhline(y=1.0, color='gray', linestyle='--', label='η = 1')
plt.title('Eficiencias del Enjambre (Aparente vs. Real)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Eficiencia')
plt.legend()
plt.grid(True)

# --- Gráfico 3: Ciclos completados por módulo ---
plt.subplot(3, 1, 3)
for i in range(N_KM):
    ciclos = [estado['modulos'][i]['ciclos'] for estado in estados]
    plt.plot(t, ciclos, label=f'Módulo {i+1}')
plt.title('Ciclos Completados por Módulo')
plt.xlabel('Tiempo (s)')
plt.ylabel('Número de Ciclos')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('enjambre_kilometros_energia_ciclos.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 🔍 ANÁLISIS DE RESULTADOS
# ========================================================================

estado_final = enjambre.get_estado_global()
print(f"🔹 Energía total generada: {estado_final['E_total_generada']:.2f} J")
print(f"🔹 Energía total consumida: {estado_final['E_total_consumida']:.2f} J")
print(f"🔹 Eficiencia aparente: {estado_final['eta_aparente']:.2f} (η > 1)")
print(f"🔹 Eficiencia real: {estado_final['eta_real']:.2f} (η < 1)")

for modulo in enjambre.modulos:
    estado = modulo.get_estado()
    print(f"🔹 Módulo {estado['id'] + 1}: {estado['ciclos']} ciclos, Stock ALTA: {estado['stock_ALTA']}, Stock BAJA: {estado['stock_BAJA']}")

# ========================================================================
# 🎯 INTEGRACIÓN CON MÓDULO AUTOREGULANTE (Kuramoto + Quijote)
# ========================================================================

# Para integrar el Gemelo 4 (Kilómetro) con el módulo autoregulante (Kuramoto + Quijote):
# 1. La energía generada por el enjambre de Kilómetros se usa para:
#    - Compensar las pérdidas en el sistema Kuramoto + Quijote.
#    - Almacenar energía excedente.
# 2. El modo pausa del Kilómetro permite resetear el potencial de flotación.

# Ejemplo de cómo sería la función integrada:
def sistema_completo(t, state, enjambre, M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V, K, K_bus, V_bus, V_ref, T):
    """
    Sistema completo: Kuramoto + Quijote + Kilómetro.
    
    Args:
        state: Estado del sistema [theta1..theta5, r1_1..r5_3, omega, v_wind, E_KM]
        enjambre: Objeto EnjambreKilometros
        ... (otros parámetros)
    """
    # Extraer estados de Kuramoto + Quijote
    theta = state[:N_molinos]
    r_q = state[N_molinos:N_molinos + N_molinos * N_blades]
    omega = state[N_molinos + N_molinos * N_blades]
    v_wind = state[N_molinos + N_molinos * N_blades + 1]
    E_KM = state[N_molinos + N_molinos * N_blades + 2]  # Energía almacenada en Kilómetros
    
    # Dinámica de Kuramoto + Quijote (igual que antes)
    dstate_dt = sistema_autoregulante(state, t, M_Q, J_G, N_blades, omega_rated, v_wind_rated, K_Q_OM, K_Q_V, K, K_bus, V_bus, V_ref, T)
    
    # Dinámica del Kilómetro: la energía generada se suma a E_KM
    enjambre.paso()
    estado_KM = enjambre.get_estado_global()
    dE_KM_dt = estado_KM['E_total_generada'] - estado_KM['E_total_consumida']
    
    # Devolver derivadas (incluyendo E_KM)
    return np.concatenate([dstate_dt, [dE_KM_dt]])

# ========================================================================
# 📌 PRÓXIMOS PASOS (Para Víctor)
# ========================================================================
# 1. Validar el enjambre de Kilómetros con datos reales (ej: tanque 40x40x100).
# 2. Ajustar parámetros (m_peso, delta_h, eta_gen, eta_lift) para optimizar la eficiencia.
# 3. Integrar con el módulo autoregulante (Kuramoto + Quijote) usando `sistema_completo`.
# 4. Prototipar un módulo Kilómetro físico para validar el reset de flotación.

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código simula un enjambre de módulos Kilómetro con reset de potencial.
# - La eficiencia aparente (η > 1) se debe a la asimetría de fase (no viola termodinámica).
# - La eficiencia real (η < 1) incluye el coste de resetear los pesos.
# - Para validar experimentalmente, necesitarás un prototipo físico (tanque + pesos).
# ========================================================================