# 🏗️ GEMelo 6: Megabatería de 400 Kilómetros + 150 Módulos de Molinos
# ========================================================================
# Descripción: Simulación completa de una megabatería con:
#             - 400 módulos Kilómetro (almacenamiento gravitacional).
#             - 150 módulos de molinos (100 con Quijote + 50 sin Quijote).
#             - Control inteligente de carga/descarga.
#             - Vientos diferentes por módulo.
# Autor: Lee Arbusto & Víctor M. A. (espiradesombra)
# Fecha: 21-08-2026
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.stats import norm

# ========================================================================
# 📌 PARÁMETROS GLOBALES (Megabatería + Molinos)
# ========================================================================

# --- Parámetros de Kilómetro ---
N_KM = 400               # Número de módulos Kilómetro
m_peso = 10.0            # Masa de cada peso (kg)
g = 9.81                # Aceleración gravitatoria (m/s²)
delta_h = 15.0           # Altura de bajada/subida (m)
eta_gen = 0.85           # Eficiencia en generación
eta_lift = 0.90          # Eficiencia en reset
E_perno = 1.5            # Energía por perno (J)
stock_ALTA_inicial = 10 # Stock inicial de pesos en ALTA (por módulo)
stock_BAJA_inicial = 0  # Stock inicial de pesos en BAJA (por módulo)

# --- Parámetros de Molinos ---
N_modulos_molinos = 150  # 100 con Quijote + 50 sin Quijote
N_molinos_por_modulo = 5 # Molinos por módulo
N_molinos_Quijote = 100  # Módulos con Quijote (500 molinos)
N_molinos_sin_Quijote = 50 # Módulos sin Quijote (250 molinos)

# Potencia nominal por molino (W)
P_nominal_Quijote = 10000.0  # 10 kW (con Quijote)
P_nominal_sin_Quijote = 8000.0 # 8 kW (sin Quijote)

# Factor de carga (viento)
factor_carga = 0.35  # 35%

# --- Parámetros de simulación ---
t_sim = 120.0          # Tiempo total (s)
dt = 0.1              # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)  # Vector de tiempo

# --- Demanda de la red ---
P_demanda_base = 2.45e6  # 2.45 MW (demanda base)
P_demanda_pico = 3.5e6   # 3.5 MW (demanda pico)
# Demanda variable (oscilante)
P_demanda = P_demanda_base + 0.5e6 * np.sin(0.1 * t)  # Oscila entre 2.45 MW y 3.5 MW

# ========================================================================
# 🔬 CLASES PARA LOS COMPONENTES
# ========================================================================

class ModuloKilometro:
    """Clase para un módulo Kilómetro individual."""
    
    def __init__(self, id_modulo, m_peso, g, delta_h, eta_gen, eta_lift, E_perno):
        self.id = id_modulo
        self.m_peso = m_peso
        self.g = g
        self.delta_h = delta_h
        self.eta_gen = eta_gen
        self.eta_lift = eta_lift
        self.E_perno = E_perno
        
        # Estado inicial
        self.cota = 'ALTA'
        self.n_pesos = 3
        self.stock_ALTA = stock_ALTA_inicial
        self.stock_BAJA = stock_BAJA_inicial
        self.energia_generada = 0.0
        self.energia_consumida = 0.0
        self.ciclos_completados = 0
        self.potencia_actual = 0.0
        
    def ciclo(self):
        """Ejecuta un ciclo del módulo Kilómetro."""
        if self.cota == 'ALTA' and self.n_pesos == 3 and self.stock_ALTA > 0:
            # Enganche: añadir 1 peso
            self.n_pesos = 4
            self.stock_ALTA -= 1
            self.energia_consumida += 4 * self.E_perno
            self.cota = 'BAJADA'
            
        elif self.cota == 'BAJADA' and self.n_pesos == 4:
            # Bajada: generar energía
            E_gen = self.eta_gen * self.m_peso * self.g * self.delta_h
            self.energia_generada += E_gen
            self.potencia_actual = E_gen / dt  # Potencia instantánea (W)
            self.cota = 'BAJA'
            self.ciclos_completados += 1
            
        elif self.cota == 'BAJA' and self.n_pesos == 4:
            # Entrega: soltar 1 peso
            self.n_pesos = 3
            self.stock_BAJA += 1
            self.energia_consumida += 4 * self.E_perno
            self.cota = 'SUBIDA'
            
        elif self.cota == 'SUBIDA' and self.n_pesos == 3:
            # Subida: flotar
            self.cota = 'ALTA'
            
        elif self.stock_ALTA == 0:
            # Modo pausa: resetear potencial
            self.cota = 'PAUSA'
            
    def reset_potencial(self):
        """Reset del potencial de flotación."""
        if self.stock_BAJA > 0 and self.cota == 'PAUSA':
            E_reset = (self.m_peso * self.g * self.delta_h) / self.eta_lift
            self.energia_consumida += E_reset
            self.stock_BAJA -= 1
            self.stock_ALTA += 1
            self.cota = 'ALTA'
            
    def get_estado(self):
        return {
            'id': self.id,
            'cota': self.cota,
            'n_pesos': self.n_pesos,
            'stock_ALTA': self.stock_ALTA,
            'stock_BAJA': self.stock_BAJA,
            'E_generada': self.energia_generada,
            'E_consumida': self.energia_consumida,
            'potencia_actual': self.potencia_actual,
            'ciclos': self.ciclos_completados
        }


class ModuloMolinos:
    """Clase para un módulo de 5 molinos (con o sin Quijote)."""
    
    def __init__(self, id_modulo, tiene_quijote, P_nominal, factor_carga):
        self.id = id_modulo
        self.tiene_quijote = tiene_quijote
        self.P_nominal = P_nominal
        self.factor_carga = factor_carga
        self.v_wind = 12.0  # Velocidad del viento inicial (m/s)
        self.P_actual = 0.0  # Potencia actual (W)
        
    def actualizar_viento(self, t):
        """Actualiza la velocidad del viento (oscilante)."""
        # Cada módulo tiene un perfil de viento diferente
        if self.tiene_quijote:
            self.v_wind = 12.0 + 3.0 * np.sin(0.1 * t + self.id * 0.2)
        else:
            self.v_wind = 8.0 + 2.0 * np.sin(0.15 * t + self.id * 0.3)
        
    def calcular_potencia(self):
        """Calcula la potencia generada por el módulo."""
        # Potencia proporcional al cubo de la velocidad del viento (ley de Betz)
        self.P_actual = self.P_nominal * (self.v_wind / 12.0)**3 * self.factor_carga
        return self.P_actual


class Megabateria:
    """Clase para la megabatería completa (Kilómetros + Molinos)."""
    
    def __init__(self, N_KM, N_modulos_molinos, N_molinos_Quijote, P_nominal_Quijote, P_nominal_sin_Quijote, factor_carga):
        # Inicializar módulos Kilómetro
        self.kilometros = [ModuloKilometro(i, m_peso, g, delta_h, eta_gen, eta_lift, E_perno) 
                          for i in range(N_KM)]
        
        # Inicializar módulos de molinos
        self.molinos = []
        for i in range(N_modulos_molinos):
            if i < N_molinos_Quijote:
                self.molinos.append(ModuloMolinos(i, True, P_nominal_Quijote, factor_carga))
            else:
                self.molinos.append(ModuloMolinos(i, False, P_nominal_sin_Quijote, factor_carga))
        
        # Estado global
        self.energia_total_generada = 0.0
        self.energia_total_consumida = 0.0
        self.potencia_total_molinos = 0.0
        self.potencia_total_KM = 0.0
        
    def paso(self, t, P_demanda):
        """Ejecuta un paso de simulación."""
        # Actualizar molinos
        self.potencia_total_molinos = 0.0
        for modulo in self.molinos:
            modulo.actualizar_viento(t)
            self.potencia_total_molinos += modulo.calcular_potencia()
        
        # Actualizar Kilómetros
        self.potencia_total_KM = 0.0
        for km in self.kilometros:
            km.ciclo()
            km.reset_potencial()
            self.potencia_total_KM += km.potencia_actual
            
        # Control inteligente: cargar/descargar Kilómetros según demanda
        if self.potencia_total_molinos > P_demanda:
            # Excedente: cargar Kilómetros (aumentar stock_ALTA)
            P_excedente = self.potencia_total_molinos - P_demanda
            # Simplificación: distribuir excedente entre Kilómetros
            for km in self.kilometros:
                if km.stock_ALTA < stock_ALTA_inicial:
                    km.stock_ALTA += P_excedente * dt / (m_peso * g * delta_h * N_KM)
        else:
            # Déficit: descargar Kilómetros (aumentar stock_BAJA)
            P_deficit = P_demanda - self.potencia_total_molinos
            for km in self.kilometros:
                if km.stock_BAJA > 0:
                    km.stock_BAJA -= P_deficit * dt / (m_peso * g * delta_h * N_KM)
        
        # Actualizar energía total
        self.energia_total_generada = sum(km.energia_generada for km in self.kilometros)
        self.energia_total_consumida = sum(km.energia_consumida for km in self.kilometros)
        
    def get_estado_global(self):
        return {
            'P_total_molinos': self.potencia_total_molinos,
            'P_total_KM': self.potencia_total_KM,
            'E_total_generada': self.energia_total_generada,
            'E_total_consumida': self.energia_total_consumida,
            'stock_ALTA_promedio': np.mean([km.stock_ALTA for km in self.kilometros]),
            'stock_BAJA_promedio': np.mean([km.stock_BAJA for km in self.kilometros])
        }

# ========================================================================
# 📊 SIMULACIÓN DE LA MEGABATERÍA
# ========================================================================

# Inicializar megabatería
megabateria = Megabateria(
    N_KM, N_modulos_molinos, N_molinos_Quijote, 
    P_nominal_Quijote, P_nominal_sin_Quijote, factor_carga
)

# Arrays para almacenar resultados
P_molinos_hist = np.zeros(len(t))
P_KM_hist = np.zeros(len(t))
P_demanda_hist = P_demanda
E_KM_hist = np.zeros(len(t))
stock_ALTA_hist = np.zeros(len(t))
stock_BAJA_hist = np.zeros(len(t))

# Simulación paso a paso
for i, tiempo in enumerate(t):
    megabateria.paso(tiempo, P_demanda[i])
    estado = megabateria.get_estado_global()
    
    P_molinos_hist[i] = estado['P_total_molinos']
    P_KM_hist[i] = estado['P_total_KM']
    E_KM_hist[i] = estado['E_total_generada']
    stock_ALTA_hist[i] = estado['stock_ALTA_promedio']
    stock_BAJA_hist[i] = estado['stock_BAJA_promedio']

# ========================================================================
# 📈 VISUALIZACIÓN
# ========================================================================

plt.figure(figsize=(20, 15))

# --- Gráfico 1: Potencia generada por molinos y Kilómetros vs. Demanda ---
plt.subplot(4, 1, 1)
plt.plot(t, P_molinos_hist / 1e6, label='Potencia Molinos (MW)', color='blue')
plt.plot(t, P_KM_hist / 1e6, label='Potencia Kilómetros (MW)', color='green')
plt.plot(t, P_demanda_hist / 1e6, label='Demanda (MW)', color='red', linestyle='--')
plt.title('Potencia Generada vs. Demanda de la Red')
plt.xlabel('Tiempo (s)')
plt.ylabel('Potencia (MW)')
plt.legend()
plt.grid(True)

# --- Gráfico 2: Proporción de aporte energético ---
plt.subplot(4, 1, 2)
P_total = P_molinos_hist + P_KM_hist
porcentaje_molinos = (P_molinos_hist / P_total) * 100
porcentaje_KM = (P_KM_hist / P_total) * 100
plt.plot(t, porcentaje_molinos, label='% Molinos', color='blue')
plt.plot(t, porcentaje_KM, label='% Kilómetros', color='green')
plt.title('Proporción de Aporte Energético')
plt.xlabel('Tiempo (s)')
plt.ylabel('Porcentaje (%)')
plt.legend()
plt.grid(True)

# --- Gráfico 3: Energía generada por Kilómetros ---
plt.subplot(4, 1, 3)
plt.plot(t, E_KM_hist / 1e6, label='Energía Generada por Kilómetros (MJ)', color='purple')
plt.title('Energía Generada por Kilómetros')
plt.xlabel('Tiempo (s)')
plt.ylabel('Energía (MJ)')
plt.legend()
plt.grid(True)

# --- Gráfico 4: Stock promedio de pesos (ALTA/BAJA) ---
plt.subplot(4, 1, 4)
plt.plot(t, stock_ALTA_hist, label='Stock ALTA (promedio)', color='orange')
plt.plot(t, stock_BAJA_hist, label='Stock BAJA (promedio)', color='cyan')
plt.title('Stock Promedio de Pesos (ALTA/BAJA)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Número de Pesos')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('megabateria_400KM_150molinos_simulacion.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 🔍 ANÁLISIS DE RESULTADOS
# ========================================================================

# Energía total generada
E_molinos_total = np.trapz(P_molinos_hist, t)  # Energía total de molinos (J)
E_KM_total = np.trapz(P_KM_hist, t)  # Energía total de Kilómetros (J)
E_total = E_molinos_total + E_KM_total  # Energía total (J)

# Proporciones de aporte
porcentaje_molinos_total = (E_molinos_total / E_total) * 100
porcentaje_KM_total = (E_KM_total / E_total) * 100

print("=" * 80)
print("📊 ANÁLISIS DE LA MEGABATERÍA (400 Kilómetros + 150 Módulos de Molinos)")
print("=" * 80)
print(f"🔹 Energía total generada por Molinos: {E_molinos_total / 1e6:.2f} MJ")
print(f"🔹 Energía total generada por Kilómetros: {E_KM_total / 1e6:.2f} MJ")
print(f"🔹 Energía total del sistema: {E_total / 1e6:.2f} MJ")
print(f"🔹 Proporción de Molinos: {porcentaje_molinos_total:.2f}%")
print(f"🔹 Proporción de Kilómetros: {porcentaje_KM_total:.2f}%")

# --- Desglose por tipo de molino ---
P_Quijote_total = N_molinos_Quijote * N_molinos_por_modulo * P_nominal_Quijote * factor_carga
P_sin_Quijote_total = N_molinos_sin_Quijote * N_molinos_por_modulo * P_nominal_sin_Quijote * factor_carga

print("\n" + "=" * 80)
print("📌 DESGLOSE POR TIPO DE MOLINO")
print("=" * 80)
print(f"🔹 Potencia total de molinos con Quijote: {P_Quijote_total / 1e6:.2f} MW")
print(f"🔹 Potencia total de molinos sin Quijote: {P_sin_Quijote_total / 1e6:.2f} MW")
print(f"🔹 Potencia total de molinos: {(P_Quijote_total + P_sin_Quijote_total) / 1e6:.2f} MW")

# --- Análisis de Kilómetros ---
E_KM_por_modulo = E_KM_total / N_KM
print("\n" + "=" * 80)
print("📌 ANÁLISIS DE KILÓMETROS")
print("=" * 80)
print(f"🔹 Energía por módulo Kilómetro: {E_KM_por_modulo:.2f} J")
print(f"🔹 Potencia pico del enjambre de Kilómetros: {np.max(P_KM_hist):.2f} W")
print(f"🔹 Stock ALTA promedio final: {stock_ALTA_hist[-1]:.2f} pesos")
print(f"🔹 Stock BAJA promedio final: {stock_BAJA_hist[-1]:.2f} pesos")

# ========================================================================
# 🎯 RECOMENDACIONES PARA OPTIMIZAR LA MEGABATERÍA
# ========================================================================

print("\n" + "=" * 80)
print("🎯 RECOMENDACIONES")
print("=" * 80)

if porcentaje_KM_total < 20:
    print("✅ El aporte de Kilómetros es bajo (<20%). Recomendación:")
    print("   - Aumentar N_KM (número de módulos Kilómetro).")
    print("   - Aumentar m_peso o delta_h para generar más energía por ciclo.")
else:
    print("✅ El aporte de Kilómetros es adecuado (>=20%).")

if porcentaje_molinos_total > 80:
    print("✅ El aporte de molinos es alto (>80%). Recomendación:")
    print("   - Reducir el número de molinos si hay excedente de energía.")
    print("   - Usar el excedente para cargar más Kilómetros.")
else:
    print("✅ El aporte de molinos es equilibrado.")

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código simula una megabatería de 400 Kilómetros + 150 módulos de molinos.
# - Los molinos con Quijote generan más potencia (10 kW) que los sin Quijote (8 kW).
# - Los Kilómetros actúan como batería de respaldo, cargando/descargando según la demanda.
# - Cada módulo de molinos tiene un perfil de viento diferente.
# - Para validar experimentalmente, necesitarás datos reales de molinos y Kilómetros.
# ========================================================================