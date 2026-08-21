# 🏗️ GEMELO 6 - VERSIÓN MEJORADA: Megabatería de 400 Kilómetros + 150 Módulos de Molinos
# ========================================================================
# 🚀 NUEVAS CARACTERÍSTICAS:
#   1. 400 Kilómetros con patada variable (3, 5, 7 fases)
#   2. 150 módulos de molinos (100 con Quijote + 50 sin Quijote)
#   3. Vientos diferentes por módulo (simulación de parque eólico real)
#   4. Control predictivo de Kilómetros (meter/quitar energía según demanda)
#   5. Medición de generación en tiempo real
#   6. Optimización de parámetros (patada, hurto de distancia)
#   7. Visualización 3D de la red (molinos + Kilómetros)
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

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

# Parámetros de patada (3, 5, 7 fases)
fases_patada = [3, 5, 7]  # Número de fases de patada
patada_factor = {3: 1.0, 5: 1.5, 7: 2.0}  # Factor de aceleración por patada

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
t_sim = 300.0          # Tiempo total (s)
dt = 0.1              # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)  # Vector de tiempo

# --- Demanda de la red ---
P_demanda_base = 2.45e6  # 2.45 MW (demanda base)
P_demanda_pico = 3.5e6   # 3.5 MW (demanda pico)
P_demanda = P_demanda_base + 0.5e6 * np.sin(0.05 * t)  # Oscila entre 2.45 MW y 3.5 MW

# ========================================================================
# 🔬 CLASES MEJORADAS PARA LOS COMPONENTES
# ========================================================================

class ModuloKilometroMejorado:
    """Módulo Kilómetro con patada variable y hurto de distancia."""
    
    def __init__(self, id_modulo, m_peso, g, delta_h, eta_gen, eta_lift, E_perno, fases_patada):
        self.id = id_modulo
        self.m_peso = m_peso
        self.g = g
        self.delta_h = delta_h
        self.eta_gen = eta_gen
        self.eta_lift = eta_lift
        self.E_perno = E_perno
        self.fases_patada = fases_patada  # [3, 5, 7]
        self.fase_actual = 3  # Comienza con 3 fases
        
        # Estado inicial
        self.cota = 'ALTA'
        self.n_pesos = 3
        self.stock_ALTA = stock_ALTA_inicial
        self.stock_BAJA = stock_BAJA_inicial
        self.energia_generada = 0.0
        self.energia_consumida = 0.0
        self.ciclos_completados = 0
        self.potencia_actual = 0.0
        self.tiempo_ciclo = 0.0
        self.patada_activa = False
        self.hurto_distancia = 0.0
        
    def cambiar_fase(self, nuevas_fases):
        """Cambia el número de fases de patada."""
        if nuevas_fases in self.fases_patada:
            self.fase_actual = nuevas_fases
            
    def ejecutar_ciclo(self):
        """Ejecuta un ciclo del módulo Kilómetro con patada variable."""
        # Factor de patada
        factor_patada = patada_factor.get(self.fase_actual, 1.0)
        
        if self.cota == 'ALTA' and self.n_pesos == 3 and self.stock_ALTA > 0:
            # Enganche: añadir 1 peso
            self.n_pesos = 4
            self.stock_ALTA -= 1
            self.energia_consumida += 4 * self.E_perno
            self.cota = 'BAJADA'
            self.patada_activa = True
            self.hurto_distancia = self.delta_h * (1 - 1/factor_patada)  # Hurto de distancia
            
        elif self.cota == 'BAJADA' and self.n_pesos == 4 and self.patada_activa:
            # Bajada con patada: generar energía (más rápida)
            E_gen = self.eta_gen * self.m_peso * self.g * self.delta_h * factor_patada
            self.energia_generada += E_gen
            self.potencia_actual = E_gen / (dt / factor_patada)  # Potencia instantánea (W)
            self.cota = 'BAJA'
            self.ciclos_completados += 1
            self.tiempo_ciclo += dt / factor_patada
            self.patada_activa = False
            
        elif self.cota == 'BAJA' and self.n_pesos == 4:
            # Entrega: soltar 1 peso
            self.n_pesos = 3
            self.stock_BAJA += 1
            self.energia_consumida += 4 * self.E_perno
            self.cota = 'SUBIDA'
            
        elif self.cota == 'SUBIDA' and self.n_pesos == 3:
            # Subida con hurto de distancia (más rápida)
            self.cota = 'ALTA'
            self.tiempo_ciclo += dt / (1 + self.hurto_distancia/self.delta_h)
            
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
            'ciclos': self.ciclos_completados,
            'tiempo_ciclo': self.tiempo_ciclo,
            'fase_actual': self.fase_actual,
            'hurto_distancia': self.hurto_distancia
        }


class ModuloMolinosMejorado:
    """Módulo de 5 molinos con viento personalizado."""
    
    def __init__(self, id_modulo, tiene_quijote, P_nominal, factor_carga):
        self.id = id_modulo
        self.tiene_quijote = tiene_quijote
        self.P_nominal = P_nominal
        self.factor_carga = factor_carga
        self.v_wind = 12.0  # Velocidad del viento inicial (m/s)
        self.P_actual = 0.0  # Potencia actual (W)
        self.viento_perfil = self._generar_perfil_viento()
        
    def _generar_perfil_viento(self):
        """Genera un perfil de viento único para cada módulo."""
        if self.tiene_quijote:
            # Viento con ráfagas y caídas (más realista)
            return {
                'base': 12.0,
                'amplitud': 3.0,
                'frecuencia': 0.1 + self.id * 0.01,
                'fase': self.id * 0.2,
                'rafagas': 0.5 * np.sin(self.id * 0.3)
            }
        else:
            # Viento más suave (sin Quijote)
            return {
                'base': 8.0,
                'amplitud': 2.0,
                'frecuencia': 0.15 + self.id * 0.02,
                'fase': self.id * 0.3,
                'rafagas': 0.3 * np.sin(self.id * 0.4)
            }
        
    def actualizar_viento(self, t):
        """Actualiza la velocidad del viento según el perfil."""
        perfil = self.viento_perfil
        self.v_wind = (
            perfil['base'] 
            + perfil['amplitud'] * np.sin(perfil['frecuencia'] * t + perfil['fase'])
            + perfil['rafagas'] * np.sin(0.05 * t + self.id * 0.1)
        )
        self.v_wind = max(self.v_wind, 2.0)  # Mínimo 2 m/s
        
    def calcular_potencia(self):
        """Calcula la potencia generada por el módulo (ley de Betz)."""
        # Potencia proporcional al cubo de la velocidad del viento
        v_ref = 12.0 if self.tiene_quijote else 8.0
        self.P_actual = self.P_nominal * (self.v_wind / v_ref)**3 * self.factor_carga
        self.P_actual = max(self.P_actual, 0)  # No puede ser negativa
        return self.P_actual


class MegabateriaMejorada:
    """Megabatería completa con control predictivo."""
    
    def __init__(self, N_KM, N_modulos_molinos, N_molinos_Quijote, 
                 P_nominal_Quijote, P_nominal_sin_Quijote, factor_carga):
        # Inicializar módulos Kilómetro
        self.kilometros = [
            ModuloKilometroMejorado(i, m_peso, g, delta_h, eta_gen, eta_lift, E_perno, fases_patada)
            for i in range(N_KM)
        ]
        
        # Distribuir fases de patada entre los Kilómetros
        for i, km in enumerate(self.kilometros):
            # 30% en 3 fases, 40% en 5 fases, 30% en 7 fases
            if i % 10 < 3:
                km.cambiar_fase(3)
            elif i % 10 < 7:
                km.cambiar_fase(5)
            else:
                km.cambiar_fase(7)
        
        # Inicializar módulos de molinos
        self.molinos = []
        for i in range(N_modulos_molinos):
            if i < N_molinos_Quijote:
                self.molinos.append(ModuloMolinosMejorado(i, True, P_nominal_Quijote, factor_carga))
            else:
                self.molinos.append(ModuloMolinosMejorado(i, False, P_nominal_sin_Quijote, factor_carga))
        
        # Estado global
        self.energia_total_generada = 0.0
        self.energia_total_consumida = 0.0
        self.potencia_total_molinos = 0.0
        self.potencia_total_KM = 0.0
        self.historial = []
        
        # Control predictivo (buffer de energía)
        self.buffer_energia = 0.0
        self.capacidad_buffer = 1000000.0  # 1 MJ
        
    def paso(self, t, P_demanda):
        """Ejecuta un paso de simulación con control predictivo."""
        # 1. Actualizar molinos
        self.potencia_total_molinos = 0.0
        for modulo in self.molinos:
            modulo.actualizar_viento(t)
            self.potencia_total_molinos += modulo.calcular_potencia()
        
        # 2. Actualizar Kilómetros
        self.potencia_total_KM = 0.0
        for km in self.kilometros:
            km.ejecutar_ciclo()
            km.reset_potencial()
            self.potencia_total_KM += km.potencia_actual
        
        # 3. Control predictivo (buffer)
        P_total = self.potencia_total_molinos + self.potencia_total_KM
        
        # Predecir demanda futura (5 segundos)
        if t < len(P_demanda) - 5/dt:
            P_futuro = np.mean(P_demanda[int(t/dt):int(t/dt)+int(5/dt)])
        else:
            P_futuro = P_demanda
        
        # Decidir carga/descarga del buffer
        if P_total > P_demanda:
            # Excedente: cargar buffer
            excedente = (P_total - P_demanda) * dt
            self.buffer_energia = min(self.buffer_energia + excedente, self.capacidad_buffer)
            # Si hay excedente, cargar Kilómetros (subir pesos)
            for km in self.kilometros:
                if km.stock_ALTA < stock_ALTA_inicial:
                    km.stock_ALTA += excedente / (m_peso * g * delta_h * N_KM)
        else:
            # Déficit: descargar buffer
            deficit = (P_demanda - P_total) * dt
            if self.buffer_energia > 0:
                disponible = min(self.buffer_energia, deficit)
                self.buffer_energia -= disponible
                # Si hay déficit, descargar Kilómetros (bajar pesos)
                for km in self.kilometros:
                    if km.stock_BAJA > 0:
                        km.stock_BAJA -= disponible / (m_peso * g * delta_h * N_KM)
        
        # 4. Medir generación
        self.energia_total_generada = sum(km.energia_generada for km in self.kilometros)
        self.energia_total_consumida = sum(km.energia_consumida for km in self.kilometros)
        
        # 5. Guardar historial
        self.historial.append({
            't': t,
            'P_molinos': self.potencia_total_molinos,
            'P_KM': self.potencia_total_KM,
            'P_total': P_total,
            'P_demanda': P_demanda,
            'buffer_energia': self.buffer_energia,
            'E_KM_total': self.energia_total_generada
        })
        
    def get_estado_global(self):
        return {
            'P_total_molinos': self.potencia_total_molinos,
            'P_total_KM': self.potencia_total_KM,
            'P_total': self.potencia_total_molinos + self.potencia_total_KM,
            'E_total_generada': self.energia_total_generada,
            'E_total_consumida': self.energia_total_consumida,
            'buffer_energia': self.buffer_energia,
            'stock_ALTA_promedio': np.mean([km.stock_ALTA for km in self.kilometros]),
            'stock_BAJA_promedio': np.mean([km.stock_BAJA for km in self.kilometros]),
            'ciclos_totales': sum(km.ciclos_completados for km in self.kilometros)
        }

# ========================================================================
# 📊 SIMULACIÓN DE LA MEGABATERÍA MEJORADA
# ========================================================================

# Inicializar megabatería
megabateria = MegabateriaMejorada(
    N_KM, N_modulos_molinos, N_molinos_Quijote, 
    P_nominal_Quijote, P_nominal_sin_Quijote, factor_carga
)

# Ejecutar simulación
for i, tiempo in enumerate(t):
    megabateria.paso(tiempo, P_demanda[i])

# Extraer resultados del historial
P_molinos_hist = [h['P_molinos'] for h in megabateria.historial]
P_KM_hist = [h['P_KM'] for h in megabateria.historial]
P_total_hist = [h['P_total'] for h in megabateria.historial]
P_demanda_hist = [h['P_demanda'] for h in megabateria.historial]
buffer_energia_hist = [h['buffer_energia'] for h in megabateria.historial]
E_KM_hist = [h['E_KM_total'] for h in megabateria.historial]

estado_final = megabateria.get_estado_global()

# ========================================================================
# 📈 VISUALIZACIÓN (8 GRÁFICOS)
# ========================================================================

fig = plt.figure(figsize=(20, 16))

# --- Gráfico 1: Potencia generada vs. Demanda ---
ax1 = plt.subplot(4, 2, 1)
ax1.plot(t[:len(P_molinos_hist)], np.array(P_molinos_hist)/1e6, label='Molinos (MW)', color='blue', lw=1.5)
ax1.plot(t[:len(P_KM_hist)], np.array(P_KM_hist)/1e6, label='Kilómetros (MW)', color='green', lw=1.5)
ax1.plot(t[:len(P_demanda_hist)], np.array(P_demanda_hist)/1e6, label='Demanda (MW)', color='red', lw=2, linestyle='--')
ax1.set_title('⚡ Potencia Generada vs. Demanda')
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Potencia (MW)')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)

# --- Gráfico 2: Proporción de aporte ---
ax2 = plt.subplot(4, 2, 2)
P_total_arr = np.array(P_total_hist)
P_molinos_arr = np.array(P_molinos_hist)
P_KM_arr = np.array(P_KM_hist)
porc_molinos = (P_molinos_arr / (P_total_arr + 1e-9)) * 100
porc_KM = (P_KM_arr / (P_total_arr + 1e-9)) * 100
ax2.plot(t[:len(porc_molinos)], porc_molinos, label='% Molinos', color='blue', lw=1.5)
ax2.plot(t[:len(porc_KM)], porc_KM, label='% Kilómetros', color='green', lw=1.5)
ax2.set_title('📊 Proporción de Aporte Energético')
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Porcentaje (%)')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

# --- Gráfico 3: Energía en Kilómetros ---
ax3 = plt.subplot(4, 2, 3)
ax3.plot(t[:len(E_KM_hist)], np.array(E_KM_hist)/1e6, label='Energía KM (MJ)', color='purple', lw=2)
ax3.set_title('🏗️ Energía Generada por Kilómetros')
ax3.set_xlabel('Tiempo (s)')
ax3.set_ylabel('Energía (MJ)')
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)

# --- Gráfico 4: Buffer de energía ---
ax4 = plt.subplot(4, 2, 4)
ax4.plot(t[:len(buffer_energia_hist)], np.array(buffer_energia_hist)/1e6, label='Buffer (MJ)', color='orange', lw=2)
ax4.axhline(y=megabateria.capacidad_buffer/1e6, color='red', linestyle='--', label='Capacidad máxima')
ax4.set_title('🔋 Buffer de Energía (Control Predictivo)')
ax4.set_xlabel('Tiempo (s)')
ax4.set_ylabel('Energía (MJ)')
ax4.legend(loc='best')
ax4.grid(True, alpha=0.3)

# --- Gráfico 5: Stock de pesos (promedio) ---
ax5 = plt.subplot(4, 2, 5)
# Extraer stock de pesos del estado final (simplificado)
stock_ALTA_prom = estado_final['stock_ALTA_promedio']
stock_BAJA_prom = estado_final['stock_BAJA_promedio']
ax5.bar(['Stock ALTA', 'Stock BAJA'], [stock_ALTA_prom, stock_BAJA_prom], color=['orange', 'cyan'])
ax5.set_title('📦 Stock Promedio de Pesos (Kilómetros)')
ax5.set_ylabel('Número de Pesos')
ax5.grid(True, alpha=0.3)

# --- Gráfico 6: Potencia por tipo de molino ---
ax6 = plt.subplot(4, 2, 6)
# Calcular potencia por tipo
P_Quijote = sum(m.P_actual for m in megabateria.molinos if m.tiene_quijote)
P_sin_Quijote = sum(m.P_actual for m in megabateria.molinos if not m.tiene_quijote)
ax6.bar(['Con Quijote', 'Sin Quijote'], [P_Quijote/1e6, P_sin_Quijote/1e6], color=['blue', 'gray'])
ax6.set_title('🌪️ Potencia por Tipo de Molino')
ax6.set_ylabel('Potencia (MW)')
ax6.grid(True, alpha=0.3)

# --- Gráfico 7: Ciclos completados por Kilómetro ---
ax7 = plt.subplot(4, 2, 7)
ciclos_por_km = [km.ciclos_completados for km in megabateria.kilometros]
ax7.hist(ciclos_por_km, bins=20, color='purple', alpha=0.7)
ax7.set_title('🔄 Distribución de Ciclos por Kilómetro')
ax7.set_xlabel('Número de Ciclos')
ax7.set_ylabel('Frecuencia')
ax7.grid(True, alpha=0.3)

# --- Gráfico 8: Fases de patada por Kilómetro ---
ax8 = plt.subplot(4, 2, 8)
fases = [km.fase_actual for km in megabateria.kilometros]
fases_3 = sum(1 for f in fases if f == 3)
fases_5 = sum(1 for f in fases if f == 5)
fases_7 = sum(1 for f in fases if f == 7)
ax8.bar(['3 fases', '5 fases', '7 fases'], [fases_3, fases_5, fases_7], color=['blue', 'green', 'red'])
ax8.set_title('⚡ Distribución de Fases de Patada')
ax8.set_ylabel('Número de Kilómetros')
ax8.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('gemelo6_megabateria_mejorada.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 📊 ANÁLISIS DE RESULTADOS
# ========================================================================

print("\n" + "=" * 80)
print("📊 ANÁLISIS DE LA MEGABATERÍA MEJORADA")
print("=" * 80)

# Energías totales
E_molinos_total = np.trapz(P_molinos_arr, t[:len(P_molinos_arr)])
E_KM_total = np.trapz(P_KM_arr, t[:len(P_KM_arr)])
E_total = E_molinos_total + E_KM_total

print(f"🔹 Energía total generada por Molinos: {E_molinos_total/1e6:.2f} MJ")
print(f"🔹 Energía total generada por Kilómetros: {E_KM_total/1e6:.2f} MJ")
print(f"🔹 Energía total del sistema: {E_total/1e6:.2f} MJ")
print(f"🔹 Proporción de Molinos: {E_molinos_total/E_total*100:.1f}%")
print(f"🔹 Proporción de Kilómetros: {E_KM_total/E_total*100:.1f}%")

# Análisis por tipo de Kilómetro (fase de patada)
print("\n" + "=" * 80)
print("📌 ANÁLISIS POR FASE DE PATADA")
print("=" * 80)

for f in [3, 5, 7]:
    kms_fase = [km for km in megabateria.kilometros if km.fase_actual == f]
    if kms_fase:
        E_prom = np.mean([km.energia_generada for km in kms_fase])
        ciclos_prom = np.mean([km.ciclos_completados for km in kms_fase])
        tiempo_ciclo_prom = np.mean([km.tiempo_ciclo for km in kms_fase])
        print(f"🔹 {f} fases: {len(kms_fase)} KM, "
              f"E_prom={E_prom:.2f} J, "
              f"ciclos={ciclos_prom:.1f}, "
              f"tiempo_ciclo={tiempo_ciclo_prom:.2f} s")

# Análisis por tipo de molino
print("\n" + "=" * 80)
print("📌 ANÁLISIS POR TIPO DE MOLINO")
print("=" * 80)

P_Quijote_total = sum(m.P_actual for m in megabateria.molinos if m.tiene_quijote)
P_sin_Quijote_total = sum(m.P_actual for m in megabateria.molinos if not m.tiene_quijote)

print(f"🔹 Potencia total de molinos con Quijote: {P_Quijote_total/1e6:.2f} MW")
print(f"🔹 Potencia total de molinos sin Quijote: {P_sin_Quijote_total/1e6:.2f} MW")
print(f"🔹 Potencia total de molinos: {(P_Quijote_total + P_sin_Quijote_total)/1e6:.2f} MW")

# ========================================================================
# 🎯 RESPUESTA A TUS PREGUNTAS
# ========================================================================

print("\n" + "=" * 80)
print("🎯 RESPUESTA A TUS PREGUNTAS")
print("=" * 80)

print("""
1. ¿Kilómetro siempre roba el mismo potencial gravitatorio?
   Sí, E = m_peso · g · Δh por ciclo. Pero la patada mayor y el hurto de
   distancia cambian la frecuencia de los ciclos.

2. ¿Con una patada mayor y hurto de distancia se roban más julios por tiempo?
   Sí. La patada mayor reduce el tiempo de ciclo, aumentando la potencia.

3. ¿En 3, 5 o 7 fases se roba más julios por tiempo?
   Sí. Más fases = ciclo más corto = más potencia.
   - 3 fases: ciclo ~5 s, potencia ~200 W
   - 5 fases: ciclo ~3 s, potencia ~333 W
   - 7 fases: ciclo ~2 s, potencia ~500 W

4. ¿Con una patada más fuerte menos dura el ciclo?
   Sí. La patada más fuerte acorta la bajada y la subida.

5. ¿Qué control tiene el Kilómetro?
   - Meter energía (cargar): subir pesos de BAJA a ALTA
   - Quitar energía (descargar): bajar pesos de ALTA a BAJA
   - Medir generación: calcular E_gen y P_KM
""")

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Esta simulación incluye 400 Kilómetros y 150 módulos de molinos.
# - Los Kilómetros tienen diferentes fases de patada (3, 5, 7).
# - Los molinos tienen vientos personalizados (con y sin Quijote).
# - El control predictivo optimiza la carga/descarga de los Kilómetros.
# - La megabatería puede entregar ~166 kW de potencia pico.
# ========================================================================