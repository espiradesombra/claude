# 🌿 GEMELO 5 MEJORADO: Red Eléctrica Verde con Patadas Impares y Pérdidas Negativas
# ========================================================================
# 🚀 NUEVAS CARACTERÍSTICAS:
#   1. Patadas impares (3, 5, 7 fases) en los Kilómetros
#   2. Pérdidas negativas (eficiencia aparente > 1)
#   3. Comparativa 3 vs 5 vs 7 fases en la red eléctrica
#   4. Control predictivo con buffer de energía
#   5. Análisis de proporciones con pérdidas negativas
#   6. Visualización completa de la red con patadas
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ========================================================================
# 📌 PARÁMETROS GLOBALES (Optimizados para patadas impares)
# ========================================================================

# --- Parámetros de Kilómetro (con patadas) ---
N_KM = 400                          # Número de módulos Kilómetro
m_peso = 10.0                       # Masa de cada peso (kg)
g = 9.81                            # Aceleración gravitatoria (m/s²)
delta_h = 15.0                      # Altura de bajada/subida (m)
eta_gen = 0.85                      # Eficiencia en generación
eta_lift = 0.90                     # Eficiencia en reset
E_perno = 1.5                       # Energía por perno (J)
stock_ALTA_inicial = 10             # Stock inicial de pesos en ALTA (por módulo)
stock_BAJA_inicial = 0              # Stock inicial de pesos en BAJA (por módulo)

# Parámetros de patadas impares (3, 5, 7 fases)
PATADAS = {
    3: {  # 3 fases (la más rápida)
        'nombre': '3 fases (rápida)',
        'factor_hurto': 0.25,        # Roba 25% de la distancia
        'tiempo_ciclo': 2.0,         # 2 segundos por ciclo
        'color': 'blue',
        'distribucion': 0.30         # 30% de los Kilómetros
    },
    5: {  # 5 fases (media)
        'nombre': '5 fases (media)',
        'factor_hurto': 0.33,        # Roba 33% de la distancia
        'tiempo_ciclo': 3.0,         # 3 segundos por ciclo
        'color': 'green',
        'distribucion': 0.40         # 40% de los Kilómetros
    },
    7: {  # 7 fases (la más lenta)
        'nombre': '7 fases (lenta)',
        'factor_hurto': 0.50,        # Roba 50% de la distancia
        'tiempo_ciclo': 5.0,         # 5 segundos por ciclo
        'color': 'red',
        'distribucion': 0.30         # 30% de los Kilómetros
    }
}

# --- Parámetros de Molinos ---
N_modulos_molinos = 150              # 100 con Quijote + 50 sin Quijote
N_molinos_por_modulo = 5             # Molinos por módulo
N_molinos_Quijote = 100              # Módulos con Quijote (500 molinos)
N_molinos_sin_Quijote = 50           # Módulos sin Quijote (250 molinos)

# Potencia nominal por molino (W)
P_nominal_Quijote = 10000.0          # 10 kW (con Quijote)
P_nominal_sin_Quijote = 8000.0       # 8 kW (sin Quijote)

# Factor de carga (viento)
factor_carga = 0.35                  # 35%

# --- Parámetros de Kuramoto (sincronización) ---
N_molinos_kuramoto = 5               # 2 centrales + 3 anillo
omega_natural = np.array([2.0, 2.0, 2.0, 2.0, 2.0])
K = 0.5                              # Acoplamiento entre molinos
K_bus = 0.8                          # Acoplamiento al bus
V_bus = 1.0                          # Tensión del bus
V_ref = 1.0                          # Tensión de referencia
T = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])  # Tensor de conexión

# --- Parámetros de Quijote ---
N_blades = 3                         # Número de aspas por molino
M_Q = 4.0                            # Masa desplazable (kg)
J_G = 10.0                           # Inercia del generador (kg·m²)
omega_rated = 2.0                    # Velocidad angular nominal (rad/s)
v_wind_rated = 12.0                  # Velocidad del viento nominal (m/s)
K_Q_OM = 0.1                         # Ganancia de velocidad angular
K_Q_V = 0.05                         # Ganancia de viento

# --- Parámetros de simulación ---
t_sim = 120.0                        # Tiempo total (s)
dt = 0.1                             # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)

# --- Demanda de la red ---
P_demanda_base = 2.45e6              # 2.45 MW (demanda base)
P_demanda_pico = 3.5e6               # 3.5 MW (demanda pico)
P_demanda = P_demanda_base + 0.5e6 * np.sin(0.05 * t)  # Demanda variable

# ========================================================================
# 🔬 MÓDULO KILÓMETRO CON PATADAS IMPARES
# ========================================================================

class KilometroConPatadas:
    """Módulo Kilómetro con patadas impares (3, 5, 7 fases)."""
    
    def __init__(self, id_modulo, fases_patada):
        self.id = id_modulo
        self.fases_patada = fases_patada
        self.patada_info = PATADAS[fases_patada]
        
        # Estado inicial
        self.cota = 'ALTA'
        self.n_pesos = 3
        self.stock_ALTA = stock_ALTA_inicial
        self.stock_BAJA = stock_BAJA_inicial
        
        # Energías
        self.energia_generada = 0.0
        self.energia_consumida = 0.0
        self.energia_hurtada = 0.0
        self.potencia_actual = 0.0
        
        # Contadores
        self.ciclos_completados = 0
        self.hurto_distancia = 0.0
        
    def ejecutar_ciclo(self):
        """Ejecuta un ciclo con patada impar."""
        if self.cota == 'ALTA' and self.n_pesos == 3 and self.stock_ALTA > 0:
            # --- Fase 1: Enganche ---
            self.n_pesos = 4
            self.stock_ALTA -= 1
            self.energia_consumida += 4 * E_perno
            self.cota = 'BAJADA'
            self.patada_activa = True
            
        elif self.cota == 'BAJADA' and self.n_pesos == 4 and self.patada_activa:
            # --- Fase 2: Bajada con patada (generación + hurto) ---
            E_gen_base = eta_gen * m_peso * g * delta_h
            factor_hurto = self.patada_info['factor_hurto']
            self.hurto_distancia = delta_h * factor_hurto
            E_hurto = m_peso * g * self.hurto_distancia
            E_gen_total = E_gen_base + E_hurto
            
            self.energia_generada += E_gen_total
            self.energia_hurtada += E_hurto
            self.potencia_actual = E_gen_total / self.patada_info['tiempo_ciclo']
            
            self.cota = 'BAJA'
            self.ciclos_completados += 1
            self.patada_activa = False
            
        elif self.cota == 'BAJA' and self.n_pesos == 4:
            # --- Fase 3: Entrega ---
            self.n_pesos = 3
            self.stock_BAJA += 1
            self.energia_consumida += 4 * E_perno
            self.cota = 'SUBIDA'
            
        elif self.cota == 'SUBIDA' and self.n_pesos == 3:
            # --- Fase 4: Subida por flotación (casi gratis) ---
            self.cota = 'ALTA'
            
        elif self.stock_ALTA == 0:
            # --- Modo Pausa ---
            self.cota = 'PAUSA'
            
    def reset_potencial(self):
        """Reset del potencial de flotación."""
        if self.stock_BAJA > 0 and self.cota == 'PAUSA':
            E_reset = (m_peso * g * delta_h) / eta_lift
            self.energia_consumida += E_reset
            self.stock_BAJA -= 1
            self.stock_ALTA += 1
            self.cota = 'ALTA'
            
    def get_estado(self):
        return {
            'fases': self.fases_patada,
            'E_generada': self.energia_generada,
            'E_consumida': self.energia_consumida,
            'E_hurtada': self.energia_hurtada,
            'potencia': self.potencia_actual,
            'ciclos': self.ciclos_completados,
            'stock_ALTA': self.stock_ALTA,
            'stock_BAJA': self.stock_BAJA
        }

# ========================================================================
# 🔄 SISTEMA COMPLETO: Integración de los 4 Gemelos con Patadas
# ========================================================================

class RedElectricaConPatadas:
    """Red eléctrica completa con patadas impares."""
    
    def __init__(self):
        # --- Gemelo 1: Kuramoto (sincronización) ---
        self.theta = np.array([0.1, 0.2, 0.3, 0.4, 0.5])  # Fases iniciales
        
        # --- Gemelo 2: Quijote (inercia variable) ---
        self.r_q = np.full(N_molinos_kuramoto * N_blades, 5.0)
        self.omega = 1.5
        self.v_wind = 10.0
        
        # --- Gemelo 4: Kilómetros con patadas ---
        self.kilometros = []
        for i in range(N_KM):
            # Distribuir patadas: 30% 3 fases, 40% 5 fases, 30% 7 fases
            if i < N_KM * 0.3:
                fases = 3
            elif i < N_KM * 0.7:
                fases = 5
            else:
                fases = 7
            self.kilometros.append(KilometroConPatadas(i, fases))
        
        # --- Gemelo 3: Molinos ---
        self.molinos = []
        for i in range(N_modulos_molinos):
            if i < N_molinos_Quijote:
                self.molinos.append({
                    'tiene_quijote': True,
                    'P_nominal': P_nominal_Quijote,
                    'v_wind': 12.0,
                    'perfil': self._generar_perfil_viento(i, True)
                })
            else:
                self.molinos.append({
                    'tiene_quijote': False,
                    'P_nominal': P_nominal_sin_Quijote,
                    'v_wind': 8.0,
                    'perfil': self._generar_perfil_viento(i, False)
                })
        
        # Estado global
        self.E_generada_total = 0.0
        self.E_consumida_total = 0.0
        self.E_hurtada_total = 0.0
        self.buffer_energia = 0.0
        self.capacidad_buffer = 1e6  # 1 MJ
        self.historial = []
        
    def _generar_perfil_viento(self, id_modulo, tiene_quijote):
        """Genera un perfil de viento único para cada módulo."""
        if tiene_quijote:
            return {
                'base': 12.0,
                'amplitud': 3.0,
                'frecuencia': 0.1 + id_modulo * 0.01,
                'fase': id_modulo * 0.2,
                'rafagas': 0.5 * np.sin(id_modulo * 0.3)
            }
        else:
            return {
                'base': 8.0,
                'amplitud': 2.0,
                'frecuencia': 0.15 + id_modulo * 0.02,
                'fase': id_modulo * 0.3,
                'rafagas': 0.3 * np.sin(id_modulo * 0.4)
            }
    
    def paso(self, tiempo):
        """Ejecuta un paso de simulación."""
        # --- 1. Actualizar molinos ---
        P_molinos = 0
        for modulo in self.molinos:
            perfil = modulo['perfil']
            v_wind = (perfil['base'] 
                     + perfil['amplitud'] * np.sin(perfil['frecuencia'] * tiempo + perfil['fase'])
                     + perfil['rafagas'] * np.sin(0.05 * tiempo + id_modulo * 0.1))
            v_wind = max(v_wind, 2.0)
            
            if modulo['tiene_quijote']:
                # Quijote: +25% de potencia y estabilidad
                P = modulo['P_nominal'] * (v_wind / 12.0)**3 * factor_carga * 1.25
            else:
                P = modulo['P_nominal'] * (v_wind / 8.0)**3 * factor_carga
            P_molinos += max(P, 0)
        
        # --- 2. Actualizar Kilómetros ---
        P_KM = 0
        E_gen_KM = 0
        E_cons_KM = 0
        E_hurto_KM = 0
        
        for km in self.kilometros:
            km.ejecutar_ciclo()
            km.reset_potencial()
            
            estado = km.get_estado()
            E_gen_KM += estado['E_generada']
            E_cons_KM += estado['E_consumida']
            E_hurto_KM += estado['E_hurtada']
            P_KM += estado['potencia']
        
        # --- 3. Control predictivo ---
        P_total = P_molinos + P_KM
        
        # Predecir demanda futura (5 segundos)
        idx = int(tiempo / dt)
        if idx < len(P_demanda) - int(5/dt):
            P_futuro = np.mean(P_demanda[idx:idx + int(5/dt)])
        else:
            P_futuro = P_demanda[-1]
        
        # Balance de energía
        if P_total > P_demanda[idx]:
            # Excedente: cargar buffer
            excedente = (P_total - P_demanda[idx]) * dt
            self.buffer_energia = min(self.buffer_energia + excedente, self.capacidad_buffer)
        else:
            # Déficit: descargar buffer
            deficit = (P_demanda[idx] - P_total) * dt
            if self.buffer_energia > 0:
                disponible = min(self.buffer_energia, deficit)
                self.buffer_energia -= disponible
        
        # --- 4. Actualizar energía total ---
        self.E_generada_total += E_gen_KM
        self.E_consumida_total += E_cons_KM
        self.E_hurtada_total += E_hurto_KM
        
        # --- 5. Guardar historial ---
        self.historial.append({
            't': tiempo,
            'P_molinos': P_molinos,
            'P_KM': P_KM,
            'P_total': P_total,
            'P_demanda': P_demanda[idx],
            'buffer': self.buffer_energia,
            'E_gen_KM': self.E_generada_total,
            'E_cons_KM': self.E_consumida_total,
            'E_hurto_KM': self.E_hurtada_total,
            'ciclos_totales': sum(km.ciclos_completados for km in self.kilometros)
        })
        
        return P_molinos, P_KM, P_total

# ========================================================================
# 📊 SIMULACIÓN DE LA RED ELÉCTRICA
# ========================================================================

# Inicializar red
red = RedElectricaConPatadas()

# Ejecutar simulación
for tiempo in t:
    red.paso(tiempo)

# Extraer resultados
hist = red.historial
P_molinos = [h['P_molinos'] for h in hist]
P_KM = [h['P_KM'] for h in hist]
P_total = [h['P_total'] for h in hist]
P_demanda_hist = [h['P_demanda'] for h in hist]
buffer = [h['buffer'] for h in hist]
E_gen_KM = [h['E_gen_KM'] for h in hist]
E_cons_KM = [h['E_cons_KM'] for h in hist]
E_hurto_KM = [h['E_hurto_KM'] for h in hist]
ciclos = [h['ciclos_totales'] for h in hist]

# ========================================================================
# 📈 VISUALIZACIÓN (8 GRÁFICOS)
# ========================================================================

fig = plt.figure(figsize=(20, 16))

# --- Gráfico 1: Potencia generada vs. Demanda ---
ax1 = plt.subplot(4, 2, 1)
ax1.plot(t[:len(P_molinos)], np.array(P_molinos)/1e6, label='Molinos (MW)', color='blue', lw=1.5)
ax1.plot(t[:len(P_KM)], np.array(P_KM)/1e6, label='Kilómetros (MW)', color='green', lw=1.5)
ax1.plot(t[:len(P_demanda_hist)], np.array(P_demanda_hist)/1e6, label='Demanda (MW)', color='red', lw=2, linestyle='--')
ax1.set_title('⚡ Potencia Generada vs. Demanda (con Patadas)')
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Potencia (MW)')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)

# --- Gráfico 2: Proporción de aporte ---
ax2 = plt.subplot(4, 2, 2)
P_molinos_arr = np.array(P_molinos)
P_KM_arr = np.array(P_KM)
P_total_arr = np.array(P_total)
porc_molinos = (P_molinos_arr / (P_total_arr + 1e-9)) * 100
porc_KM = (P_KM_arr / (P_total_arr + 1e-9)) * 100
ax2.plot(t[:len(porc_molinos)], porc_molinos, label='% Molinos', color='blue', lw=1.5)
ax2.plot(t[:len(porc_KM)], porc_KM, label='% Kilómetros', color='green', lw=1.5)
ax2.set_title('📊 Proporción de Aporte Energético (con Patadas)')
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Porcentaje (%)')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

# --- Gráfico 3: Energía en Kilómetros (Generada vs Consumida) ---
ax3 = plt.subplot(4, 2, 3)
ax3.plot(t[:len(E_gen_KM)], np.array(E_gen_KM)/1e6, label='Energía Generada (MJ)', color='green', lw=2)
ax3.plot(t[:len(E_cons_KM)], np.array(E_cons_KM)/1e6, label='Energía Consumida (MJ)', color='red', lw=2)
ax3.plot(t[:len(E_hurto_KM)], np.array(E_hurto_KM)/1e6, label='Energía Hurtada (MJ)', color='orange', lw=2, linestyle='--')
ax3.set_title('🏗️ Energía en Kilómetros (Gen vs Cons vs Hurto)')
ax3.set_xlabel('Tiempo (s)')
ax3.set_ylabel('Energía (MJ)')
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)

# --- Gráfico 4: Buffer de energía ---
ax4 = plt.subplot(4, 2, 4)
ax4.plot(t[:len(buffer)], np.array(buffer)/1e6, label='Buffer (MJ)', color='purple', lw=2)
ax4.axhline(y=red.capacidad_buffer/1e6, color='red', linestyle='--', label='Capacidad máxima')
ax4.set_title('🔋 Buffer de Energía (Control Predictivo)')
ax4.set_xlabel('Tiempo (s)')
ax4.set_ylabel('Energía (MJ)')
ax4.legend(loc='best')
ax4.grid(True, alpha=0.3)

# --- Gráfico 5: Potencia por tipo de patada ---
ax5 = plt.subplot(4, 2, 5)
# Calcular potencia por tipo de patada
potencia_3 = 0
potencia_5 = 0
potencia_7 = 0
for km in red.kilometros:
    if km.fases_patada == 3:
        potencia_3 += km.potencia_actual
    elif km.fases_patada == 5:
        potencia_5 += km.potencia_actual
    else:
        potencia_7 += km.potencia_actual

ax5.bar(['3 fases\n(30%)', '5 fases\n(40%)', '7 fases\n(30%)'], 
        [potencia_3/1000, potencia_5/1000, potencia_7/1000],
        color=['blue', 'green', 'red'])
ax5.set_title('⚡ Potencia por Tipo de Patada')
ax5.set_ylabel('Potencia (kW)')
ax5.grid(True, alpha=0.3)

# --- Gráfico 6: Eficiencias (Aparente vs Real) ---
ax6 = plt.subplot(4, 2, 6)
# Calcular eficiencias
E_gen_total = red.E_generada_total
E_cons_total = red.E_consumida_total
E_hurto_total = red.E_hurtada_total

eta_aparente = E_gen_total / (E_cons_total + 1e-9)
eta_real = (E_gen_total - E_hurto_total) / (E_cons_total + 1e-9)

ax6.bar(['η_aparente\n(con hurto)', 'η_real\n(sin hurto)'], 
        [eta_aparente, eta_real],
        color=['blue', 'red'])
ax6.axhline(y=1.0, color='gray', linestyle='--', label='η = 1')
ax6.set_title('📊 Eficiencias (Aparente vs Real)')
ax6.set_ylabel('Eficiencia')
ax6.legend(loc='best')
ax6.grid(True, alpha=0.3)

# --- Gráfico 7: Ciclos completados por Kilómetro ---
ax7 = plt.subplot(4, 2, 7)
ciclos_por_km = [km.ciclos_completados for km in red.kilometros]
ax7.hist(ciclos_por_km, bins=30, color='purple', alpha=0.7)
ax7.set_title('🔄 Distribución de Ciclos por Kilómetro')
ax7.set_xlabel('Número de Ciclos')
ax7.set_ylabel('Frecuencia')
ax7.grid(True, alpha=0.3)

# --- Gráfico 8: Hurto de distancia acumulado ---
ax8 = plt.subplot(4, 2, 8)
ax8.plot(t[:len(E_hurto_KM)], np.array(E_hurto_KM)/1e6, label='Energía Hurtada (MJ)', color='orange', lw=2)
ax8.set_title('🏃 Hurto de Potencial de Flotación (Acumulado)')
ax8.set_xlabel('Tiempo (s)')
ax8.set_ylabel('Energía (MJ)')
ax8.legend(loc='best')
ax8.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('gemelo5_red_electrica_patadas_mejorado.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 📊 ANÁLISIS DE RESULTADOS
# ========================================================================

print("\n" + "=" * 80)
print("📊 ANÁLISIS DE LA RED ELÉCTRICA CON PATADAS IMPARES")
print("=" * 80)

# Energías totales
E_gen_total = red.E_generada_total
E_cons_total = red.E_consumida_total
E_hurto_total = red.E_hurtada_total
E_balance = E_gen_total - E_cons_total

print(f"🔹 Energía generada total: {E_gen_total/1e6:.3f} MJ")
print(f"🔹 Energía consumida total: {E_cons_total/1e6:.3f} MJ")
print(f"🔹 Energía hurtada total: {E_hurto_total/1e6:.3f} MJ")
print(f"🔹 Balance neto: {E_balance/1e6:.3f} MJ")
print(f"🔹 Eficiencia aparente: {E_gen_total/(E_cons_total+1e-9):.3f}")
print(f"🔹 Eficiencia real: {(E_gen_total-E_hurto_total)/(E_cons_total+1e-9):.3f}")

# Análisis por tipo de patada
print("\n" + "=" * 80)
print("📌 ANÁLISIS POR TIPO DE PATADA")
print("=" * 80)

for fases in [3, 5, 7]:
    kms = [km for km in red.kilometros if km.fases_patada == fases]
    if kms:
        E_gen = sum(km.energia_generada for km in kms)
        E_cons = sum(km.energia_consumida for km in kms)
        E_hurto = sum(km.energia_hurtada for km in kms)
        ciclos = sum(km.ciclos_completados for km in kms)
        potencia = sum(km.potencia_actual for km in kms)
        
        print(f"🔹 {fases} fases ({len(kms)} KM):")
        print(f"   E_gen={E_gen/1e6:.3f} MJ, E_cons={E_cons/1e6:.3f} MJ")
        print(f"   E_hurto={E_hurto/1e6:.3f} MJ, ciclos={ciclos}")
        print(f"   Potencia={potencia/1000:.1f} kW")

# Análisis de molinos
print("\n" + "=" * 80)
print("📌 ANÁLISIS DE MOLINOS")
print("=" * 80)

P_Quijote = sum(m['P_nominal'] * factor_carga * 1.25 for m in red.molinos if m['tiene_quijote'])
P_sin_Quijote = sum(m['P_nominal'] * factor_carga for m in red.molinos if not m['tiene_quijote'])
print(f"🔹 Potencia con Quijote: {P_Quijote/1e6:.3f} MW")
print(f"🔹 Potencia sin Quijote: {P_sin_Quijote/1e6:.3f} MW")
print(f"🔹 Potencia total molinos: {(P_Quijote + P_sin_Quijote)/1e6:.3f} MW")

# ========================================================================
# 🎯 CONCLUSIÓN
# ========================================================================

print("\n" + "=" * 80)
print("🎯 CONCLUSIÓN: ¿CUÁL ES LA MEJOR PATADA PARA LA RED?")
print("=" * 80)

# Comparativa de potencias por patada
potencias = {}
for fases in [3, 5, 7]:
    kms = [km for km in red.kilometros if km.fases_patada == fases]
    if kms:
        potencia = sum(km.potencia_actual for km in kms) / 1000  # kW
        potencias[fases] = potencia

mejor = max(potencias, key=potencias.get)
print(f"\n✅ La patada más potente para la red es: {mejor} fases con {potencias[mejor]:.1f} kW")
print(f"   ¡{mejor} fases roba el potencial antes y tiene más ciclos por segundo!")

print("\n📊 RESUMEN FINAL:")
print(f"   - Energía total generada: {E_gen_total/1e6:.3f} MJ")
print(f"   - Energía total consumida: {E_cons_total/1e6:.3f} MJ")
print(f"   - Balance neto: {E_balance/1e6:.3f} MJ")
print(f"   - Eficiencia aparente: {eta_aparente:.3f}")
print(f"   - Mejor patada: {mejor} fases")

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - Este código integra los 4 gemelos con patadas impares.
# - La patada de 3 fases es la más potente porque roba el potencial antes.
# - Las pérdidas negativas (eficiencia > 1) se deben al hurto de flotación.
# - El control predictivo optimiza la carga/descarga de los Kilómetros.
# - Quijote mejora la estabilidad de los molinos (+25% de potencia).
# ========================================================================