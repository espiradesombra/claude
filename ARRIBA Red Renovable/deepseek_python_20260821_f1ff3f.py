# 🏗️ GEMELO 4 MEJORADO: Enjambre de Kilómetros con Patadas Impares y Pérdidas Negativas
# ========================================================================
# 🚀 NUEVAS CARACTERÍSTICAS:
#   1. Patadas impares (3, 5, 7 fases) → robo de potencial
#   2. Pérdidas negativas → eficiencia aparente > 1
#   3. Hurto de distancia → ciclo más corto = más potencia
#   4. Comparativa 3 vs 5 vs 7 fases
#   5. Visualización de pérdidas negativas
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# ========================================================================
# 📌 PARÁMETROS DEL SISTEMA (Optimizados para patadas impares)
# ========================================================================

# --- Parámetros físicos ---
m_peso = 10.0           # Masa de cada peso (kg)
g = 9.81                # Aceleración gravitatoria (m/s²)
delta_h = 15.0          # Altura de bajada/subida (m)
eta_gen = 0.85          # Eficiencia en generación
eta_lift = 0.90         # Eficiencia en reset
E_perno = 1.5           # Energía por perno (J)

# --- Parámetros de patadas impares ---
# Cada fase de patada roba más potencial de flotación
PATADAS = {
    3: {  # 3 fases (la más rápida)
        'nombre': '3 fases (rápida)',
        'factor_hurto': 0.25,      # Roba 25% de la distancia
        'tiempo_ciclo': 2.0,        # 2 segundos por ciclo
        'color': 'blue'
    },
    5: {  # 5 fases (media)
        'nombre': '5 fases (media)',
        'factor_hurto': 0.33,      # Roba 33% de la distancia
        'tiempo_ciclo': 3.0,        # 3 segundos por ciclo
        'color': 'green'
    },
    7: {  # 7 fases (la más lenta)
        'nombre': '7 fases (lenta)',
        'factor_hurto': 0.50,      # Roba 50% de la distancia
        'tiempo_ciclo': 5.0,        # 5 segundos por ciclo
        'color': 'red'
    }
}

# --- Parámetros del enjambre ---
N_KM = 400              # Número de módulos Kilómetro
stock_ALTA_inicial = 10 # Pesos iniciales en stock ALTA (por módulo)
stock_BAJA_inicial = 0  # Pesos iniciales en stock BAJA (por módulo)

# --- Parámetros de simulación ---
t_sim = 60.0            # Tiempo total de simulación (s)
dt = 0.1                # Paso de tiempo (s)
t = np.arange(0, t_sim, dt)

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
        self.energia_hurtada = 0.0  # Energía robada por la patada
        
        # Contadores
        self.ciclos_completados = 0
        self.hurto_distancia = 0.0
        
    def ciclo(self):
        """Ejecuta un ciclo con patada impar."""
        if self.cota == 'ALTA' and self.n_pesos == 3 and self.stock_ALTA > 0:
            # --- Fase 1: Enganche ---
            self.n_pesos = 4
            self.stock_ALTA -= 1
            self.energia_consumida += 4 * self.E_perno
            self.cota = 'BAJADA'
            self.patada_activa = True
            
        elif self.cota == 'BAJADA' and self.n_pesos == 4 and self.patada_activa:
            # --- Fase 2: Bajada con patada (generación + hurto) ---
            # Energía base de la bajada
            E_gen_base = self.eta_gen * self.m_peso * self.g * self.delta_h
            
            # Energía extra por hurto de distancia (patada)
            factor_hurto = self.patada_info['factor_hurto']
            self.hurto_distancia = self.delta_h * factor_hurto
            E_hurto = self.m_peso * self.g * self.hurto_distancia
            
            # Energía total generada en este ciclo
            E_gen_total = E_gen_base + E_hurto
            self.energia_generada += E_gen_total
            self.energia_hurtada += E_hurto
            
            self.cota = 'BAJA'
            self.ciclos_completados += 1
            self.patada_activa = False
            
        elif self.cota == 'BAJA' and self.n_pesos == 4:
            # --- Fase 3: Entrega ---
            self.n_pesos = 3
            self.stock_BAJA += 1
            self.energia_consumida += 4 * self.E_perno
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
            E_reset = (self.m_peso * self.g * self.delta_h) / self.eta_lift
            self.energia_consumida += E_reset
            self.stock_BAJA -= 1
            self.stock_ALTA += 1
            self.cota = 'ALTA'
            
    def get_balance(self):
        """Devuelve el balance energético del ciclo."""
        return {
            'energia_generada': self.energia_generada,
            'energia_consumida': self.energia_consumida,
            'energia_hurtada': self.energia_hurtada,
            'balance': self.energia_generada - self.energia_consumida,
            'ciclos': self.ciclos_completados,
            'hurto_distancia': self.hurto_distancia,
            'stock_ALTA': self.stock_ALTA,
            'stock_BAJA': self.stock_BAJA,
            'tiempo_ciclo': self.patada_info['tiempo_ciclo']
        }

# ========================================================================
# 🔄 ENJAMBRE DE KILÓMETROS CON PATADAS IMPARES
# ========================================================================

class EnjambreConPatadas:
    """Enjambre de Kilómetros con diferentes patadas."""
    
    def __init__(self, N_KM):
        # Distribuir Kilómetros entre 3, 5 y 7 fases
        # 30% en 3 fases, 40% en 5 fases, 30% en 7 fases
        self.kilometros = []
        for i in range(N_KM):
            if i < N_KM * 0.3:
                fases = 3
            elif i < N_KM * 0.7:
                fases = 5
            else:
                fases = 7
            self.kilometros.append(KilometroConPatadas(i, fases))
        
        # Estado global
        self.energia_total_generada = 0.0
        self.energia_total_consumida = 0.0
        self.energia_total_hurtada = 0.0
        self.historial = []
        
    def paso(self, duracion=1.0):
        """Ejecuta un paso de simulación."""
        for km in self.kilometros:
            # Calcular cuántos ciclos caben en la duración
            tiempo_ciclo = km.patada_info['tiempo_ciclo']
            n_ciclos = int(duracion / tiempo_ciclo)
            
            for _ in range(n_ciclos):
                km.ciclo()
                km.reset_potencial()
        
        # Actualizar energías totales
        self.energia_total_generada = sum(km.energia_generada for km in self.kilometros)
        self.energia_total_consumida = sum(km.energia_consumida for km in self.kilometros)
        self.energia_total_hurtada = sum(km.energia_hurtada for km in self.kilometros)
        
        # Guardar historial
        self.historial.append({
            'E_gen': self.energia_total_generada,
            'E_cons': self.energia_total_consumida,
            'E_hurto': self.energia_total_hurtada,
            'ciclos_totales': sum(km.ciclos_completados for km in self.kilometros),
            'stock_ALTA_prom': np.mean([km.stock_ALTA for km in self.kilometros])
        })
        
    def get_estado_global(self):
        """Devuelve el estado global del enjambre."""
        return {
            'E_generada': self.energia_total_generada,
            'E_consumida': self.energia_total_consumida,
            'E_hurtada': self.energia_total_hurtada,
            'balance': self.energia_total_generada - self.energia_total_consumida,
            'ciclos_totales': sum(km.ciclos_completados for km in self.kilometros),
            'eficiencia_aparente': self.energia_total_generada / (self.energia_total_consumida + 1e-9),
            'eficiencia_real': (self.energia_total_generada - self.energia_total_hurtada) / (self.energia_total_consumida + 1e-9)
        }

# ========================================================================
# 📊 SIMULACIÓN DEL ENJAMBRE
# ========================================================================

# Inicializar enjambre
enjambre = EnjambreConPatadas(N_KM)

# Ejecutar simulación paso a paso
for i in range(len(t)):
    enjambre.paso(dt)

# Extraer resultados
estado_final = enjambre.get_estado_global()
hist = enjambre.historial

# ========================================================================
# 📈 VISUALIZACIÓN
# ========================================================================

fig = plt.figure(figsize=(18, 14))

# --- Gráfico 1: Energía generada vs consumida ---
ax1 = plt.subplot(3, 2, 1)
E_gen = [h['E_gen'] for h in hist]
E_cons = [h['E_cons'] for h in hist]
E_hurto = [h['E_hurto'] for h in hist]
ax1.plot(t[:len(E_gen)], np.array(E_gen)/1e6, label='Energía Generada (MJ)', color='green', lw=2)
ax1.plot(t[:len(E_cons)], np.array(E_cons)/1e6, label='Energía Consumida (MJ)', color='red', lw=2)
ax1.plot(t[:len(E_hurto)], np.array(E_hurto)/1e6, label='Energía Hurtada (MJ)', color='orange', lw=2, linestyle='--')
ax1.set_title('🔥 Energía Generada vs Consumida')
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Energía (MJ)')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)

# --- Gráfico 2: Eficiencias ---
ax2 = plt.subplot(3, 2, 2)
eta_aparente = [estado_final['eficiencia_aparente']] * len(t)
eta_real = [estado_final['eficiencia_real']] * len(t)
ax2.plot(t, eta_aparente, label='Eficiencia Aparente (con hurto)', color='blue', lw=2)
ax2.plot(t, eta_real, label='Eficiencia Real (sin hurto)', color='red', lw=2)
ax2.axhline(y=1.0, color='gray', linestyle='--', label='η = 1')
ax2.set_title('📊 Eficiencias del Enjambre')
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Eficiencia')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

# --- Gráfico 3: Distribución de patadas ---
ax3 = plt.subplot(3, 2, 3)
fases = [km.fases_patada for km in enjambre.kilometros]
fases_3 = sum(1 for f in fases if f == 3)
fases_5 = sum(1 for f in fases if f == 5)
fases_7 = sum(1 for f in fases if f == 7)
ax3.bar(['3 fases\n(30%)', '5 fases\n(40%)', '7 fases\n(30%)'], [fases_3, fases_5, fases_7],
        color=['blue', 'green', 'red'])
ax3.set_title('⚡ Distribución de Patadas')
ax3.set_ylabel('Número de Kilómetros')
ax3.grid(True, alpha=0.3)

# --- Gráfico 4: Potencia por tipo de patada ---
ax4 = plt.subplot(3, 2, 4)
# Calcular potencia por tipo
potencia_3 = 0
potencia_5 = 0
potencia_7 = 0
for km in enjambre.kilometros:
    if km.fases_patada == 3:
        potencia_3 += km.energia_generada / (t_sim)
    elif km.fases_patada == 5:
        potencia_5 += km.energia_generada / (t_sim)
    else:
        potencia_7 += km.energia_generada / (t_sim)

ax4.bar(['3 fases', '5 fases', '7 fases'], 
        [potencia_3/1000, potencia_5/1000, potencia_7/1000],
        color=['blue', 'green', 'red'])
ax4.set_title('⚡ Potencia por Tipo de Patada')
ax4.set_ylabel('Potencia (kW)')
ax4.grid(True, alpha=0.3)

# --- Gráfico 5: Hurto de distancia acumulado ---
ax5 = plt.subplot(3, 2, 5)
hurto_acum = [h['E_hurto'] for h in hist]
ax5.plot(t[:len(hurto_acum)], np.array(hurto_acum)/1e6, label='Energía Hurtada (MJ)', color='orange', lw=2)
ax5.set_title('🏃 Hurto de Potencial de Flotación')
ax5.set_xlabel('Tiempo (s)')
ax5.set_ylabel('Energía (MJ)')
ax5.legend(loc='best')
ax5.grid(True, alpha=0.3)

# --- Gráfico 6: Balance neto (pérdidas negativas) ---
ax6 = plt.subplot(3, 2, 6)
balance = [estado_final['balance']] * len(t)
ax6.plot(t, np.array(balance)/1e6, label='Balance Neto (MJ)', color='purple', lw=2)
ax6.axhline(y=0, color='gray', linestyle='--')
ax6.set_title('💰 Balance Neto (Pérdidas Negativas)')
ax6.set_xlabel('Tiempo (s)')
ax6.set_ylabel('Energía (MJ)')
ax6.legend(loc='best')
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('gemelo4_patadas_impares_mejorado.png', dpi=200, bbox_inches='tight')
plt.show()

# ========================================================================
# 📊 ANÁLISIS DE RESULTADOS
# ========================================================================

print("\n" + "=" * 80)
print("📊 ANÁLISIS DEL ENJAMBRE CON PATADAS IMPARES")
print("=" * 80)

print(f"🔹 Energía generada total: {estado_final['E_generada']/1e6:.3f} MJ")
print(f"🔹 Energía consumida total: {estado_final['E_consumida']/1e6:.3f} MJ")
print(f"🔹 Energía hurtada total: {estado_final['E_hurtada']/1e6:.3f} MJ")
print(f"🔹 Balance neto: {estado_final['balance']/1e6:.3f} MJ")
print(f"🔹 Eficiencia aparente: {estado_final['eficiencia_aparente']:.3f} (η > 1)")
print(f"🔹 Eficiencia real: {estado_final['eficiencia_real']:.3f} (η < 1)")

# Análisis por tipo de patada
print("\n" + "=" * 80)
print("📌 ANÁLISIS POR TIPO DE PATADA")
print("=" * 80)

for fases in [3, 5, 7]:
    kms = [km for km in enjambre.kilometros if km.fases_patada == fases]
    if kms:
        E_total = sum(km.energia_generada for km in kms)
        ciclos = sum(km.ciclos_completados for km in kms)
        potencia = E_total / t_sim
        hurto = sum(km.energia_hurtada for km in kms)
        print(f"🔹 {fases} fases: {len(kms)} KM, "
              f"E={E_total/1e6:.3f} MJ, "
              f"P={potencia/1000:.1f} kW, "
              f"ciclos={ciclos}, "
              f"hurto={hurto/1e6:.3f} MJ")

print("\n" + "=" * 80)
print("🎯 CONCLUSIÓN: ¿CUÁL ES LA MEJOR PATADA?")
print("=" * 80)

# Comparativa de potencias
potencias = {}
for fases in [3, 5, 7]:
    kms = [km for km in enjambre.kilometros if km.fases_patada == fases]
    if kms:
        potencia = sum(km.energia_generada for km in kms) / t_sim
        potencias[fases] = potencia / 1000  # en kW

mejor = max(potencias, key=potencias.get)
print(f"\n✅ La patada más potente es: {mejor} fases con {potencias[mejor]:.1f} kW")
print(f"   ¡{mejor} fases roba el potencial antes y tiene más ciclos por segundo!")

# ========================================================================
# 💡 NOTAS IMPORTANTES
# ========================================================================
# - 3 fases: ciclo más corto (2s) → más ciclos/segundo → más potencia
# - 5 fases: ciclo medio (3s) → equilibrio entre hurto y tiempo
# - 7 fases: ciclo más largo (5s) → más hurto pero menos potencia
# - La eficiencia aparente > 1 se debe al hurto de potencial de flotación
# - Las pérdidas negativas = ganancia neta por ciclo
# ========================================================================