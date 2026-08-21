# FILE: gemelo_1_dac.py
# Objetivo: Simular la captura de CO2 del escape del motor.
# Modelo: Sistema de captura directa (DAC) con zeolita.
# Entrada: Consumo de combustible. Salida: CO2 capturado, almacenado y emitido.

import numpy as np
import matplotlib.pyplot as plt
import json

# --- PARÁMETROS DEL MOTOR ---
CONSUMO_L_100KM = 5.0
CO2_KG_POR_LITRO = 3.15
EFICIENCIA_CAPTURA_DAC = 0.88  # 88% de eficiencia
PRESION_ALMACENAMIENTO_BAR = 150

# --- PARÁMETROS DEL DEPÓSITO DE CO2 ---
VOLUMEN_DEPOSITO_LITROS = 15.0
DENSIDAD_CO2_A_150BAR = 500  # kg/m3 (aprox.)
CAPACIDAD_MAXIMA_KG = (VOLUMEN_DEPOSITO_LITROS / 1000) * DENSIDAD_CO2_A_150BAR

# --- SIMULACIÓN DE UNA JORNADA DE CONDUCCIÓN (8 HORAS) ---
def simular_captura_co2(km_recorridos, duracion_horas):
    """Simula la captura de CO2 durante una jornada de conducción."""
    
    # 1. Cálculo de Consumo y Emisiones Totales
    litros_consumidos = (CONSUMO_L_100KM / 100) * km_recorridos
    co2_total_producido_kg = litros_consumidos * CO2_KG_POR_LITRO
    
    # 2. Captura y Emisión
    co2_capturado_kg = co2_total_producido_kg * EFICIENCIA_CAPTURA_DAC
    co2_emitido_kg = co2_total_producido_kg - co2_capturado_kg
    
    # 3. Almacenamiento en el Depósito
    co2_almacenado_kg = min(co2_capturado_kg, CAPACIDAD_MAXIMA_KG)
    co2_desbordado_kg = co2_capturado_kg - co2_almacenado_kg
    
    # El CO2 desbordado se emite a la atmósfera
    co2_emitido_final_kg = co2_emitido_kg + co2_desbordado_kg
    
    # 4. Energía Requerida para el DAC
    # El calor residual del motor cubre la energía de activación
    energia_dac_kwh = 0.35 * duracion_horas  # 0.35 kW de consumo medio
    
    # 5. Resultados
    resultados = {
        "km_recorridos": km_recorridos,
        "litros_consumidos": litros_consumidos,
        "co2_total_producido_kg": co2_total_producido_kg,
        "co2_capturado_kg": co2_capturado_kg,
        "co2_almacenado_deposito_kg": co2_almacenado_kg,
        "co2_emitido_final_kg": co2_emitido_final_kg,
        "eficiencia_captura": EFICIENCIA_CAPTURA_DAC,
        "capacidad_deposito_kg": CAPACIDAD_MAXIMA_KG,
        "porcentaje_deposito_lleno": (co2_almacenado_kg / CAPACIDAD_MAXIMA_KG) * 100,
        "energia_dac_kwh": energia_dac_kwh,
        "reduccion_co2_pct": (1 - (co2_emitido_final_kg / co2_total_producido_kg)) * 100,
    }
    
    return resultados

# --- EJECUTAR SIMULACIÓN ---
km_totales = 400  # Ejemplo: 400 km en un día
horas_conduccion = 8

resultados = simular_captura_co2(km_totales, horas_conduccion)

# --- IMPRIMIR Y GUARDAR RESULTADOS ---
print("="*60)
print("RESULTADOS DE LA SIMULACIÓN DE CAPTURA DE CO2 (DAC)")
print("="*60)
print(f"Km recorridos:                 {resultados['km_recorridos']:.0f} km")
print(f"Litros de gasoil consumidos:   {resultados['litros_consumidos']:.2f} L")
print(f"CO2 total producido:           {resultados['co2_total_producido_kg']:.2f} kg")
print(f"CO2 capturado por el DAC:      {resultados['co2_capturado_kg']:.2f} kg ({resultados['eficiencia_captura']*100:.0f}%)")
print(f"CO2 almacenado en depósito:    {resultados['co2_almacenado_deposito_kg']:.2f} kg")
print(f"Capacidad del depósito:        {resultados['capacidad_deposito_kg']:.2f} kg")
print(f"Depósito lleno al:             {resultados['porcentaje_deposito_lleno']:.1f}%")
print(f"CO2 emitido final:             {resultados['co2_emitido_final_kg']:.2f} kg")
print(f"REDUCCIÓN DE EMISIONES:        {resultados['reduccion_co2_pct']:.1f}%")
print(f"Energía consumida por DAC:     {resultados['energia_dac_kwh']:.2f} kWh")
print("="*60)

# Guardar resultados en JSON
with open('gemelo_1_dac_resultados.json', 'w') as f:
    json.dump(resultados, f, indent=4)

# --- GRÁFICA ---
categorias = ['CO2 Producido', 'CO2 Capturado', 'CO2 Almacenado', 'CO2 Emitido']
valores = [
    resultados['co2_total_producido_kg'],
    resultados['co2_capturado_kg'],
    resultados['co2_almacenado_deposito_kg'],
    resultados['co2_emitido_final_kg']
]

plt.figure(figsize=(10, 6))
plt.bar(categorias, valores, color=['red', 'green', 'blue', 'orange'])
plt.ylabel('Masa de CO2 (kg)')
plt.title('Balance de CO2 en el Sistema de Captura (DAC)')
plt.grid(axis='y', linestyle='--', alpha=0.7)

for i, v in enumerate(valores):
    plt.text(i, v + 0.5, f"{v:.1f} kg", ha='center', va='bottom')

plt.tight_layout()
plt.savefig('gemelo_1_dac_grafico.png')
# plt.show()
print("Gráfico guardado como 'gemelo_1_dac_grafico.png'")