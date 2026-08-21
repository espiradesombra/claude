# FILE: gemelo_3_ciclo_completo.py
# Objetivo: Simular el ciclo completo durante una semana de conducción.
# Integra los gemelos 1 y 2.

import json
from gemelo_1_dac import simular_captura_co2
from gemelo_2_sintesis import sintetizar_gasoil

# --- PARÁMETROS DE CONDUCCIÓN SEMANAL ---
DIAS_LABORALES = 5
KM_POR_DIA_LABORAL = 60
KM_SABADO = 100
KM_DOMINGO = 0
HORAS_CONDUCCION_POR_DIA = 8  # Estimación

# --- FUNCIÓN PRINCIPAL ---
def simular_ciclo_completo_semanal():
    """Simula una semana completa de conducción y regeneración."""
    
    # 1. Calcular kilometraje total
    km_totales = (DIAS_LABORALES * KM_POR_DIA_LABORAL) + KM_SABADO + KM_DOMINGO
    horas_totales = (DIAS_LABORALES * HORAS_CONDUCCION_POR_DIA) + 8  # Aprox.
    
    # 2. Captura de CO2 (Gemelo 1)
    co2_captura = simular_captura_co2(km_totales, horas_totales)
    co2_capturado_semana = co2_captura['co2_capturado_kg']
    co2_almacenado_semana = co2_captura['co2_almacenado_deposito_kg']
    co2_emitido_semana = co2_captura['co2_emitido_final_kg']
    
    # 3. Síntesis de Gasoil (Gemelo 2)
    # Usamos el CO2 capturado pero no almacenado (para no limitar)
    co2_a_sintetizar = co2_capturado_semana
    sintesis = sintetizar_gasoil(co2_a_sintetizar, horas_totales)
    gasoil_regenerado_litros = sintesis['gasoil_sintetizado_litros']
    co2_realmente_sintetizado = sintesis['co2_sintetizado_kg']
    
    # 4. Balance Final
    co2_emitido_final = co2_emitido_semana + (co2_capturado_semana - co2_realmente_sintetizado)
    
    consumo_litros = (CONSUMO_L_100KM / 100) * km_totales
    porcentaje_regenerado = (gasoil_regenerado_litros / consumo_litros) * 100
    
    # 5. Resultados consolidados
    resultados = {
        "km_totales_semana": km_totales,
        "litros_consumidos_semana": consumo_litros,
        "co2_producido_semana": co2_captura['co2_total_producido_kg'],
        "co2_capturado_semana": co2_capturado_semana,
        "co2_almacenado_semana": co2_almacenado_semana,
        "co2_emitido_semana": co2_emitido_final,
        "gasoil_regenerado_litros": gasoil_regenerado_litros,
        "porcentaje_gasoil_regenerado": porcentaje_regenerado,
        "reduccion_co2_total": (1 - (co2_emitido_final / co2_captura['co2_total_producido_kg'])) * 100,
        "energia_kilometro_kwh": sintesis['energia_disponible_kwh'],
        "energia_sintesis_kwh": sintesis['energia_necesaria_kwh'],
    }
    
    return resultados

# --- EJECUTAR SIMULACIÓN ---
resultados_semana = simular_ciclo_completo_semanal()

# --- IMPRIMIR RESULTADOS ---
print("="*60)
print("RESULTADOS DEL CICLO COMPLETO (SIMULACIÓN SEMANAL)")
print("="*60)
print(f"Kilómetros totales:             {resultados_semana['km_totales_semana']:.0f} km")
print(f"Gasoil consumido:               {resultados_semana['litros_consumidos_semana']:.2f} L")
print(f"CO2 producido:                  {resultados_semana['co2_producido_semana']:.2f} kg")
print(f"CO2 capturado:                  {resultados_semana['co2_capturado_semana']:.2f} kg")
print(f"CO2 emitido final:              {resultados_semana['co2_emitido_semana']:.2f} kg")
print(f"Gasoil regenerado:              {resultados_semana['gasoil_regenerado_litros']:.2f} L")
print(f"Porcentaje de gasoil regenerado: {resultados_semana['porcentaje_gasoil_regenerado']:.1f}%")
print(f"REDUCCIÓN DE EMISIONES TOTAL:    {resultados_semana['reduccion_co2_total']:.1f}%")
print(f"Energía Kilómetro (disponible):  {resultados_semana['energia_kilometro_kwh']:.2f} kWh")
print(f"Energía Síntesis (necesaria):   {resultados_semana['energia_sintesis_kwh']:.2f} kWh")
print("="*60)

with open('gemelo_3_ciclo_completo_resultados.json', 'w') as f:
    json.dump(resultados_semana, f, indent=4)