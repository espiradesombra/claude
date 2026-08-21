# FILE: gemelo_2_sintesis.py
# Objetivo: Simular la síntesis de gasoil a partir de CO2 capturado.
# Modelo: Electrólisis + Reactor Sabatier + Fischer-Tropsch.
# Entrada: CO2 disponible y energía del Kilómetro.
# Salida: Litros de gasoil sintético regenerado.

import json

# --- CONSTANTES QUÍMICAS Y ENERGÉTICAS ---
# Masa molar: CO2 = 44 kg/kmol, CH4 = 16 kg/kmol, H2 = 2 kg/kmol
# 1 kmol CO2 + 4 kmol H2 -> 1 kmol CH4 + 2 kmol H2O
# 4 CH4 -> 1 C12H23 (gasoil) (aproximación)

CO2_KG_POR_KMOL = 44
CH4_KG_POR_KMOL = 16
H2_KG_POR_KMOL = 2
GASOIL_KG_POR_KMOL = 170  # Masa molar aproximada (C12H23)

# Energía requerida para la síntesis
# 46 MJ/mol CO2 = 12.8 kWh/mol CO2 (teórico)
ENERGIA_KWH_POR_KMOL_CO2 = 12.8 * 1000  # 12,800 kWh/kmol CO2
ENERGIA_KWH_POR_KG_CO2 = ENERGIA_KWH_POR_KMOL_CO2 / CO2_KG_POR_KMOL
# = 12,800 / 44 ≈ 290.9 kWh/kg CO2

# --- PARÁMETROS DEL KILÓMETRO ---
ENERGIA_DISPONIBLE_KWH_POR_HORA = 1.5  # Genera 1-2 kWh por hora

# --- FUNCIÓN DE SÍNTESIS ---
def sintetizar_gasoil(co2_disponible_kg, horas_conduccion):
    """Calcula cuánto gasoil se puede sintetizar."""
    
    # 1. Energía disponible del Kilómetro
    energia_total_kwh = ENERGIA_DISPONIBLE_KWH_POR_HORA * horas_conduccion
    
    # 2. CO2 que se puede sintetizar con la energía disponible
    co2_sintetizable_kg = energia_total_kwh / ENERGIA_KWH_POR_KG_CO2
    co2_a_sintetizar_kg = min(co2_disponible_kg, co2_sintetizable_kg)
    
    # 3. Química de la síntesis
    # Convertir kg CO2 a kmol
    co2_kmol = co2_a_sintetizar_kg / CO2_KG_POR_KMOL
    
    # 3a. Reactor Sabatier: CO2 + 4H2 -> CH4 + 2H2O
    ch4_kmol = co2_kmol  # Relación 1:1
    h2_necesario_kmol = co2_kmol * 4
    # 3b. Fischer-Tropsch: 4CH4 -> C12H23
    gasoil_kmol = ch4_kmol / 4  # Relación 4:1
    gasoil_kg = gasoil_kmol * GASOIL_KG_POR_KMOL
    gasoil_litros = gasoil_kg / 0.85  # Densidad del gasoil ≈ 0.85 kg/L
    
    # 4. Energía realmente necesaria vs disponible
    energia_necesaria_kwh = co2_kmol * ENERGIA_KWH_POR_KMOL_CO2
    
    # 5. Balance
    balance_energetico_kwh = energia_total_kwh - energia_necesaria_kwh
    es_autosuficiente = balance_energetico_kwh >= 0
    
    # 6. Resultados
    resultados = {
        "co2_disponible_kg": co2_disponible_kg,
        "co2_sintetizado_kg": co2_a_sintetizar_kg,
        "energia_disponible_kwh": energia_total_kwh,
        "energia_necesaria_kwh": energia_necesaria_kwh,
        "balance_energetico_kwh": balance_energetico_kwh,
        "es_autosuficiente": es_autosuficiente,
        "gasoil_sintetizado_kg": gasoil_kg,
        "gasoil_sintetizado_litros": gasoil_litros,
        "porcentaje_gasoil_regenerado": (gasoil_litros / 30) * 100,  # 30L depósito
        "ch4_producido_kg": ch4_kmol * CH4_KG_POR_KMOL,
        "h2_necesario_kg": h2_necesario_kmol * H2_KG_POR_KMOL,
    }
    
    return resultados

# --- EJECUTAR SIMULACIÓN ---
co2_capturado_ejemplo = 12.0  # kg (resultado típico de gemelo_1)
horas_viaje = 8

resultados = sintetizar_gasoil(co2_capturado_ejemplo, horas_viaje)

# --- IMPRIMIR Y GUARDAR RESULTADOS ---
print("="*60)
print("RESULTADOS DE LA SIMULACIÓN DE SÍNTESIS DE GASOIL")
print("="*60)
print(f"CO2 disponible para síntesis:   {resultados['co2_disponible_kg']:.2f} kg")
print(f"CO2 realmente sintetizado:      {resultados['co2_sintetizado_kg']:.2f} kg")
print(f"Energía disponible (Kilómetro): {resultados['energia_disponible_kwh']:.2f} kWh")
print(f"Energía necesaria (Síntesis):   {resultados['energia_necesaria_kwh']:.2f} kWh")
print(f"Balance energético:             {resultados['balance_energetico_kwh']:.2f} kWh")
print(f"¿Es autosuficiente?             {'SÍ' if resultados['es_autosuficiente'] else 'NO'}")
print(f"Gasoil sintetizado:             {resultados['gasoil_sintetizado_litros']:.3f} L")
print(f"Porcentaje de depósito (30L):   {resultados['porcentaje_gasoil_regenerado']:.1f}%")
print("="*60)

# Guardar resultados en JSON
with open('gemelo_2_sintesis_resultados.json', 'w') as f:
    json.dump(resultados, f, indent=4)