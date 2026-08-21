# FILE: analisis_roi.py
# Objetivo: Analizar la viabilidad financiera del proyecto.

import matplotlib.pyplot as plt

# --- COSTES DEL SISTEMA (por vehículo) ---
COSTES_SISTEMA = {
    "DAC System": 3000,
    "Depósito CO2": 2000,
    "Reactor Sabatier": 5000,
    "Electrolizador": 3000,
    "Kilómetro (Batería)": 8000,
    "Sensores y Control": 1500,
    "Integración y Montaje": 2500,
}
COSTE_TOTAL_SISTEMA = sum(COSTES_SISTEMA.values())  # 25,000€

# --- AHORROS E INGRESOS ANUALES ---
# Caso: Coche Privado (15,000 km/año)
CONSUMO_ANUAL_COCHE = 1000  # Litros/año
PRECIO_GASOIL = 1.40  # €/litro

# Con el sistema (20% menos consumo)
ahorro_combustible = (CONSUMO_ANUAL_COCHE * 0.2) * PRECIO_GASOIL
# Por el gasoil sintético premium
ingreso_premium = (CONSUMO_ANUAL_COCHE * 0.2) * 0.30  # 0.30€ extra por litro verde
# Créditos de carbono
ahorro_co2_ton = (CONSUMO_ANUAL_COCHE * 3.15 * 0.80) / 1000  # 0.8 es la eficiencia de captura
ingreso_carbono = ahorro_co2_ton * 100  # 100€/ton

BENEFICIO_ANUAL_COCHE = ahorro_combustible + ingreso_premium + ingreso_carbono

# Caso: Camión (80,000 km/año)
CONSUMO_ANUAL_CAMION = 5000  # Litros/año
COSTE_SISTEMA_CAMION = COSTE_TOTAL_SISTEMA * 1.5  # 50% más caro

ahorro_combustible_camion = (CONSUMO_ANUAL_CAMION * 0.2) * PRECIO_GASOIL
ingreso_premium_camion = (CONSUMO_ANUAL_CAMION * 0.2) * 0.30
ahorro_co2_ton_camion = (CONSUMO_ANUAL_CAMION * 3.15 * 0.80) / 1000
ingreso_carbono_camion = ahorro_co2_ton_camion * 100

BENEFICIO_ANUAL_CAMION = ahorro_combustible_camion + ingreso_premium_camion + ingreso_carbono_camion

# --- CÁLCULO DEL PAYBACK ---
payback_coche = COSTE_TOTAL_SISTEMA / BENEFICIO_ANUAL_COCHE
payback_camion = COSTE_SISTEMA_CAMION / BENEFICIO_ANUAL_CAMION

# --- PROYECCIÓN DE VENTAS DE REPSOL ---
# Unidades vendidas (millones de vehículos/año) - Escenario optimista
años = list(range(2027, 2038))
ventas_coches = [0.01, 0.05, 0.15, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0]  # Millones
ventas_camiones = [0.001, 0.005, 0.015, 0.03, 0.05, 0.08, 0.1, 0.12, 0.15, 0.18, 0.2]  # Millones

# Ingresos y costes para Repsol
INGRESO_POR_COCHE = 5000  # Margen de 5k€ por coche
INGRESO_POR_CAMION = 15000
COSTE_RD = [150, 200, 150, 100, 50, 50, 50, 50, 50, 50, 50]  # Millones €

ingresos_anuales = []
costes_anuales = []
flujo_caja = []

for i in range(len(años)):
    ing_coches = ventas_coches[i] * 1_000_000 * INGRESO_POR_COCHE
    ing_camiones = ventas_camiones[i] * 1_000_000 * INGRESO_POR_CAMION
    ingresos_totales = (ing_coches + ing_camiones) / 1_000_000  # Millones €
    coste_anual = COSTE_RD[i]
    flujo = ingresos_totales - coste_anual
    ingresos_anuales.append(ingresos_totales)
    costes_anuales.append(coste_anual)
    flujo_caja.append(flujo)

# Acumulado
flujo_acumulado = []
acum = 0
for f in flujo_caja:
    acum += f
    flujo_acumulado.append(acum)

# --- IMPRIMIR RESULTADOS ---
print("="*60)
print("ANÁLISIS FINANCIERO DEL PROYECTO GASOIL VERDE")
print("="*60)

print("\n--- COSTES DEL SISTEMA (por vehículo) ---")
for componente, coste in COSTES_SISTEMA.items():
    print(f"{componente:20s}: {coste:>5,} €")
print(f"\nCOSTE TOTAL SISTEMA: {COSTE_TOTAL_SISTEMA:>18,} €")
print(f"COSTE TOTAL CAMIÓN:  {COSTE_SISTEMA_CAMION:>18,} €")

print("\n--- BENEFICIOS ANUALES (por vehículo) ---")
print(f"Beneficio anual (Coche): €{BENEFICIO_ANUAL_COCHE:.0f}")
print(f"Payback (Coche): {payback_coche:.1f} años")
print(f"Beneficio anual (Camión): €{BENEFICIO_ANUAL_CAMION:.0f}")
print(f"Payback (Camión): {payback_camion:.1f} años")

print("\n--- PROYECCIÓN REPSOL (Millones €) ---")
for i, año in enumerate(años):
    print(f"{año}: Ingresos {ingresos_anuales[i]:>5.0f} | Costes {costes_anuales[i]:>5.0f} | Flujo {flujo_caja[i]:>5.0f} | Acumulado {flujo_acumulado[i]:>6.0f}")

print(f"\nVPN estimado (10 años): ~€{sum(flujo_acumulado)/len(flujo_acumulado):.0f}M")

# --- GRÁFICA ---
plt.figure(figsize=(14, 8))

plt.subplot(2, 2, 1)
plt.bar(COSTES_SISTEMA.keys(), COSTES_SISTEMA.values(), color='orange')
plt.title('Coste del Sistema por Componente')
plt.ylabel('Coste (€)')
plt.xticks(rotation=45, ha='right')

plt.subplot(2, 2, 2)
plt.bar(['Coche', 'Camión'], [BENEFICIO_ANUAL_COCHE, BENEFICIO_ANUAL_CAMION], color=['blue', 'green'])
plt.title('Beneficio Anual por Tipo de Vehículo')
plt.ylabel('Beneficio (€/año)')

plt.subplot(2, 2, 3)
plt.plot(años, flujo_caja, marker='o', label='Flujo de Caja Anual')
plt.plot(años, flujo_acumulado, marker='s', label='Flujo Acumulado')
plt.axhline(y=0, color='r', linestyle='--')
plt.title('Flujo de Caja del Proyecto (Repsol)')
plt.xlabel('Año')
plt.ylabel('Millones €')
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 4)
plt.pie([COSTE_TOTAL_SISTEMA, BENEFICIO_ANUAL_COCHE*5, BENEFICIO_ANUAL_CAMION*2],
        labels=['Coste Sistema', 'Beneficio 5 Años (Coche)', 'Beneficio 2 Años (Camión)'],
        autopct='%1.1f%%', startangle=90)
plt.title('Distribución de Costes vs Beneficios')

plt.tight_layout()
plt.savefig('analisis_financiero.png')
print("\nGráfico guardado como 'analisis_financiero.png'")