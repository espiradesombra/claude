import math
import random
import time
from dataclasses import dataclass
from typing import List, Tuple

# ============================================================
# 1. ESTRUCTURA DE INFERENCIA (AntiPC para vehículo)
# ============================================================

@dataclass
class MarchaVirtual:
    """Representa una marcha virtual con su perfil de par y eficiencia."""
    nombre: str
    par_min: float      # Par mínimo (Nm)
    par_max: float      # Par máximo (Nm)
    eficiencia: float   # 0.0 a 1.0
    rango_velocidad: Tuple[float, float]  # (min_kmh, max_kmh)

class AntiPC_Vehiculo:
    """
    Motor de inferencia para control de vehículo eléctrico.
    Calcula la combinación óptima de par y marcha virtual.
    """
    def __init__(self):
        # Definición de marchas virtuales (5 marchas)
        self.marchas = [
            MarchaVirtual("M1", 400, 600, 0.92, (0, 30)),
            MarchaVirtual("M2", 300, 450, 0.95, (20, 60)),
            MarchaVirtual("M3", 200, 350, 0.97, (40, 90)),
            MarchaVirtual("M4", 150, 250, 0.98, (60, 120)),
            MarchaVirtual("M5", 100, 180, 0.99, (90, 200)),
        ]
        
        # Buffer de inferencias (sobrantes)
        self.buffer = {
            'ultimo_par': 0,
            'ultima_marcha': self.marchas[0],
            'eficiencia_acumulada': 0,
            'combinaciones_probadas': 0
        }
        
        # Umbrales para inferencia mínima absurda
        self.umbral_cambio_marcha = 0.15   # 15% de cambio para cambiar marcha
        self.umbral_par = 20               # Nm mínimo de cambio

    def inferir_par_optimo(self, velocidad: float, pendiente: float, carga: float, temperatura: float) -> float:
        """
        Infiere el par óptimo usando el principio de "cálculo mínimo absurdo".
        """
        # 1. Calcular el par base (mínimo absurdo)
        par_base = 200 + (velocidad * 0.5) + (pendiente * 50) + (carga * 0.1) - (temperatura * 2)
        par_base = max(50, min(600, par_base))  # Limitar entre 50 y 600 Nm

        # 2. Buscar en el buffer si ya tenemos un par similar (sobrante)
        if abs(par_base - self.buffer['ultimo_par']) < self.umbral_par:
            # Reutilizar el par anterior (ahorro de cálculo)
            par_inferido = self.buffer['ultimo_par']
            self.buffer['combinaciones_probadas'] += 1
        else:
            # 3. Inferir el par mediante mapeo (no cálculo directo)
            # Combinación de factores: velocidad, pendiente, carga, temperatura
            factor_vel = velocidad / 200          # 0 a 1
            factor_pend = pendiente / 30          # 0 a 1 (pendiente máxima 30%)
            factor_carga = carga / 5000           # 0 a 1 (carga máxima 5000 kg)
            factor_temp = (temperatura - 20) / 60 # -0.33 a 0.67 (20°C a 80°C)

            # Combinación no lineal (polinomio de 2º grado)
            combinacion = (0.3 * factor_vel + 0.3 * factor_pend + 
                          0.2 * factor_carga - 0.2 * factor_temp)
            combinacion = max(0, min(1, combinacion))

            # Mapeo: combinación → par (relación no armónica)
            par_inferido = 100 + combinacion * 500  # 100-600 Nm
            par_inferido = max(50, min(600, par_inferido))

            # Guardar en buffer (permutar = guardar)
            self.buffer['ultimo_par'] = par_inferido
            self.buffer['combinaciones_probadas'] += 1

        return par_inferido

    def inferir_marcha_virtual(self, velocidad: float, par_inferido: float) -> MarchaVirtual:
        """
        Infiere la marcha virtual óptima.
        """
        # 1. Buscar marcha por velocidad
        marcha_candidata = self.marchas[0]
        for m in self.marchas:
            if m.rango_velocidad[0] <= velocidad <= m.rango_velocidad[1]:
                marcha_candidata = m
                break

        # 2. Ajustar por par (si el par es alto, bajar marcha)
        if par_inferido > 400 and marcha_candidata.nombre in ["M4", "M5"]:
            # Buscar marcha inferior con mayor par
            for m in self.marchas:
                if m.par_max >= par_inferido:
                    marcha_candidata = m
                    break

        # 3. Verificar si hay un cambio significativo (umbral)
        if self.buffer['ultima_marcha'].nombre != marcha_candidata.nombre:
            # Calcular "distancia" entre marchas (índice)
            idx_actual = self.marchas.index(self.buffer['ultima_marcha'])
            idx_nueva = self.marchas.index(marcha_candidata)
            distancia = abs(idx_actual - idx_nueva) / len(self.marchas)

            if distancia < self.umbral_cambio_marcha:
                # No cambiar (sobrante de la decisión anterior)
                marcha_candidata = self.buffer['ultima_marcha']
            else:
                # Guardar nueva marcha en buffer
                self.buffer['ultima_marcha'] = marcha_candidata

        return marcha_candidata

    def ejecutar_ciclo(self, velocidad: float, pendiente: float, carga: float, temperatura: float) -> dict:
        """
        Ejecuta un ciclo completo de inferencia.
        """
        # 1. Inferir par
        par = self.inferir_par_optimo(velocidad, pendiente, carga, temperatura)

        # 2. Inferir marcha
        marcha = self.inferir_marcha_virtual(velocidad, par)

        # 3. Calcular eficiencia (combinatoria)
        eficiencia_base = marcha.eficiencia
        factor_par = 1.0 - abs(par - marcha.par_min) / (marcha.par_max - marcha.par_min + 1)
        eficiencia = eficiencia_base * (0.8 + 0.2 * factor_par)
        eficiencia = max(0.5, min(1.0, eficiencia))

        # 4. Actualizar buffer de eficiencia
        self.buffer['eficiencia_acumulada'] = (self.buffer['eficiencia_acumulada'] + eficiencia) / 2

        return {
            'velocidad': velocidad,
            'pendiente': pendiente,
            'carga': carga,
            'temperatura': temperatura,
            'par_inferido': round(par, 1),
            'marcha': marcha.nombre,
            'eficiencia': round(eficiencia * 100, 1),
            'eficiencia_media': round(self.buffer['eficiencia_acumulada'] * 100, 1),
            'combinaciones_probadas': self.buffer['combinaciones_probadas']
        }


# ============================================================
# 2. SIMULACIÓN DE CONDUCCIÓN
# ============================================================

def simular_conduccion(duracion_segundos: int = 60, intervalo: float = 0.5):
    """
    Simula un trayecto con cambios de velocidad, pendiente, carga y temperatura.
    """
    anti_pc = AntiPC_Vehiculo()
    historial = []

    velocidad = 0
    pendiente = 0
    carga = 2000  # kg
    temperatura = 25  # °C

    print("=" * 70)
    print("🚗 SIMULACIÓN DE INFERENCIA ANTI-PC EN VEHÍCULO ELÉCTRICO")
    print("=" * 70)
    print(f"{'Tiempo':<8} | {'Vel':<6} | {'Pend':<6} | {'Carga':<7} | {'Temp':<5} | {'Par':<8} | {'Marcha':<6} | {'Efic':<6} | {'Media':<6}")
    print("-" * 70)

    for t in range(0, duracion_segundos, int(intervalo * 10)):
        # Variar condiciones de conducción (simulación)
        velocidad = 30 + 70 * (0.5 + 0.5 * math.sin(t / 10))  # 30 a 100 km/h
        pendiente = 5 * math.sin(t / 15)  # -5% a 5%
        carga = 2000 + 500 * math.sin(t / 20)  # 1500 a 2500 kg
        temperatura = 25 + 10 * math.sin(t / 30)  # 15 a 35 °C

        resultado = anti_pc.ejecutar_ciclo(velocidad, pendiente, carga, temperatura)
        historial.append(resultado)

        print(f"{t/10:>6.1f}s | {velocidad:>5.1f} | {pendiente:>5.1f} | {carga:>7.0f} | {temperatura:>4.1f} | {resultado['par_inferido']:>8.1f} | {resultado['marcha']:>6} | {resultado['eficiencia']:>5.1f}% | {resultado['eficiencia_media']:>5.1f}%")

    # Resumen final
    print("-" * 70)
    eficiencia_media_total = sum(h['eficiencia'] for h in historial) / len(historial)
    print(f"📊 Eficiencia media total: {eficiencia_media_total:.1f}%")
    print(f"🔄 Combinaciones inferidas: {anti_pc.buffer['combinaciones_probadas']}")
    print("=" * 70)

    return historial


# ============================================================
# 3. EJECUCIÓN
# ============================================================

if __name__ == "__main__":
    simular_conduccion(duracion_segundos=120, intervalo=0.5)