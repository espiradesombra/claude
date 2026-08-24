import math
import random
import time
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from collections import deque
import threading
import queue

# ============================================================
# 1. NODO ANTI-PC (Cada robotaxi es un nodo)
# ============================================================

@dataclass
class NodoRobotaxi:
    """
    Cada robotaxi es un nodo de inferencia en la red AntiPC.
    """
    id: str
    ubicacion: Tuple[float, float]  # (lat, lon)
    estado: str = "disponible"      # disponible, en_ruta, recargando
    bateria: float = 100.0          # 0-100%
    velocidad: float = 0.0          # km/h
    destino: Optional[Tuple[float, float]] = None
    ruta: List[Tuple[float, float]] = field(default_factory=list)
    
    # Buffer local (sobrantes)
    buffer_inferencias: deque = field(default_factory=lambda: deque(maxlen=50))
    
    # Mapeos locales (radios y órbitas)
    mapeo_estado: Dict[str, float] = field(default_factory=dict)
    orbitas_activas: List[Dict] = field(default_factory=list)

    def inferir_estado(self, entrada: Dict) -> Dict:
        """
        Inferencia local del robotaxi.
        Aplica el principio de "cálculo mínimo absurdo".
        """
        # 1. Calcular estado base (mínimo absurdo)
        bateria_abs = min(100, self.bateria + entrada.get('carga_prevista', 0))
        demanda = entrada.get('demanda', 0.5)  # 0-1
        
        # 2. Mapeo de fase (no armónico)
        fase_bateria = (bateria_abs / 100) ** 1.5
        fase_demanda = demanda ** 0.8
        
        # 3. Combinación (órbita)
        estado_inferido = {
            'bateria_efectiva': 50 + 45 * fase_bateria - 10 * fase_demanda,
            'disponibilidad': max(0, min(1, fase_bateria - 0.3 * fase_demanda)),
            'prioridad': 0.5 + 0.4 * fase_demanda - 0.2 * (1 - fase_bateria),
            'tiempo_estimado': 10 + 20 * (1 - fase_bateria) + 5 * (1 - fase_demanda)
        }
        
        # 4. Guardar en buffer (sobrante)
        self.buffer_inferencias.append(estado_inferido)
        
        # 5. Actualizar mapeo local
        self.mapeo_estado = estado_inferido
        
        return estado_inferido


# ============================================================
# 2. RED DE INFERENCIA (La nube AntiPC)
# ============================================================

class RedAntiPC_Nube:
    """
    Red distribuida de inferencia AntiPC.
    Coordina múltiples nodos (robotaxis) sin servidor central.
    """
    def __init__(self, num_nodos: int = 10):
        self.nodos: Dict[str, NodoRobotaxi] = {}
        self.buffer_global = deque(maxlen=200)
        self.mapeo_global = {}
        self.orbitas_globales = []
        self.cola_inferencias = queue.Queue()
        self.resultados = {}
        
        # Inicializar nodos
        for i in range(num_nodos):
            nodo_id = f"taxi_{i:03d}"
            lat = 40.0 + random.uniform(-0.5, 0.5)
            lon = -3.7 + random.uniform(-0.5, 0.5)
            self.nodos[nodo_id] = NodoRobotaxi(
                id=nodo_id,
                ubicacion=(lat, lon),
                bateria=random.uniform(60, 100)
            )
    
    def inferir_flota(self, demanda_global: float = 0.5) -> Dict:
        """
        Inferencia global: coordina toda la flota.
        """
        # 1. Recoger inferencias de todos los nodos
        estados_nodos = {}
        for nodo_id, nodo in self.nodos.items():
            entrada = {
                'demanda': demanda_global,
                'carga_prevista': random.uniform(-5, 15)
            }
            estado = nodo.inferir_estado(entrada)
            estados_nodos[nodo_id] = estado
        
        # 2. Mapeo global (fase de toda la flota)
        baterias = [e['bateria_efectiva'] for e in estados_nodos.values()]
        disponibilidades = [e['disponibilidad'] for e in estados_nodos.values()]
        
        fase_global_bateria = sum(baterias) / len(baterias) / 100
        fase_global_disp = sum(disponibilidades) / len(disponibilidades)
        
        # 3. Combinación global (órbita de la flota)
        eficiencia_flota = (0.6 * fase_global_bateria + 0.4 * fase_global_disp) * 100
        disponibilidad_flota = sum(1 for d in disponibilidades if d > 0.5) / len(disponibilidades)
        
        # 4. Identificar nodos críticos (baja batería, baja disponibilidad)
        nodos_criticos = [
            nid for nid, e in estados_nodos.items()
            if e['bateria_efectiva'] < 30 or e['disponibilidad'] < 0.3
        ]
        
        # 5. Generar recomendaciones (salidas múltiples)
        resultado = {
            'eficiencia_flota': round(eficiencia_flota, 1),
            'disponibilidad_flota': round(disponibilidad_flota * 100, 1),
            'nodos_activos': len([n for n in self.nodos.values() if n.estado == "disponible"]),
            'nodos_criticos': nodos_criticos,
            'recomendacion': self._generar_recomendacion(estados_nodos),
            'timestamp': time.time()
        }
        
        # 6. Guardar en buffer global
        self.buffer_global.append(resultado)
        
        return resultado
    
    def _generar_recomendacion(self, estados: Dict) -> str:
        """
        Genera una recomendación basada en la inferencia global.
        """
        baterias = [e['bateria_efectiva'] for e in estados.values()]
        disp = [e['disponibilidad'] for e in estados.values()]
        
        if min(baterias) < 20:
            return "⚠️ Recarga urgente: nodos con batería crítica"
        elif sum(disp) / len(disp) < 0.4:
            return "🔋 Aumentar disponibilidad: redistribuir carga"
        elif max(baterias) > 90 and min(baterias) > 50:
            return "✅ Flota estable: operación normal"
        else:
            return "🔄 Balancear carga: optimizar asignación de rutas"
    
    def asignar_ruta(self, nodo_id: str, destino: Tuple[float, float]) -> bool:
        """
        Asigna una ruta a un robotaxi usando inferencia.
        """
        if nodo_id not in self.nodos:
            return False
        
        nodo = self.nodos[nodo_id]
        if nodo.estado != "disponible":
            return False
        
        # Inferir viabilidad
        entrada = {'demanda': 0.8, 'carga_prevista': -5}
        estado = nodo.inferir_estado(entrada)
        
        if estado['bateria_efectiva'] < 20 or estado['disponibilidad'] < 0.3:
            return False
        
        # Asignar ruta
        nodo.destino = destino
        nodo.estado = "en_ruta"
        nodo.ruta = self._calcular_ruta(nodo.ubicacion, destino)
        
        return True
    
    def _calcular_ruta(self, origen: Tuple[float, float], destino: Tuple[float, float]) -> List[Tuple[float, float]]:
        """
        Calcula una ruta simple (simulación).
        """
        # Ruta directa con puntos intermedios (simulación)
        lat0, lon0 = origen
        lat1, lon1 = destino
        puntos = []
        for i in range(1, 6):
            t = i / 6
            lat = lat0 + (lat1 - lat0) * t + random.uniform(-0.001, 0.001)
            lon = lon0 + (lon1 - lon0) * t + random.uniform(-0.001, 0.001)
            puntos.append((lat, lon))
        return puntos


# ============================================================
# 3. SIMULADOR DE FLOTA DE ROBOTAXIS
# ============================================================

class SimuladorRobotaxis:
    """
    Simula el comportamiento de la flota de robotaxis con AntiPC.
    """
    def __init__(self, num_taxis: int = 10):
        self.red = RedAntiPC_Nube(num_taxis)
        self.historial = []
        self.demanda_global = 0.5
        self.tiempo = 0
    
    def paso_simulacion(self):
        """
        Avanza un paso de la simulación.
        """
        self.tiempo += 1
        self.demanda_global = 0.3 + 0.5 * (0.5 + 0.5 * math.sin(self.tiempo / 20))
        
        # 1. Inferir estado global
        resultado = self.red.inferir_flota(self.demanda_global)
        
        # 2. Asignar rutas aleatorias a nodos disponibles
        nodos_disponibles = [nid for nid, n in self.red.nodos.items() if n.estado == "disponible"]
        if nodos_disponibles and random.random() < 0.3:
            nodo_id = random.choice(nodos_disponibles)
            destino = (40.0 + random.uniform(-0.5, 0.5), -3.7 + random.uniform(-0.5, 0.5))
            self.red.asignar_ruta(nodo_id, destino)
        
        # 3. Actualizar estados de los nodos en ruta
        for nodo in self.red.nodos.values():
            if nodo.estado == "en_ruta":
                nodo.bateria -= random.uniform(0.1, 0.5)
                if nodo.bateria < 10 or random.random() < 0.05:
                    nodo.estado = "disponible"
                    nodo.destino = None
                    nodo.ruta = []
        
        # 4. Recargar nodos con batería baja
        for nodo in self.red.nodos.values():
            if nodo.bateria < 20 and nodo.estado == "disponible":
                nodo.bateria += random.uniform(5, 15)
                if nodo.bateria > 80:
                    nodo.bateria = 80 + random.uniform(0, 10)
        
        # 5. Guardar en historial
        self.historial.append({
            'tiempo': self.tiempo,
            'demanda': self.demanda_global,
            'resultado': resultado
        })
        
        return resultado
    
    def ejecutar_simulacion(self, pasos: int = 50):
        """
        Ejecuta la simulación completa.
        """
        print("=" * 80)
        print("🚕 SIMULACIÓN DE FLOTA DE ROBOTAXIS CON ANTI-PC")
        print("=" * 80)
        print(f"{'Paso':<6} | {'Demanda':<8} | {'Eficiencia':<10} | {'Disponibilidad':<14} | {'Activos':<8} | {'Críticos':<8} | {'Recomendación'}")
        print("-" * 80)
        
        for i in range(pasos):
            resultado = self.paso_simulacion()
            
            print(f"{i+1:>4} | {self.demanda_global:>7.2f} | {resultado['eficiencia_flota']:>9.1f}% | {resultado['disponibilidad_flota']:>13.1f}% | {resultado['nodos_activos']:>6} | {len(resultado['nodos_criticos']):>7} | {resultado['recomendacion'][:30]}...")
            
            time.sleep(0.1)  # Para ver la evolución
        
        print("-" * 80)
        print("📊 SIMULACIÓN COMPLETADA")
        print(f"   - Total de inferencias: {len(self.historial)}")
        print(f"   - Nodos en la flota: {len(self.red.nodos)}")
        print(f"   - Buffer global: {len(self.red.buffer_global)} entradas")
        print("=" * 80)
        
        return self.historial


# ============================================================
# 4. EJECUCIÓN PRINCIPAL
# ============================================================

if __name__ == "__main__":
    # Crear simulador con 15 robotaxis
    simulador = SimuladorRobotaxis(num_taxis=15)
    
    # Ejecutar simulación de 60 pasos (≈ 60 minutos de operación)
    historial = simulador.ejecutar_simulacion(pasos=60)
    
    # Mostrar resumen final de la flota
    print("\n📋 RESUMEN FINAL DE LA FLOTA:")
    print("-" * 80)
    for nodo_id, nodo in simulador.red.nodos.items():
        estado_str = f"{nodo.estado} {'→' + str(nodo.destino) if nodo.destino else ''}"
        print(f"  {nodo_id}: Batería={nodo.bateria:.1f}% | Estado={estado_str}")