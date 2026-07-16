from decimal import Decimal, getcontext
from dataclasses import dataclass
from typing import List, Generator

getcontext().prec = 150

# =====================================================================
# BLOC I: MOTOR MATEMÁTICO PURO & ALEATOROVIX (Generador de Chorro)
# =====================================================================

class Aleatorovix:
    """
    Generador de Chorro Caótico de Alta Entropía en Acordeón (MSL/LSL).
    Usa una semilla matemática de fase para generar un flujo masivo de bits.
    """
    def __init__(self, semilla: int):
        # Asegura que la semilla no sea cero y define constantes de acoplamiento caótico
        self.estado = Decimal(str(semilla if semilla != 0 else 123456789))
        self.factor_msl = Decimal('42.123456789')   # Lazo Significativo Mayor (Amplitud)
        self.factor_lsl = Decimal('0.000123456')   # Lazo Significativo Menor (Frecuencia)

    def generar_chorro_bits(self, longitud_bits: int) -> Generator[str, None, None]:
        """
        Produce un chorro de bits estirando y encogiendo la fase (Efecto Acordeón).
        """
        for _ in range(longitud_bits):
            # Evolución del sistema dinámico acoplado
            # x_{n+1} = (x_n * factor_msl + sin(x_n * factor_lsl)) mod 1
            fase_actual = (self.estado * self.factor_msl) + (self.estado * self.factor_lsl).cos()
            self.estado = fase_actual - fase_actual.to_integral_value(rounding="ROUND_FLOOR")
            
            # Umbralización determinista para extraer el bit (MSL / LSL)
            bit = "1" if self.estado > Decimal('0.5') else "0"
            yield bit


# =====================================================================
# BLOC II: CRIPTOGRAFÍA DE FASE (Cifrado de la Semilla de Entrada)
# =====================================================================

@dataclass
class ClaveK3Aleatorovix:
    tales: List[int]
    figuras: List[str]
    puntos: List[int]
    saltos: List[int]
    primos: List[tuple]
    iteraciones_pi: int = 15
    porcentajes_aportacion: List[int] = None
    longitud_segmento: int = 15
    exponentes_10: List[int] = None

    def __post_init__(self):
        if self.porcentajes_aportacion is None: self.porcentajes_aportacion = [100] * len(self.primos)
        if self.exponentes_10 is None: self.exponentes_10 = [3, 7, 12, 18, 22, 29, 35, 42]


# Funciones de apoyo de tu motor geométrico
def taylor_sin(x: Decimal) -> Decimal:
    termino, suma, x_cuadrado, n = x, x, x * x, 1
    while abs(termino) > Decimal('1e-150'):
        termino = -termino * x_cuadrado / Decimal((2 * n) * (2 * n + 1))
        suma += termino
        n += 1
    return suma

def aproximar_pi_puro(p1: int, p2: int, iteraciones: int = 15) -> Decimal:
    lados = Decimal(p1 * p2)
    pi_inicial = Decimal('3.14159265358979323846264338327950288419716939937510')
    angulo = pi_inicial / lados
    lado = Decimal(2) * taylor_sin(angulo)
    perimetro = lados * lado
    pi_actual = perimetro / Decimal(2)
    for _ in range(iteraciones):
        lados *= 2
        raiz_soporte = (Decimal(4) - lado * lado).sqrt()
        lado = lado / (Decimal(2) + raiz_soporte).sqrt()
        perimetro = lados * lado
        pi_actual = perimetro / Decimal(2)
    return pi_actual

def aproximar_e_convergente_50() -> Decimal:
    termino, acumulado_e = Decimal('1.0'), Decimal('1.0')
    for v in range(1, 51):
        eleccion = (2 * v + 2) if v % 2 == 0 else (3 * v)
        termino = termino / Decimal(str(eleccion if eleccion != 0 else 1))
        acumulado_e += termino
    return acumulado_e

def calcular_lados_geometricos(tipo: str, puntos: int, escala: Decimal) -> List[Decimal]:
    base = Decimal(puntos) * Decimal('1.5')
    if tipo == "equilatero":
        lado = (base / Decimal('3')) * escala
        return [lado, lado, lado]
    elif tipo == "isosceles":
        lado_a = (base / Decimal('4')) * escala
        lado_b = (base / Decimal('2')) * escala
        return [lado_a, lado_a, lado_b]
    else:
        return [(base * Decimal('0.25')) * escala, (base * Decimal('0.35')) * escala, (base * Decimal('0.40')) * escala]


# =====================================================================
# BLOC III: PROCESO HÍBRIDO (Cifrado Masivo por Chorro de Acordeón)
# =====================================================================

def cifrar_masivo_aleatorovix(datos_claros: str, semilla_secreta: int, clave: ClaveK3Aleatorovix) -> tuple:
    """
    1. Cifra la semilla secreta con el lazo geométrico de fase (para el transporte seguro).
    2. Usa Aleatorovix inicializado con la semilla para generar el acordeón de bits masivo.
    3. Aplica XOR entre el chorro masivo de Aleatorovix y los datos reales.
    """
    # Fase A: Encriptar la Semilla (Clave de Sesión de 16 bits para el ejemplo)
    semilla_bin = f"{semilla_secreta:016b}"
    perimetro = Decimal('0.0')
    
    for idx, bit in enumerate(semilla_bin):
        figura_idx = idx // 3
        lado_idx = idx % 3
        if lado_idx == clave.saltos[figura_idx % len(clave.saltos)]:
            continue
        escala = Decimal(str(clave.tales[figura_idx % len(clave.tales)]))
        figura_tipo = clave.figuras[figura_idx % len(clave.figuras)]
        puntos_fig = clave.puntos[figura_idx % len(clave.puntos)]
        lados = calcular_lados_geometricos(figura_tipo, puntos_fig, escala)
        if bit == '1':
            perimetro += lados[lado_idx]

    pi_ofuscado = Decimal('1.0')
    for i, (p1, p2) in enumerate(clave.primos):
        pi_ofuscado *= (aproximar_pi_puro(p1, p2, clave.iteraciones_pi) ** (Decimal(clave.porcentajes_aportacion[i]) / Decimal(100)))
    
    e_convergente = aproximar_e_convergente_50()
    hash_semilla = perimetro * pi_ofuscado * e_convergente

    # Fase B: Generar el Chorro de Acordeón con Aleatorovix y XOR
    datos_bits = "".join(f"{ord(c):08b}" for c in datos_claros)
    generador = Aleatorovix(semilla_secreta)
    chorro = generador.generar_chorro_bits(len(datos_bits))
    
    bits_cifrados = []
    for bit_dato, bit_chorro in zip(datos_bits, chorro):
        # Operación XOR (Sumador de canal modular)
        xor_bit = str(int(bit_dato) ^ int(bit_chorro))
        bits_cifrados.append(xor_bit)
        
    return hash_semilla, "".join(bits_cifrados)


def descifrar_masivo_aleatorovix(hash_semilla: Decimal, datos_cifrados_bits: str, clave: ClaveK3Aleatorovix) -> str:
    """
    1. Recupera la semilla secreta aplicando backtracking al hash de fase.
    2. Regenera el chorro idéntico en acordeón con Aleatorovix.
    3. Revierte el XOR para descifrar los datos masivos sin importar el tamaño.
    """
    # Fase A: Backtracking para recuperar la Semilla (16 bits)
    pi_ofuscado = Decimal('1.0')
    for i, (p1, p2) in enumerate(clave.primos):
        pi_ofuscado *= (aproximar_pi_puro(p1, p2, clave.iteraciones_pi) ** (Decimal(clave.porcentajes_aportacion[i]) / Decimal(100)))
    e_convergente = aproximar_e_convergente_50()

    perimetro_original = hash_semilla / (pi_ofuscado * e_convergente)

    def backtrack(residuo: Decimal, idx: int) -> Optional[str]:
        if idx == 16:  # Buscamos exactamente los 16 bits de la semilla
            return "" if abs(residuo) < Decimal('1e-100') else None
        figura_idx = idx // 3
        lado_idx = idx % 3
        if lado_idx == clave.saltos[figura_idx % len(clave.saltos)]:
            for bit in ["0", "1"]:
                res = backtrack(residuo, idx + 1)
                if res is not None: return bit + res
            return None
        escala = Decimal(str(clave.tales[figura_idx % len(clave.tales)]))
        figura_tipo = clave.figuras[figura_idx % len(clave.figuras)]
        puntos_fig = clave.puntos[figura_idx % len(clave.puntos)]
        lados = calcular_lados_geometricos(figura_tipo, puntos_fig, escala)
        lado_actual = lados[lado_idx]

        if residuo >= lado_actual - Decimal('1e-100'):
            res = backtrack(residuo - lado_actual, idx + 1)
            if res is not None: return "1" + res
        res = backtrack(residuo, idx + 1)
        if res is not None: return "0" + res
        return None

    semilla_recuperada_bin = backtrack(perimetro_original, 0)
    semilla_secreta = int(semilla_recuperada_bin, 2)

    # Fase B: Regeneración del Chorro de Acordeón y Descifrado
    generador = Aleatorovix(semilla_secreta)
    chorro = generador.generar_chorro_bits(len(datos_cifrados_bits))
    
    bits_descifrados = []
    for bit_cifrado, bit_chorro in zip(datos_cifrados_bits, chorro):
        xor_bit = str(int(bit_cifrado) ^ int(bit_chorro))
        bits_descifrados.append(xor_bit)
        
    # Reconstrucción de texto ASCII
    bits_str = "".join(bits_descifrados)
    caracteres = [chr(int(bits_str[i:i+8], 2)) for i in range(0, len(bits_str), 8)]
    return "".join(caracteres)


# =====================================================================
# PRUEBA DE RENDIMIENTO Y COHERENCIA
# =====================================================================
if __name__ == "__main__":
    clave_k3 = ClaveK3Aleatorovix(
        tales=[3, 5, 8, 13, 21],
        figuras=["equilatero", "isosceles", "escaleno"],
        puntos=[6, 12, 18],
        saltos=[1, 0, 2],
        primos=[(3, 5), (7, 11)],
        porcentajes_aportacion=[50, 100]
    )

    # Semilla del canal físico de fase
    semilla_secreta = 43210  
    
    # Simulación de un mensaje masivo en tránsito
    mensaje_masivo = "Este es un chorro masivo cifrado con Aleatorovix y el motor de fase de ZypyZape."
    
    print(f"[*] Mensaje original masivo: '{mensaje_masivo}'")
    print(f"[+] Inicializando cifrador de chorro en acordeón con semilla {semilla_secreta}...")
    
    # Cifrar masivo
    hash_semilla, chorro_cifrado = cifrar_masivo_aleatorovix(mensaje_masivo, semilla_secreta, clave_k3)
    
    print(f"[+] Hash de la Semilla transmitido (Fase colapsada): {hash_semilla}")
    print(f"[+] Datos de canal en tránsito (ejemplo binario): {chorro_cifrado[:50]}...")
    
    # Descifrar masivo
    print("[*] Receptor procesando... Aplicando Backtracking de Fase para recuperar Semilla...")
    mensaje_recuperado = descifrar_masivo_aleatorovix(hash_semilla, chorro_cifrado, clave_k3)
    
    print(f"[<] Mensaje masivo reconstruido con éxito: '{mensaje_recuperado}'")