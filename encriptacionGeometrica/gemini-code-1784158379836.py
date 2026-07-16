import sys
import os
from decimal import Decimal, getcontext

# Forzamos la precisión del canal de fase para el backtracking
getcontext().prec = 150

# Añadimos tu carpeta nucleo al path para poder importar tus archivos reales
sys.path.append(os.path.abspath("./nucleo"))

try:
    # Simulamos la integración con tus archivos reales del repositorio
    import aleatorovix as ax
    from mascara_lila import aplicar_mascara_lila  # Tu lógica de enmascaramiento
    USAR_MODULOS_REALES = True
except ImportError:
    # Si lo corres fuera de tu estructura local, usamos un fallback compatible
    USAR_MODULOS_REALES = False

# =====================================================================
# BLOC I: BACKTRACKING GEOMÉTRICO (La llave de paso para tu semilla)
# =====================================================================

def recuperar_semilla_fase(hash_semilla: Decimal, clave_k3, dimension_bits: int = 16) -> int:
    """
    Usa el motor matemático de fase para recuperar la semilla que activará
    tu algoritmo de Aleatorovix en tu script local.
    """
    from __main__ import aproximar_pi_puro, aproximar_e_convergente_50, calcular_lados_geometricos
    
    # Desofuscar constante pi de primos y constante e
    pi_ofuscado = Decimal('1.0')
    for i, (p1, p2) in enumerate(clave_k3.primos):
        pi_ofuscado *= (aproximar_pi_puro(p1, p2, clave_k3.iteraciones_pi) ** (Decimal(clave_k3.porcentajes_aportacion[i]) / Decimal(100)))
    e_convergente = aproximar_e_convergente_50()

    perimetro_original = hash_semilla / (pi_ofuscado * e_convergente)

    def backtrack(residuo: Decimal, idx: int) -> str | None:
        if idx == dimension_bits:
            return "" if abs(residuo) < Decimal('1e-100') else None
        
        figura_idx = idx // 3
        lado_idx = idx % 3
        if lado_idx == clave_k3.saltos[figura_idx % len(clave_k3.saltos)]:
            for bit in ["0", "1"]:
                res = backtrack(residuo, idx + 1)
                if res is not None: return bit + res
            return None

        escala = Decimal(str(clave_k3.tales[figura_idx % len(clave_k3.tales)]))
        figura_tipo = clave_k3.figuras[figura_idx % len(clave_k3.figuras)]
        puntos_fig = clave_k3.puntos[figura_idx % len(clave_k3.puntos)]
        lados = calcular_lados_geometricos(figura_tipo, puntos_fig, escala)
        lado_actual = lados[lado_idx]

        if residuo >= lado_actual - Decimal('1e-100'):
            res = backtrack(residuo - lado_actual, idx + 1)
            if res is not None: return "1" + res
            
        res = backtrack(residuo, idx + 1)
        if res is not None: return "0" + res
        return None

    bits_recuperados = backtrack(perimetro_original, 0)
    if bits_recuperados is None:
        raise ValueError("[❌] Error catastrófico: Sincronización de fase perdida.")
    
    return int(bits_recuperados, 2)


# =====================================================================
# BLOC II: LA UNIÓN (Fase + Tu Aleatorovix de nucleo/)
# =====================================================================

def ejecutar_encriptacion_total(mensaje_masivo: str, semilla_secreta: int, clave_k3) -> tuple:
    """
    1. Encripta la semilla en el espacio de fase de los polígonos del Libro 4.
    2. Ejecuta tu módulo Aleatorovix para generar el chorro masivo de bits (MSL/LSL).
    """
    # 1. Encriptación geométrica de la Semilla
    from __main__ import encriptar_token_k3  # Reutiliza tu motor del Bloque anterior
    hash_fase_semilla = encriptar_token_k3(f"{semilla_secreta:016b}", clave_k3)
    
    # 2. Generación del chorro usando tus archivos reales
    if USAR_MODULOS_REALES:
        print("[!] Utilizando tu módulo nativo 'aleatorovix.py' y 'mascara_lila'...")
        # Llama a tu función de máscara sobre el mensaje usando tu generador
        mensaje_cifrado = aplicar_mascara_lila(mensaje_masivo, semilla_secreta)
    else:
        # Fallback local simulado en acordeón MSL/LSL si se corre aislado
        print("[i] Corriendo simulación local de Aleatorovix (Acordeón MSL/LSL)...")
        from __main__ import Aleatorovix
        gen = Aleatorovix(semilla_secreta)
        bits_mensaje = "".join(f"{ord(c):08b}" for c in mensaje_masivo)
        chorro = gen.generar_chorro_bits(len(bits_mensaje))
        mensaje_cifrado = "".join(str(int(b) ^ int(c)) for b, c in zip(bits_mensaje, chorro))
        
    return hash_fase_semilla, mensaje_cifrado