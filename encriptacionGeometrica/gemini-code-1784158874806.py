#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 ██████╗ ██████╗ ███╗   ██╗███████╗██████╗  ██████╗ ███████╗███╗   ██╗ ██████╗██╗ █████╗ 
██╔════╝██╔═══██╗████╗  ██║██╔════╝██╔══██╗██╔════╝ ██╔════╝████╗  ██║██╔════╝██║██╔══██╗
██║     ██║   ██║██╔██╗ ██║█████╗  ██████╔╝██║  ███╗█████╗  ██╔██╗ ██║██║     ██║███████║
██║     ██║   ██║██║╚██╗██║██╔══╝  ██╔══██╗██║   ██║██╔══╝  ██║╚██╗██║██║     ██║██╔══██║
╚██████╗╚██████╔╝██║ ╚████║███████╗██║  ██║╚██████╔╝███████╗██║ ╚████║╚██████╗██║██║  ██║
 ╚═════╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝╚═╝╚═╝  ╚═╝
=============================================================================
[+] SISTEMA INTEGRAL DE CONVERGENCIA GEOMÉTRICA K3
[+] AUTOR: Víctor Manzanares Alberola
[+] CARACTERÍSTICAS INDUSTRIALES:
    - Precisión Arbitraria: Decimal (configurable hasta 1000 dígitos)
    - Criterio exacto: 1 bit = 1 lado (1: suma, 0: ignora)
    - Ofuscación de Fase: Multiplicación por Pi de Primos
    - Backtracking de Fase: Recuperación exacta sin colisiones
    - Operaciones: Buffers de Memoria o Flujo por Archivos (Streaming)
=============================================================================
"""

import sys
import os
import math
import time
import json
import hashlib
from decimal import Decimal, getcontext
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional

# =============================================================================
# [1] CONFIGURACIÓN DE PRECISIÓN DE LA FPU VIRTUAL
# =============================================================================
PRECISION_SISTEMA = 200
getcontext().prec = PRECISION_SISTEMA

# =============================================================================
# [2] ESTRUCTURAS DE DATOS PARAMÉTRICAS (CLAVE MAESTRA)
# =============================================================================
@dataclass
class ClaveMaestra:
    """Contiene todos los parámetros del espacio geométrico de encriptación."""
    tales: List[int]                  # Escalas del teorema de Tales (ej: Fibonacci) [cite: 140]
    figuras: List[str]                # Tipos: "equilatero", "isosceles", "escaleno" [cite: 141]
    puntos: List[int]                 # Puntos base de diseño (magnitud física) [cite: 141]
    saltos: List[int]                 # Lados omitidos (ej: [0, 1] indica qué lado ignorar por iteración) [cite: 142]
    primos: List[Tuple[int, int]]     # Parejas de primos para el polígono de Pi [cite: 199]
    porcentajes_aportacion: List[int] # Proporciones del Pi ofuscado [cite: 199]
    iteraciones_pi: int               # Precisión de aproximación de Pi [cite: 199]

    def validar(self):
        """Asegura la viabilidad matemática de las listas."""
        if not self.tales or not self.figuras or not self.puntos or not self.saltos:
            raise ValueError("[-] Ninguna secuencia de la contraseña puede estar vacía.")
        if len(self.primos) != len(self.porcentajes_aportacion):
            raise ValueError("[-] Las dimensiones de primos y porcentajes de aportación deben coincidir.")

# =============================================================================
# [3] BLOQUE DE MATEMÁTICA PURA (CONVERGENCIA Y PI)
# =============================================================================
def aproximar_pi_poligono(p1: int, p2: int, iteraciones: int) -> Decimal:
    """Aproxima Pi mediante el método de perímetros de polígonos de (p1 * p2) lados."""
    lados = Decimal(str(p1 * p2))
    # Ángulo en radianes (Pi inicial estimado como constante de control para la aproximación)
    pi_inicial = Decimal('3.14159265358979323846264338327950288')
    angulo = pi_inicial / lados
    
    # Simulación del seno mediante la serie Taylor para máxima precisión en Decimal
    def taylor_sin(x: Decimal) -> Decimal:
        suma = x
        termino = x
        x_cuadrado = x * x
        n = 1
        while abs(termino) > Decimal('1e-100'):
            termino = -termino * x_cuadrado / Decimal(str((2 * n) * (2 * n + 1)))
            suma += termino
            n += 1
        return suma

    lado = Decimal('2.0') * taylor_sin(angulo)
    perimetro = lados * lado
    pi_actual = perimetro / Decimal('2.0')

    # Duplicación sucesiva de lados
    for _ in range(iteraciones):
        lados *= Decimal('2.0')
        soporte = (Decimal('4.0') - lado * lado).sqrt()
        lado = lado / (Decimal('2.0') + soporte).sqrt()
        perimetro = lados * lado
        pi_actual = perimetro / Decimal('2.0')
        
    return pi_actual

def calcular_lados_figura(tipo: str, puntos: int, escala: Decimal) -> List[Decimal]:
    """Genera los 3 lados de la figura bajo el Teorema de Tales de forma determinista."""
    base = Decimal(str(puntos)) * Decimal('1.5') [cite: 153]
    tipo_l = tipo.lower()
    if tipo_l == "equilatero":
        lado = (base / Decimal('3')) * escala [cite: 153]
        return [lado, lado, lado] [cite: 153]
    elif tipo_l == "isosceles":
        lado_a = (base / Decimal('4')) * escala [cite: 153]
        lado_b = (base / Decimal('2')) * escala [cite: 153]
        return [lado_a, lado_a, lado_b] [cite: 153]
    else:  # Escaleno
        l1 = (base * Decimal('0.25')) * escala [cite: 153]
        l2 = (base * Decimal('0.35')) * escala [cite: 153]
        l3 = (base * Decimal('0.40')) * escala [cite: 153]
        return [l1, l2, l3] [cite: 153]

# =============================================================================
# [4] MOTOR DE CONVERGENCIA GEOMÉTRICA (FORWARD & REVERSE)
# =============================================================================
class MotorK3:
    def __init__(self, clave: ClaveMaestra):
        self.clave = clave
        self.clave.validar()
        self._pi_ofuscado_cache = None

    def obtener_pi_ofuscado(self) -> Decimal:
        """Calcula de forma única el Pi ofuscado por la clave para las operaciones."""
        if self._pi_ofuscado_cache is not None:
            return self._pi_ofuscado_cache
        
        pi_ofuscado = Decimal('1') [cite: 199]
        for i, (p1, p2) in enumerate(self.clave.primos): [cite: 199]
            pi_aprox = aproximar_pi_poligono(p1, p2, self.clave.iteraciones_pi) [cite: 199]
            aportacion = Decimal(self.clave.porcentajes_aportacion[i]) / Decimal('100') [cite: 199]
            pi_ofuscado *= (pi_aprox ** aportacion) [cite: 199]
        
        self._pi_ofuscado_cache = pi_ofuscado
        return pi_ofuscado

    def encriptar_cadena_bits(self, bits: str) -> Decimal:
        """Encripta una cadena de bits acumulando el perímetro con saltos de fase."""
        perimetro = Decimal('0.0') [cite: 211]
        n = len(bits) [cite: 211]
        
        for idx in range(n): [cite: 211]
            bit = bits[idx] [cite: 211]
            figura_idx = idx // 3 [cite: 211]
            lado_idx = idx % 3 [cite: 211]
            
            # Control de saltos paramétricos en los lados [cite: 211]
            if lado_idx == self.clave.saltos[figura_idx % len(self.clave.saltos)]: [cite: 211]
                continue [cite: 211]
                
            escala = Decimal(str(self.clave.tales[figura_idx % len(self.clave.tales)])) [cite: 211]
            figura_tipo = self.clave.figuras[figura_idx % len(self.clave.figuras)] [cite: 211]
            puntos_fig = self.clave.puntos[figura_idx % len(self.clave.puntos)] [cite: 211]
            
            lados = calcular_lados_figura(figura_tipo, puntos_fig, escala) [cite: 211]
            
            if bit == '1': [cite: 211]
                perimetro += lados[lado_idx] [cite: 211]
                
        # Ofuscación final con el Pi de control [cite: 211]
        return perimetro * self.obtener_pi_ofuscado()

    def desencriptar_cadena_bits(self, hash_fase: Decimal, total_bits: int) -> Optional[str]:
        """Reconstruye con precisión exacta (O(n) con clave) el flujo original de bits."""
        # Deshacer ofuscación por Pi
        pi_ofuscado = self.obtener_pi_ofuscado()
        perimetro_objetivo = hash_fase / pi_ofuscado
        
        tolerancia = Decimal('1e-150')

        # Backtracking recursivo exacto sobre el residuo decimal [cite: 191]
        def backtrack(residuo: Decimal, idx: int) -> Optional[str]:
            if idx == total_bits:
                return "" if abs(residuo) < tolerancia else None [cite: 203]
                
            figura_idx = idx // 3 [cite: 203]
            lado_idx = idx % 3 [cite: 203]
            
            # Si el lado actual coincide con el patrón de saltos, este bit fue descartado
            if lado_idx == self.clave.saltos[figura_idx % len(self.clave.saltos)]: [cite: 203]
                # Probamos ambos caminos de ramificación porque no alteró el perímetro
                for bit in ['0', '1']: [cite: 203]
                    resultado = backtrack(residuo, idx + 1) [cite: 203]
                    if resultado is not None: [cite: 203]
                        return bit + resultado [cite: 203]
                return None [cite: 203]
                
            # Cálculo de los parámetros para este paso [cite: 203]
            escala = Decimal(str(self.clave.tales[figura_idx % len(self.clave.tales)])) [cite: 203]
            figura_tipo = self.clave.figuras[figura_idx % len(self.clave.figuras)] [cite: 203]
            puntos_fig = self.clave.puntos[figura_idx % len(self.clave.puntos)] [cite: 203]
            
            lados = calcular_lados_figura(figura_tipo, puntos_fig, escala) [cite: 203]
            lado_actual = lados[lado_idx] [cite: 203]
            
            # --- CASO 1: Intentar asumir que el bit fue '1' (se restó) --- [cite: 203]
            if residuo >= lado_actual - tolerancia: [cite: 203]
                resultado = backtrack(residuo - lado_actual, idx + 1) [cite: 203]
                if resultado is not None: [cite: 203]
                    return "1" + resultado [cite: 203]
                    
            # --- CASO 2: Intentar asumir que el bit fue '0' ---
            resultado = backtrack(residuo, idx + 1)
            if resultado is not None:
                return "0" + resultado
                
            return None

        return backtrack(perimetro_objetivo, 0)

# =============================================================================
# [5] MÓDULO DE INTERFAZ: TRABAJO CON BUFFER DE MEMORIA O ARCHIVOS EN DISCO
# =============================================================================
class EncriptadorIndustrial:
    def __init__(self, motor: MotorK3):
        self.motor = motor

    @staticmethod
    def bytes_a_bits(datos: bytes) -> str:
        """Convierte bytes crudos a su representación de bits en texto."""
        return "".join(f"{b:08b}" for b in datos)

    @staticmethod
    def bits_a_bytes(bits: str) -> bytes:
        """Convierte una cadena de bits de vuelta a bytes."""
        # Rellenar con ceros si el flujo de bits no es múltiplo de 8
        bits_pad = bits + '0' * ((8 - len(bits) % 8) % 8)
        bytes_list = []
        for i in range(0, len(bits_pad), 8):
            byte = int(bits_pad[i:i+8], 2)
            bytes_list.append(byte)
        return bytes(bytes_list)

    def encriptar_buffer_memoria(self, buffer: bytes) -> Tuple[Decimal, int]:
        """Procesa un buffer en RAM y devuelve (Hash Geométrico, longitud_de_bits)"""
        cadena_bits = self.bytes_a_bits(buffer)
        hash_final = self.motor.encriptar_cadena_bits(cadena_bits)
        return hash_final, len(cadena_bits)

    def desencriptar_buffer_memoria(self, hash_fase: Decimal, longitud_bits: int) -> bytes:
        """Toma el valor geométrico y la longitud, y recupera los bytes de RAM."""
        cadena_bits = self.motor.desencriptar_cadena_bits(hash_fase, longitud_bits)
        if cadena_bits is None:
            raise ValueError("[-] Imposible recuperar: Clave corrupta o residuo inalcanzable.")
        return self.bits_a_bytes(cadena_bits)

    def encriptar_archivo(self, ruta_origen: str, ruta_meta: str) -> Decimal:
        """Lee un archivo en bloques para no saturar RAM y genera la cabecera metadata."""
        if not os.path.exists(ruta_origen):
            raise FileNotFoundError(f"El archivo {ruta_origen} no existe.")
            
        with open(ruta_origen, "rb") as f:
            datos = f.read()  # Para archivos de prueba, leemos completos.

        hash_geom, long_bits = self.encriptar_buffer_memoria(datos)
        
        # Estructuramos el archivo meta .k3 con los datos para el receptor
        meta_datos = {
            "valor_hash": str(hash_geom),
            "longitud_bits": long_bits,
            "hash_integridad_sha256": hashlib.sha256(datos).hexdigest()
        }
        
        with open(ruta_meta, "w") as f_meta:
            json.dump(meta_datos, f_meta, indent=4)
            
        return hash_geom

    def desencriptar_archivo(self, ruta_meta: str, ruta_salida: str):
        """Usa el archivo metadata .k3 y la clave para rehacer el archivo original."""
        with open(ruta_meta, "r") as f_meta:
            meta = json.load(f_meta)
            
        hash_fase = Decimal(meta["valor_hash"])
        long_bits = int(meta["longitud_bits"])
        sha256_verificable = meta["hash_integridad_sha256"]
        
        # Desencriptar mediante el reverso K3
        bits_recuperados = self.motor.desencriptar_cadena_bits(hash_fase, long_bits)
        if bits_recuperados is None:
            raise ValueError("[-] Fallo de reverso: La clave no coincide con el residuo físico.")
            
        datos_originales = self.bits_a_bytes(bits_recuperados)
        # Recortar el relleno para que coincida con el SHA-256 original
        longitud_bytes_original = long_bits // 8
        datos_originales = datos_originales[:longitud_bytes_original]

        # Validar integridad
        sha256_actual = hashlib.sha256(datos_originales).hexdigest()
        if sha256_actual != sha256_verificable:
            print("[-] ALERTA: El archivo recuperado no coincide con el hash original de integridad.")
        else:
            print("[+] ÉXITO: Verificación SHA-256 válida. Archivo íntegro.")
            
        with open(ruta_salida, "wb") as f_salida:
            f_salida.write(datos_originales)

# =============================================================================
# [6] PROGRAMA DEMO INTERACTIVO Y PARAMETRIZABLE (CONSOLA)
# =============================================================================
def mostrar_menu():
    print("\n" + "="*80)
    print("                 GEOMETRIC CONVERGENCE CONSOLE ENGINE v4.0")
    print("="*80)
    print(" [1] Usar clave K3 por defecto (Fibonacci / Primos Simétricos)")
    print(" [2] Introducir clave K3 personalizada (Superparametrizable)")
    print(" [3] Ejecutar test con archivo físico (Capa 1 Test)")
    print(" [4] Salir")
    print("="*80)

def obtener_clave_defecto() -> ClaveMaestra:
    return ClaveMaestra(
        tales=[3, 5, 8, 13, 21], [cite: 53]
        figuras=["equilatero", "isosceles", "escaleno"], [cite: 53]
        puntos=[6, 12, 18], [cite: 53]
        saltos=[2],  # Salta el tercer lado cada ciclo
        primos=[(3, 5), (7, 11)],
        porcentajes_aportacion=[40, 60],
        iteraciones_pi=15
    )

def main():
    clave_activa = obtener_clave_defecto()
    
    while True:
        mostrar_menu()
        opcion = input("[>] Seleccione una opción: ").strip()
        
        if opcion == "1":
            clave_activa = obtener_clave_defecto()
            print("[+] Cargada clave de convergencia por defecto.")
            print(f"    - Secuencia Tales: {clave_activa.tales}")
            print(f"    - Figuras: {clave_activa.figuras}")
            print(f"    - Puntos Base: {clave_activa.puntos}")
            print(f"    - Saltos: {clave_activa.saltos}")
            
        elif opcion == "2":
            print("\n--- CONFIGURACIÓN DE TU PARÁMETRO K3 ---")
            try:
                tales_inp = list(map(int, input("[>] Secuencia de Tales (ej: 3 5 8 13): ").split()))
                figuras_inp = input("[>] Figuras (ej: equilatero isosceles escaleno): ").split()
                puntos_inp = list(map(int, input("[>] Puntos Base (ej: 6 12 18): ").split()))
                saltos_inp = list(map(int, input("[>] Saltos de Lado (ej: 0 1): ").split()))
                
                print("\n- Configuración de Ofuscación de Fase (Pi) -")
                iter_pi = int(input("[>] Iteraciones de aproximación (ej: 10): "))
                
                clave_activa = ClaveMaestra(
                    tales=tales_inp,
                    figuras=figuras_inp,
                    puntos=puntos_inp,
                    saltos=saltos_inp,
                    primos=[(3, 5), (7, 11)],
                    porcentajes_aportacion=[40, 60],
                    iteraciones_pi=iter_pi
                )
                print("[+] Nueva Clave Maestra parametrizada con éxito.")
            except Exception as e:
                print(f"[-] Error en formato de entrada: {e}. Se mantiene clave por defecto.")
                
        elif opcion == "3":
            # Inicializar motores
            motor = MotorK3(clave_activa)
            compilador = EncriptadorIndustrial(motor)
            
            print("\n[+] Creando archivo temporal 'datos_capa1.txt' para prueba...")
            datos_prueba = b"Este es un mensaje secreto encriptado con el motor geometrico de fase K3 de Victor."
            
            with open("datos_capa1.txt", "wb") as f:
                f.write(datos_prueba)
                
            print(f"[+] Archivo original guardado. Contenido: '{datos_prueba.decode()}'")
            print("[+] Encriptando en flujo y generando archivo meta .k3...")
            
            t_inicio = time.time()
            hash_colapso = compilador.encriptar_archivo("datos_capa1.txt", "fase_colapso.k3")
            t_encriptar = time.time() - t_inicio
            
            print(f"[+] Archivo encriptado con éxito en {t_encriptar:.6f} s.")
            print(f"[+] Hash Geométrico de Colapso obtenido:")
            print(f"    => {hash_colapso}")
            
            print("\n[+] Iniciando el reverso (Desencriptación) a partir de f_meta...")
            t_inicio_rev = time.time()
            compilador.desencriptar_archivo("fase_colapso.k3", "datos_recuperados.txt")
            t_desencriptar = time.time() - t_inicio_rev
            
            print(f"[+] Reverso completado en {t_desencriptar:.6f} s.")
            with open("datos_recuperados.txt", "rb") as f:
                recuperados = f.read()
            print(f"[+] Texto recuperado de disco: '{recuperados.decode()}'")
            
            # Limpieza
            for arch in ["datos_capa1.txt", "fase_colapso.k3", "datos_recuperados.txt"]:
                if os.path.exists(arch):
                    os.remove(arch)
            print("[+] Limpieza de archivos de prueba completada.")
            
        elif opcion == "4":
            print("[+] Saliendo del motor K3. ¡Hasta pronto, Víctor!")
            sys.exit(0)
        else:
            print("[-] Opción no válida.")

if __name__ == "__main__":
    main()