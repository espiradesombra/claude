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
[+] SISTEMA INTEGRAL DE CONVERGENCIA GEOMÉTRICA K3 - CON RUIDO DE FASE Y RED
[+] AUTOR: Víctor Manzanares Alberola
[+] VERSIÓN: 2.0 - Industrial (con Aleatorovix y Simulación de Red)
=============================================================================
"""

import sys
import os
import math
import time
import json
import hashlib
import random
import socket
import threading
from decimal import Decimal, getcontext
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

# =============================================================================
# [1] CONFIGURACIÓN DE PRECISIÓN
# =============================================================================
PRECISION_SISTEMA = 200
getcontext().prec = PRECISION_SISTEMA

# =============================================================================
# [2] ESTRUCTURAS DE DATOS PARAMÉTRICAS (CLAVE MAESTRA)
# =============================================================================
@dataclass
class ClaveMaestra:
    tales: List[int]
    figuras: List[str]
    puntos: List[int]
    saltos: List[int]
    primos: List[Tuple[int, int]]
    porcentajes_aportacion: List[int]
    iteraciones_pi: int

    def validar(self):
        if not self.tales or not self.figuras or not self.puntos or not self.saltos:
            raise ValueError("[-] Ninguna secuencia de la contraseña puede estar vacía.")
        if len(self.primos) != len(self.porcentajes_aportacion):
            raise ValueError("[-] Las dimensiones de primos y porcentajes de aportación deben coincidir.")

# =============================================================================
# [3] MATEMÁTICA PURA: PI Y LADOS GEOMÉTRICOS
# =============================================================================
def aproximar_pi_poligono(p1: int, p2: int, iteraciones: int) -> Decimal:
    lados = Decimal(str(p1 * p2))
    pi_inicial = Decimal('3.14159265358979323846264338327950288')
    angulo = pi_inicial / lados

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

    for _ in range(iteraciones):
        lados *= Decimal('2.0')
        soporte = (Decimal('4.0') - lado * lado).sqrt()
        lado = lado / (Decimal('2.0') + soporte).sqrt()
        perimetro = lados * lado
        pi_actual = perimetro / Decimal('2.0')

    return pi_actual

def calcular_lados_figura(tipo: str, puntos: int, escala: Decimal) -> List[Decimal]:
    base = Decimal(str(puntos)) * Decimal('1.5')
    tipo_l = tipo.lower()
    if tipo_l == "equilatero":
        lado = (base / Decimal('3')) * escala
        return [lado, lado, lado]
    elif tipo_l == "isosceles":
        lado_a = (base / Decimal('4')) * escala
        lado_b = (base / Decimal('2')) * escala
        return [lado_a, lado_a, lado_b]
    else:  # escaleno
        l1 = (base * Decimal('0.25')) * escala
        l2 = (base * Decimal('0.35')) * escala
        l3 = (base * Decimal('0.40')) * escala
        return [l1, l2, l3]

# =============================================================================
# [4] GENERADOR DE RUIDO DE FASE: ALEATOROVIX MEJORADO
# =============================================================================
class Aleatorovix:
    """
    Generador de flujo de bits caótico determinista basado en la semilla de fase.
    La semilla se deriva del hash geométrico, lo que hace que el flujo sea único para cada mensaje.
    """
    def __init__(self, semilla: Decimal):
        self.estado = semilla
        self.factor_msl = Decimal('42.123456789')
        self.factor_lsl = Decimal('0.000123456789')

    def siguiente_bit(self) -> int:
        # Sincronización caótica: x_{n+1} = (x_n * factor_msl + cos(x_n * factor_lsl)) mod 1
        fase = (self.estado * self.factor_msl) + self.estado.cos() * self.factor_lsl
        self.estado = fase - fase.to_integral_value()
        if self.estado < 0:
            self.estado += Decimal('1.0')
        return 1 if self.estado > Decimal('0.5') else 0

    def generar_bytes(self, longitud: int) -> bytes:
        bytes_out = bytearray()
        for _ in range(longitud):
            byte = 0
            for bit_pos in range(8):
                byte = (byte << 1) | self.siguiente_bit()
            bytes_out.append(byte)
        return bytes(bytes_out)

# =============================================================================
# [5] MOTOR DE CONVERGENCIA GEOMÉTRICA CON RUIDO DE FASE
# =============================================================================
class MotorK3:
    def __init__(self, clave: ClaveMaestra):
        self.clave = clave
        self.clave.validar()
        self._pi_ofuscado_cache = None

    def obtener_pi_ofuscado(self) -> Decimal:
        if self._pi_ofuscado_cache is not None:
            return self._pi_ofuscado_cache
        pi_ofuscado = Decimal('1')
        for i, (p1, p2) in enumerate(self.clave.primos):
            pi_aprox = aproximar_pi_poligono(p1, p2, self.clave.iteraciones_pi)
            aportacion = Decimal(self.clave.porcentajes_aportacion[i]) / Decimal('100')
            pi_ofuscado *= (pi_aprox ** aportacion)
        self._pi_ofuscado_cache = pi_ofuscado
        return pi_ofuscado

    def encriptar_cadena_bits(self, bits: str) -> Decimal:
        perimetro = Decimal('0.0')
        n = len(bits)

        for idx in range(n):
            bit = bits[idx]
            figura_idx = idx // 3
            lado_idx = idx % 3

            if lado_idx == self.clave.saltos[figura_idx % len(self.clave.saltos)]:
                continue

            escala = Decimal(str(self.clave.tales[figura_idx % len(self.clave.tales)]))
            figura_tipo = self.clave.figuras[figura_idx % len(self.clave.figuras)]
            puntos_fig = self.clave.puntos[figura_idx % len(self.clave.puntos)]

            lados = calcular_lados_figura(figura_tipo, puntos_fig, escala)

            if bit == '1':
                perimetro += lados[lado_idx]

        return perimetro * self.obtener_pi_ofuscado()

    def desencriptar_cadena_bits(self, hash_fase: Decimal, total_bits: int) -> Optional[str]:
        pi_ofuscado = self.obtener_pi_ofuscado()
        perimetro_objetivo = hash_fase / pi_ofuscado
        tolerancia = Decimal('1e-150')

        def backtrack(residuo: Decimal, idx: int) -> Optional[str]:
            if idx == total_bits:
                return "" if abs(residuo) < tolerancia else None

            figura_idx = idx // 3
            lado_idx = idx % 3

            if lado_idx == self.clave.saltos[figura_idx % len(self.clave.saltos)]:
                for bit in ['0', '1']:
                    resultado = backtrack(residuo, idx + 1)
                    if resultado is not None:
                        return bit + resultado
                return None

            escala = Decimal(str(self.clave.tales[figura_idx % len(self.clave.tales)]))
            figura_tipo = self.clave.figuras[figura_idx % len(self.clave.figuras)]
            puntos_fig = self.clave.puntos[figura_idx % len(self.clave.puntos)]

            lados = calcular_lados_figura(figura_tipo, puntos_fig, escala)
            lado_actual = lados[lado_idx]

            if residuo >= lado_actual - tolerancia:
                resultado = backtrack(residuo - lado_actual, idx + 1)
                if resultado is not None:
                    return "1" + resultado

            resultado = backtrack(residuo, idx + 1)
            if resultado is not None:
                return "0" + resultado

            return None

        return backtrack(perimetro_objetivo, 0)

    def generar_flujo_ruido(self, hash_fase: Decimal, longitud: int) -> bytes:
        """Genera un flujo de bytes caótico derivado del hash de fase."""
        aleatorovix = Aleatorovix(hash_fase)
        return aleatorovix.generar_bytes(longitud)

# =============================================================================
# [6] INTERFAZ DE RED Y PAQUETES
# =============================================================================
class PaqueteK3:
    """Representa un paquete de datos con su hash de fase y metadatos."""
    def __init__(self, datos: bytes, hash_fase: Decimal, longitud_bits: int):
        self.datos = datos
        self.hash_fase = hash_fase
        self.longitud_bits = longitud_bits

    def serializar(self) -> bytes:
        """Convierte el paquete a bytes para transmisión por red."""
        meta = {
            "hash_fase": str(self.hash_fase),
            "longitud_bits": self.longitud_bits
        }
        meta_json = json.dumps(meta).encode('utf-8')
        meta_len = len(meta_json).to_bytes(4, 'big')
        return meta_len + meta_json + self.datos

    @staticmethod
    def deserializar(paquete_bytes: bytes) -> 'PaqueteK3':
        """Reconstruye un paquete desde bytes recibidos."""
        meta_len = int.from_bytes(paquete_bytes[:4], 'big')
        meta_json = paquete_bytes[4:4+meta_len].decode('utf-8')
        meta = json.loads(meta_json)
        datos = paquete_bytes[4+meta_len:]
        return PaqueteK3(datos, Decimal(meta["hash_fase"]), meta["longitud_bits"])

class ServidorK3:
    """Simula un servidor que recibe paquetes y los desencripta con la clave maestra."""
    def __init__(self, motor: MotorK3, puerto: int = 9999):
        self.motor = motor
        self.puerto = puerto
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.socket.bind(('localhost', puerto))
        self.socket.listen(1)

    def manejar_cliente(self, conn):
        print("[*] Cliente conectado. Recibiendo paquete...")
        datos_completos = b""
        while True:
            fragmento = conn.recv(4096)
            if not fragmento:
                break
            datos_completos += fragmento
        paquete = PaqueteK3.deserializar(datos_completos)

        print(f"[*] Paquete recibido. Hash: {paquete.hash_fase}, Longitud: {paquete.longitud_bits} bits.")
        try:
            bits_recuperados = self.motor.desencriptar_cadena_bits(paquete.hash_fase, paquete.longitud_bits)
            if bits_recuperados is None:
                print("[!] Error: No se pudo desencriptar el paquete.")
                conn.send(b"ERROR")
            else:
                # Recuperar bytes originales
                bytes_recuperados = bytes(int(bits_recuperados[i:i+8], 2) for i in range(0, len(bits_recuperados), 8))
                print(f"[+] Datos recuperados: {bytes_recuperados}")
                conn.send(b"OK")
        except Exception as e:
            print(f"[-] Error durante la desencriptación: {e}")
            conn.send(b"ERROR")
        conn.close()

    def iniciar(self):
        print(f"[+] Servidor K3 escuchando en el puerto {self.puerto}...")
        while True:
            conn, addr = self.socket.accept()
            threading.Thread(target=self.manejar_cliente, args=(conn,)).start()

class ClienteK3:
    """Simula un cliente que envía un paquete encriptado a un servidor."""
    def __init__(self, motor: MotorK3, servidor_host: str = 'localhost', servidor_puerto: int = 9999):
        self.motor = motor
        self.servidor_host = servidor_host
        self.servidor_puerto = servidor_puerto

    def enviar_datos(self, datos: bytes):
        # Convertir a bits
        bits = "".join(f"{b:08b}" for b in datos)
        longitud_bits = len(bits)

        # Generar hash de fase
        hash_fase = self.motor.encriptar_cadena_bits(bits)

        # Generar flujo de ruido (para simular tráfico adicional en red)
        ruido = self.motor.generar_flujo_ruido(hash_fase, len(datos))

        # Crear paquete (datos originales + metadatos)
        paquete = PaqueteK3(datos, hash_fase, longitud_bits)
        paquete_bytes = paquete.serializar()

        # Conectar al servidor y enviar
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect((self.servidor_host, self.servidor_puerto))
        sock.sendall(paquete_bytes)
        respuesta = sock.recv(1024)
        print(f"[+] Respuesta del servidor: {respuesta.decode('utf-8')}")
        sock.close()

# =============================================================================
# [7] PROGRAMA PRINCIPAL (DEMO CON RED)
# =============================================================================
def obtener_clave_defecto() -> ClaveMaestra:
    return ClaveMaestra(
        tales=[3, 5, 8, 13, 21],
        figuras=["equilatero", "isosceles", "escaleno"],
        puntos=[6, 12, 18],
        saltos=[2],
        primos=[(3, 5), (7, 11)],
        porcentajes_aportacion=[40, 60],
        iteraciones_pi=15
    )

def main():
    print("=" * 80)
    print("   SISTEMA K3 CON RUIDO DE FASE Y SIMULACIÓN DE RED")
    print("=" * 80)

    clave = obtener_clave_defecto()
    motor = MotorK3(clave)

    print("[+] Motor K3 inicializado con clave por defecto.")
    print("[*] Seleccione una opción:")
    print("  1. Ejecutar servidor de prueba (escucha paquetes)")
    print("  2. Enviar paquete de prueba a servidor local")
    print("  3. Encriptar archivo (flujo local)")
    print("  4. Salir")

    opcion = input("[>] Opción: ").strip()

    if opcion == "1":
        servidor = ServidorK3(motor)
        servidor.iniciar()

    elif opcion == "2":
        mensaje = input("[>] Introduzca el mensaje a enviar: ").encode('utf-8')
        cliente = ClienteK3(motor)
        cliente.enviar_datos(mensaje)

    elif opcion == "3":
        ruta = input("[>] Ruta del archivo a encriptar: ").strip()
        try:
            with open(ruta, "rb") as f:
                datos = f.read()
            bits = "".join(f"{b:08b}" for b in datos)
            hash_fase = motor.encriptar_cadena_bits(bits)
            print(f"[+] Hash de fase generado: {hash_fase}")
            print("[+] Archivo encriptado (solo hash generado).")
        except Exception as e:
            print(f"[-] Error: {e}")

    elif opcion == "4":
        print("[+] Saliendo...")
        sys.exit(0)
    else:
        print("[-] Opción no válida.")

if __name__ == "__main__":
    main()