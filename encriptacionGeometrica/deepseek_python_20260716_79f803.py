#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess
import sys
import os
import getpass
import json
from decimal import Decimal

def compilar_motor():
    """Compila el motor en C si no existe el binario."""
    if not os.path.exists("./k3_engine"):
        print("[*] Compilando el motor K3 en C...")
        try:
            subprocess.run(["gcc", "-O3", "-o", "k3_engine", "k3_engine.c", "-lm"], check=True)
            print("[+] Motor compilado con éxito.")
        except subprocess.CalledProcessError as e:
            print(f"[-] Error al compilar: {e}")
            sys.exit(1)

def ejecutar_motor(comando, semilla, datos, clave_params):
    """Ejecuta el binario del motor y captura la salida."""
    # Crear un archivo temporal con los parámetros de la clave
    with open("clave_temp.json", "w") as f:
        json.dump(clave_params, f)
    
    # Construir el comando del motor
    # Nota: El motor en C debe ser modificado para leer desde stdin o archivos.
    # Por ahora, pasamos los datos como argumento.
    proceso = subprocess.Popen(
        ["./k3_engine", comando, str(semilla), datos, "clave_temp.json"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    salida, error = proceso.communicate()
    
    if proceso.returncode != 0:
        print(f"[-] Error en el motor: {error.decode()}")
        sys.exit(1)
    
    return salida.decode()

def main():
    print("=== PROYECTO 33x1: CONTROLADOR INDUSTRIAL K3 ===")
    
    # 1. Compilar el motor si es necesario
    compilar_motor()
    
    # 2. Configuración de parámetros (por defecto 33x1)
    config = {"base": 33, "rel": 1}
    print(f"[!] Modo de convergencia por defecto: {config['base']}x{config['rel']}")
    
    opcion = input("[?] ¿Cambiar los factores de convergencia? (s/n): ").lower()
    if opcion == 's':
        config['base'] = int(input("Nuevo factor base: "))
        config['rel'] = int(input("Nuevo factor relación: "))
    
    # 3. Captura de la Clave Maestra (con getpass para seguridad)
    clave_maestra = getpass.getpass("[?] Introduce la Clave Maestra de Convergencia: ")
    
    # 4. Derivación de la semilla a partir de la clave maestra (usando SHA-256)
    import hashlib
    semilla = int(hashlib.sha256(clave_maestra.encode()).hexdigest(), 16) % 65536
    
    # 5. Selección de la operación
    print("\n[1] Cifrar un archivo")
    print("[2] Descifrar un archivo")
    operacion = input("[?] Selecciona una opción: ").strip()
    
    if operacion == '1':
        ruta_archivo = input("[?] Ruta del archivo a cifrar: ")
        if not os.path.exists(ruta_archivo):
            print("[-] Archivo no encontrado.")
            return
        
        # Parámetros de la clave K3 (esto debería derivarse de la clave maestra)
        clave_params = {
            "tales": [5, 8, 13, 21, 34],
            "figuras": ["escaleno", "isosceles", "equilatero"],
            "puntos": [8, 16, 24],
            "saltos": [2, 0, 1],
            "primos": [[11, 13], [17, 19]],
            "porcentajes_aportacion": [40, 60],
            "iteraciones_pi": 20
        }
        
        # Ejecutar el motor de cifrado
        salida = ejecutar_motor("cifrar", semilla, ruta_archivo, clave_params)
        print(f"[+] Archivo cifrado con éxito. Hash de fase: {salida}")
        
    elif operacion == '2':
        ruta_hash = input("[?] Ruta del archivo .k3 (hash de fase): ")
        if not os.path.exists(ruta_hash):
            print("[-] Archivo no encontrado.")
            return
        
        # Aquí iría la lógica para descifrar usando el motor
        print("[+] Funcionalidad de descifrado en desarrollo...")
    
    else:
        print("[-] Opción no válida.")

if __name__ == "__main__":
    main()