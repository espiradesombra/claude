# k3_tool.py
import getpass
import sys
from k3_engine import MotorK3

def main():
    print("=== PROYECTO 33x1: CLI DE CONVERGENCIA INDUSTRIAL ===")
    
    # Valores por defecto del sistema
    config = {"base": 33, "rel": 1}
    
    # Interfaz de usuario para confirmación de parámetros
    print(f"[!] Configuración actual: {config['base']}x{config['rel']}")
    opcion = input("[?] ¿Deseas modificar los factores? (s/n): ").lower()
    
    if opcion == 's':
        config['base'] = int(input("Introduce nuevo factor base: "))
        config['rel'] = int(input("Introduce nuevo factor relación: "))
    
    # Scanf de seguridad (no visible en terminal)
    clave = getpass.getpass("[?] Introduce la Clave Maestra: ")
    
    # Instanciación del motor
    motor = MotorK3(config['base'], config['rel'])
    
    # Ejemplo de uso: Cifrar archivo
    ruta = input("[?] Ruta del archivo a encriptar: ")
    try:
        with open(ruta, "rb") as f:
            datos = f.read()
        
        print(f"[+] Iniciando convergencia en modo {config['base']}x{config['rel']}...")
        hash_salida = motor.encriptar(datos)
        
        with open(f"{ruta}.k3", "w") as f_out:
            f_out.write(hash_salida)
            
        print(f"[+] Archivo encriptado exitosamente: {ruta}.k3")
        
    except FileNotFoundError:
        print("[!] Error: Archivo no encontrado.")

if __name__ == "__main__":
    main()