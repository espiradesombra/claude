import sys
from k3_core.cipher import EncriptadorIndustrial

def main():
    print("--- PROYECTO 33x1: MOTOR DE CONVERGENCIA INDUSTRIAL ---")
    
    # 1. Configuración de parámetros con valores de fábrica (33:1)
    # El usuario puede aceptar el 33:1 o personalizarlo.
    config = {
        "factor_base": 33,
        "factor_relacion": 1
    }
    
    opcion = input(f"[?] Configuración actual {config['factor_base']}x{config['factor_relacion']}. ¿Cambiar? (s/n): ")
    if opcion.lower() == 's':
        config['factor_base'] = int(input("Introduce nuevo factor base: "))
        config['factor_relacion'] = int(input("Introduce nuevo factor relación: "))
    
    # 2. Captura de contraseña (scanf)
    import getpass
    pwd = getpass.getpass("[!] Introduce tu clave maestra: ")
    
    # 3. Ejecución
    # Aquí instanciamos el motor usando el modo 33x1
    engine = EncriptadorIndustrial(base=config['factor_base'], rel=config['factor_relacion'])
    
    # ... lógica de comandos (cifrar/descifrar) ...
    print(f"[+] Sistema operando en modo {config['factor_base']}x{config['factor_relacion']}")