# aleatorovix_final.py
from k3_industrial import EncriptadorIndustrial
import getpass

def ejecutar():
    print("--- PROYECTO 33x1: MOTOR DE CONVERGENCIA INDUSTRIAL ---")
    
    # 1. Configuración de parámetros
    factor_base = 33
    factor_relacion = 1
    
    # 2. Captura segura
    pwd = getpass.getpass("[!] Introduce tu clave maestra: ")
    
    # 3. Lanzamiento
    engine = EncriptadorIndustrial(base=factor_base, rel=factor_relacion)
    engine.procesar("datos_capa1.txt")
    print("[+] Sistema operando. Memoria purgada.")

if __name__ == "__main__":
    ejecutar()