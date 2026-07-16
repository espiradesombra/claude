import ctypes
import os
import subprocess
import getpass

# Compilación automática del núcleo
subprocess.run(["gcc", "-O3", "-shared", "-fPIC", "-o", "k3_core.so", "k3_core.c"])

class EncriptadorIndustrial:
    def __init__(self, base, rel):
        self.lib = ctypes.CDLL("./k3_core.so")
        self.motor = (ctypes.c_uint64 * 2)(base, rel)
        
    def procesar(self, ruta_archivo):
        with open(ruta_archivo, "rb") as f:
            datos = bytearray(f.read())
        
        # Inyectar al motor C
        ptr = (ctypes.c_ubyte * len(datos)).from_buffer(datos)
        self.lib.k3_cifrar_bloque(ptr, len(datos), self.motor)
        
        # Eliminar archivo original tras cifrar
        os.remove(ruta_archivo)
        print(f"[+] Archivo {ruta_archivo} cifrado y rastro eliminado.")

if __name__ == "__main__":
    print("--- MOTOR DE CONVERGENCIA K3 INDUSTRIAL ---")
    pwd = getpass.getpass("Clave de fase: ")
    engine = EncriptadorIndustrial(base=33, rel=1)
    engine.procesar("datos_capa1.txt")