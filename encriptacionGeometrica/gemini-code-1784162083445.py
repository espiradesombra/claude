import ctypes
import os

class EncriptadorIndustrial:
    def __init__(self, base, rel):
        self.lib = ctypes.CDLL("./k3_core.so")
        self.motor = (ctypes.c_uint64 * 2)(base, rel)
        
    def procesar(self, archivo_ruta):
        with open(archivo_ruta, "rb") as f:
            datos = bytearray(f.read())
        
        # Llamada al núcleo C
        ptr = (ctypes.c_ubyte * len(datos)).from_buffer(datos)
        self.lib.k3_cifrar_bloque(ptr, len(datos), self.motor)
        
        print(f"[+] Bloque procesado con pinza {self.motor[0]}x{self.motor[1]}")