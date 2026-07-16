import ctypes
import os
from ctypes import c_uint16, c_uint8, c_uint32, c_int, c_longdouble, c_bool, POINTER, Structure, Array, byref

# Cargar la biblioteca
_lib = None
_lib_path = os.path.join(os.path.dirname(__file__), '..', 'build', 'libk3.so')
if not os.path.exists(_lib_path):
    # Intentar en el directorio actual
    _lib_path = './libk3.so'
_lib = ctypes.CDLL(_lib_path)

# Definir estructuras
class ClaveK3(Structure):
    _fields_ = [
        ("tales", c_int * 5),
        ("tales_count", c_int),
        ("figuras", (c_char * 16) * 3),
        ("figuras_count", c_int),
        ("puntos", c_int * 3),
        ("puntos_count", c_int),
        ("saltos", c_int * 3),
        ("saltos_count", c_int),
        ("primos", (c_int * 2) * 2),
        ("porcentajes_aportacion", c_int * 2),
        ("iteraciones_pi", c_int),
    ]

# Prototipos de funciones
_lib.k3_encriptar_semilla.argtypes = [c_uint16, POINTER(ClaveK3)]
_lib.k3_encriptar_semilla.restype = c_longdouble

_lib.k3_desencriptar_semilla.argtypes = [c_longdouble, POINTER(ClaveK3)]
_lib.k3_desencriptar_semilla.restype = c_int

_lib.k3_aplicar_mascara.argtypes = [POINTER(c_uint8), POINTER(c_uint8), c_uint32, c_uint16]
_lib.k3_aplicar_mascara.restype = None

# Wrapper en Python
class K3Engine:
    def __init__(self, clave: dict):
        self.clave_c = ClaveK3()
        # Rellenar la estructura C desde el diccionario
        self.clave_c.tales = (c_int * 5)(*clave['tales'])
        self.clave_c.tales_count = len(clave['tales'])
        for i, f in enumerate(clave['figuras']):
            self.clave_c.figuras[i] = f.encode('utf-8')
        self.clave_c.figuras_count = len(clave['figuras'])
        self.clave_c.puntos = (c_int * 3)(*clave['puntos'])
        self.clave_c.puntos_count = len(clave['puntos'])
        self.clave_c.saltos = (c_int * 3)(*clave['saltos'])
        self.clave_c.saltos_count = len(clave['saltos'])
        self.clave_c.primos = (c_int * 2)(*clave['primos'][0]), (c_int * 2)(*clave['primos'][1])
        self.clave_c.porcentajes_aportacion = (c_int * 2)(*clave['porcentajes_aportacion'])
        self.clave_c.iteraciones_pi = clave['iteraciones_pi']

    def encriptar(self, semilla: int) -> float:
        """Colapsa una semilla de 16 bits en un número de fase (hash)."""
        return _lib.k3_encriptar_semilla(c_uint16(semilla), byref(self.clave_c))

    def desencriptar(self, hash_fase: float) -> int:
        """Recupera la semilla original desde el hash de fase."""
        return _lib.k3_desencriptar_semilla(c_longdouble(hash_fase), byref(self.clave_c))

    def aplicar_mascara(self, datos: bytes, semilla: int) -> bytes:
        """Aplica la máscara lila (XOR con flujo caótico) a los datos."""
        longitud = len(datos)
        origen = (c_uint8 * longitud).from_buffer_copy(datos)
        destino = (c_uint8 * longitud)()
        _lib.k3_aplicar_mascara(origen, destino, c_uint32(longitud), c_uint16(semilla))
        return bytes(destino)