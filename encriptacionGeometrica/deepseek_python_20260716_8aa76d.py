import ctypes
import os

_lib = ctypes.CDLL(os.path.join(os.path.dirname(__file__), 'libk3.so'))

class ClaveK3(ctypes.Structure):
    _fields_ = [
        ("tales", ctypes.c_int * 5),
        ("tales_count", ctypes.c_int),
        ("figuras", (ctypes.c_char * 16) * 3),
        ("figuras_count", ctypes.c_int),
        ("puntos", ctypes.c_int * 3),
        ("puntos_count", ctypes.c_int),
        ("saltos", ctypes.c_int * 3),
        ("saltos_count", ctypes.c_int),
        ("primos", (ctypes.c_int * 2) * 2),
        ("porcentajes_aportacion", ctypes.c_int * 2),
        ("iteraciones_pi", ctypes.c_int),
    ]

_lib.k3_encriptar.argtypes = [ctypes.c_uint16, ctypes.POINTER(ClaveK3)]
_lib.k3_encriptar.restype = ctypes.c_longdouble
_lib.k3_desencriptar.argtypes = [ctypes.c_longdouble, ctypes.POINTER(ClaveK3)]
_lib.k3_desencriptar.restype = ctypes.c_int
_lib.k3_mascara.argtypes = [ctypes.POINTER(ctypes.c_uint8), ctypes.POINTER(ctypes.c_uint8), ctypes.c_uint32, ctypes.c_uint16]
_lib.k3_mascara.restype = None

class K3:
    def __init__(self, clave_dict):
        self.clave = ClaveK3()
        self.clave.tales = (ctypes.c_int * 5)(*clave_dict['tales'])
        self.clave.tales_count = len(clave_dict['tales'])
        for i, f in enumerate(clave_dict['figuras']):
            self.clave.figuras[i] = f.encode('utf-8')[:15]
        self.clave.figuras_count = len(clave_dict['figuras'])
        self.clave.puntos = (ctypes.c_int * 3)(*clave_dict['puntos'])
        self.clave.puntos_count = len(clave_dict['puntos'])
        self.clave.saltos = (ctypes.c_int * 3)(*clave_dict['saltos'])
        self.clave.saltos_count = len(clave_dict['saltos'])
        self.clave.primos[0] = (ctypes.c_int * 2)(*clave_dict['primos'][0])
        self.clave.primos[1] = (ctypes.c_int * 2)(*clave_dict['primos'][1])
        self.clave.porcentajes_aportacion = (ctypes.c_int * 2)(*clave_dict['porcentajes_aportacion'])
        self.clave.iteraciones_pi = clave_dict['iteraciones_pi']

    def encriptar(self, semilla):
        return _lib.k3_encriptar(ctypes.c_uint16(semilla), ctypes.byref(self.clave))

    def desencriptar(self, hash_val):
        return _lib.k3_desencriptar(ctypes.c_longdouble(hash_val), ctypes.byref(self.clave))

    def mascara(self, data, semilla):
        buf = (ctypes.c_uint8 * len(data)).from_buffer_copy(data)
        out = (ctypes.c_uint8 * len(data))()
        _lib.k3_mascara(buf, out, ctypes.c_uint32(len(data)), ctypes.c_uint16(semilla))
        return bytes(out)