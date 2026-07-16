"""
K3 Cross Verifier — auditoría industrial de memoria (Theorem K3 v0).

Motor en anillo con firma activa anti-tamper, stride configurable
y cola de marcas finales con banderitas.
"""

from __future__ import annotations

K3_FIRMA_LOGICA = b"x^=(B*0x9E3779B1);x+=rotar(B,6);estado[0]=x;"
MASK_32 = 0xFFFFFFFF
MAGICO = 0x9E3779B1
OFFSETS = (5, 7, 9, 11, 13, 17, 19, 23)
MAX_MARCAS = 4


def _rotl32(valor: int, posiciones: int) -> int:
    posiciones &= 31
    if posiciones == 0:
        return valor & MASK_32
    return ((valor << posiciones) | (valor >> (32 - posiciones))) & MASK_32


def _bytes_to_block(fragmento: bytes, bytes_ancho: int) -> int:
    valor = 0
    for b in fragmento:
        valor = (valor << 8) | b
    if len(fragmento) < bytes_ancho:
        valor <<= 8 * (bytes_ancho - len(fragmento))
    return valor & MASK_32


def _ejecutar_compresion(estado: list[int], valor_bloque: int, num_registros: int) -> list[int]:
    n = num_registros
    b = valor_bloque & MASK_32

    x = estado[0] ^ (b & estado[1 % n])
    for j in range(2, n):
        x ^= _rotl32(estado[j], OFFSETS[j % 8])
    x ^= _rotl32(estado[1 % n], OFFSETS[4])
    x = (x + _rotl32(estado[0], OFFSETS[5])) & MASK_32
    x ^= (b * MAGICO) & MASK_32

    x ^= b
    x = (x + _rotl32(b, OFFSETS[6])) & MASK_32
    x ^= _rotl32(b, OFFSETS[7])
    x ^= _rotl32(x, OFFSETS[2])
    x = (x * MAGICO) & MASK_32
    x ^= x >> 15

    nuevo = [0] * n
    prev = estado[0]
    nuevo[0] = x
    temp_prev = prev
    for i in range(1, n):
        temp_actual = estado[i]
        nuevo[i] = (_rotl32(estado[i], OFFSETS[0] + i) ^ temp_prev) & MASK_32
        temp_prev = temp_actual
        if i == n - 1:
            nuevo[i] ^= b
    return nuevo


def _squeeze(estado: list[int], num_registros: int) -> int:
    acumulador = estado[0]
    for i in range(1, num_registros):
        acumulador ^= _rotl32(estado[i], 5 + i)
    acumulador ^= acumulador >> 16
    acumulador = (acumulador * 0x85EBCA6B) & MASK_32
    acumulador ^= acumulador >> 13
    acumulador = (acumulador * 0xC2B2AE35) & MASK_32
    acumulador ^= acumulador >> 16
    return acumulador & MASK_32


class K3CrossVerifier:
    """Verificador cruzado universal de volcados RAM/Flash."""

    def __init__(
        self,
        bits_ancho: int = 32,
        desfase_stride: int = 0,
        num_registros: int = 3,
        semilla: int = 0x1F2E3D4C,
    ) -> None:
        if bits_ancho not in (8, 16, 32):
            bits_ancho = 32
        self.bits_ancho = bits_ancho
        self.bytes_ancho = bits_ancho // 8
        self.desfase_stride = max(0, desfase_stride)
        self.num_registros = max(2, num_registros)
        self.semilla_inicial = semilla & MASK_32
        self.firma_logica = K3_FIRMA_LOGICA
        self.cola_marcas: list[dict[str, int]] = []

    def encolar_marca_final(self, valor_bits: int, banderita: int) -> None:
        if len(self.cola_marcas) >= MAX_MARCAS:
            return
        self.cola_marcas.append({
            "valor": valor_bits & MASK_32,
            "flag": banderita & 0xFF,
        })

    def verificar_memoria(self, mapa_memoria: bytes) -> str:
        estado = [
            (self.semilla_inicial ^ (i * MAGICO)) & MASK_32
            for i in range(self.num_registros)
        ]

        idx_f = 0
        while idx_f < len(self.firma_logica):
            fragmento = self.firma_logica[idx_f : idx_f + self.bytes_ancho]
            if not fragmento:
                break
            bloque = _bytes_to_block(fragmento, self.bytes_ancho)
            estado = _ejecutar_compresion(estado, bloque, self.num_registros)
            idx_f += self.bytes_ancho

        indice = 0
        longitud = len(mapa_memoria)
        while indice < longitud:
            fragmento = mapa_memoria[indice : indice + self.bytes_ancho]
            if not fragmento:
                break
            bloque = _bytes_to_block(fragmento, self.bytes_ancho)
            estado = _ejecutar_compresion(estado, bloque, self.num_registros)
            indice += self.bytes_ancho + self.desfase_stride

        for marca in self.cola_marcas:
            compuesto = (marca["valor"] ^ (marca["flag"] << 24)) & MASK_32
            estado = _ejecutar_compresion(estado, compuesto, self.num_registros)

        return f"0x{_squeeze(estado, self.num_registros):08X}"

    def verificar_memoria_int(self, mapa_memoria: bytes) -> int:
        return int(self.verificar_memoria(mapa_memoria), 16)