"""
K3 hash — port Python del k3hash.c oficial en HASHTOOLCODE zip.

Fuente autoritativa:
  HASHTOOLCODE(l)L()(L ) (L).zip → k3hash/k3hash/src/k3hash.c
  (copia en native/k3hash/src/k3hash.c)

Usar hash_engine.py como entrada principal (DLL nativa o este fallback).
"""

from __future__ import annotations

from dataclasses import dataclass

OFFSETS = (5, 7, 9, 11, 13, 17, 19, 23)
MAGICO = 0x9E3779B1
FINAL_MUL_A = 0x85EBCA6B
FINAL_MUL_B = 0xC2B2AE35
MAX_REGISTROS = 10


def _rotl32(value: int, positions: int) -> int:
    positions &= 31
    if positions == 0:
        return value & 0xFFFFFFFF
    value &= 0xFFFFFFFF
    return ((value << positions) | (value >> (32 - positions))) & 0xFFFFFFFF


def _ejecutar_compresion_k3(estado: list[int], valor_bloque: int, n: int) -> None:
    b = valor_bloque & 0xFFFFFFFF
    x = (estado[0] ^ (b & estado[1 % n])) & 0xFFFFFFFF

    for j in range(2, n):
        x ^= _rotl32(estado[j], OFFSETS[j % 8])
    x ^= _rotl32(estado[1 % n], OFFSETS[4])
    x = (x + _rotl32(estado[0], OFFSETS[5])) & 0xFFFFFFFF
    x ^= (b * MAGICO) & 0xFFFFFFFF
    x ^= b
    x = (x + _rotl32(b, OFFSETS[6])) & 0xFFFFFFFF
    x ^= _rotl32(b, OFFSETS[7])
    x ^= _rotl32(x, OFFSETS[2])
    x = (x * MAGICO) & 0xFFFFFFFF
    x ^= (x >> 15)

    estado_anterior_0 = estado[0]
    estado[0] = x & 0xFFFFFFFF

    temp_prev = estado_anterior_0
    for i in range(1, n):
        temp_actual = estado[i]
        estado[i] = (_rotl32(estado[i], OFFSETS[0] + i) ^ temp_prev) & 0xFFFFFFFF
        temp_prev = temp_actual
        if i == n - 1:
            estado[i] ^= b
            estado[i] &= 0xFFFFFFFF


def _extraer_hash_final(registros: list[int], n: int) -> int:
    acumulador = registros[0]
    for i in range(1, n):
        acumulador ^= _rotl32(registros[i], 5 + i)
        acumulador &= 0xFFFFFFFF
    acumulador ^= acumulador >> 16
    acumulador = (acumulador * FINAL_MUL_A) & 0xFFFFFFFF
    acumulador ^= acumulador >> 13
    acumulador = (acumulador * FINAL_MUL_B) & 0xFFFFFFFF
    acumulador ^= acumulador >> 16
    return acumulador & 0xFFFFFFFF


@dataclass
class K3HashConfig:
    bits_ancho: int = 32
    num_registros: int = 4
    semilla_inicial: int = 0x1F2E3D4C


def k3_config_default() -> K3HashConfig:
    return K3HashConfig()


def _normalize_config(config: K3HashConfig) -> int:
    if config.bits_ancho not in (8, 16, 32):
        config.bits_ancho = 32
    config.num_registros = max(2, min(config.num_registros, MAX_REGISTROS))
    return config.bits_ancho // 8


class K3HashCtx:
    def __init__(self, config: K3HashConfig | None = None) -> None:
        self.config = config or k3_config_default()
        self.ancho_bytes = _normalize_config(self.config)
        self.registros = [
            (self.config.semilla_inicial ^ (i * MAGICO)) & 0xFFFFFFFF
            for i in range(self.config.num_registros)
        ]
        self.buffer_parcial = bytearray(4)
        self.bytes_en_buffer = 0

    def update(self, data: bytes) -> None:
        n = self.config.num_registros
        ancho = self.ancho_bytes
        i = 0
        while i < len(data):
            while self.bytes_en_buffer < ancho and i < len(data):
                self.buffer_parcial[self.bytes_en_buffer] = data[i]
                self.bytes_en_buffer += 1
                i += 1
            if self.bytes_en_buffer == ancho:
                bloque = 0
                for b in range(ancho):
                    bloque = (bloque << 8) | self.buffer_parcial[b]
                _ejecutar_compresion_k3(self.registros, bloque, n)
                self.bytes_en_buffer = 0

    def final(self) -> int:
        if self.bytes_en_buffer > 0:
            bloque = 0
            for b in range(self.ancho_bytes):
                byte = self.buffer_parcial[b] if b < self.bytes_en_buffer else 0
                bloque = (bloque << 8) | byte
            _ejecutar_compresion_k3(self.registros, bloque, self.config.num_registros)
            self.bytes_en_buffer = 0
        return _extraer_hash_final(self.registros, self.config.num_registros)


def k3_hash_buffer(data: bytes, config: K3HashConfig | None = None) -> int:
    ctx = K3HashCtx(config)
    ctx.update(data)
    return ctx.final()


def k3_hash(data: bytes, seed: int | None = None) -> int:
    """Convenience wrapper used by AntiPC pipelines."""
    cfg = k3_config_default()
    if seed is not None:
        cfg.semilla_inicial = seed & 0xFFFFFFFF
    return k3_hash_buffer(data, cfg)