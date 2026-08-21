"""K3 clásico 32-bit (fingerprint; no reversible solo)."""

from __future__ import annotations

OFFSETS = (5, 7, 9, 11, 13, 17, 19, 23)
MAGICO = 0x9E3779B1
FINAL_MUL_A = 0x85EBCA6B
FINAL_MUL_B = 0xC2B2AE35
SEED = 0x1F2E3D4C
N_REGS = 4


def _rotl32(v: int, n: int) -> int:
    n &= 31
    v &= 0xFFFFFFFF
    if n == 0:
        return v
    return ((v << n) | (v >> (32 - n))) & 0xFFFFFFFF


def _compress(estado: list[int], b: int) -> None:
    n = N_REGS
    b &= 0xFFFFFFFF
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
    x ^= x >> 15
    prev0 = estado[0]
    estado[0] = x & 0xFFFFFFFF
    temp_prev = prev0
    for i in range(1, n):
        temp_actual = estado[i]
        estado[i] = (_rotl32(estado[i], OFFSETS[0] + i) ^ temp_prev) & 0xFFFFFFFF
        temp_prev = temp_actual
        if i == n - 1:
            estado[i] = (estado[i] ^ b) & 0xFFFFFFFF


def _finalize(regs: list[int]) -> int:
    acc = regs[0]
    for i in range(1, N_REGS):
        acc ^= _rotl32(regs[i], 5 + i)
        acc &= 0xFFFFFFFF
    acc ^= acc >> 16
    acc = (acc * FINAL_MUL_A) & 0xFFFFFFFF
    acc ^= acc >> 13
    acc = (acc * FINAL_MUL_B) & 0xFFFFFFFF
    acc ^= acc >> 16
    return acc & 0xFFFFFFFF


def k3_bytes(data: bytes, seed: int = SEED) -> int:
    regs = [(seed ^ (i * MAGICO)) & 0xFFFFFFFF for i in range(N_REGS)]
    if len(data) % 4:
        data = data + b"\x00" * (4 - len(data) % 4)
    for i in range(0, len(data), 4):
        bloque = (data[i] << 24) | (data[i + 1] << 16) | (data[i + 2] << 8) | data[i + 3]
        _compress(regs, bloque)
    return _finalize(regs)


def hash_k3_plain(data: bytes) -> dict:
    d = k3_bytes(data)
    return {
        "family": "k3",
        "reversible": False,
        "digest": d,
        "digest_hex": f"{d:08x}",
        "meta": {"note": "fingerprint only"},
    }
