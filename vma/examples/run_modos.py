#!/usr/bin/env python3
"""Demostración de los 3 modos operativos K3."""

from vma.k3 import K3CrossVerifier

SEMILLA = 0x1F2E3D4C


def modo_usual() -> None:
    volcado = b"\x00\x12\x7F\x8A" * 16
    v = K3CrossVerifier(bits_ancho=32, desfase_stride=0, num_registros=3, semilla=SEMILLA)
    print(f"[MODO USUAL]     {v.verificar_memoria(volcado)}")


def modo_propio_a() -> None:
    volcado = b"\xAA\xBB\x12\x34\xCC\xDD\x56\x78" * 4
    v = K3CrossVerifier(bits_ancho=16, desfase_stride=2, num_registros=3, semilla=SEMILLA)
    print(f"[MODO PROPIO A]  {v.verificar_memoria(volcado)}")


def modo_propio_b() -> None:
    datos = b"MUESTRA_ESTABLE_DE_VOLTAJE_NODO_4"
    v = K3CrossVerifier(bits_ancho=32, desfase_stride=0, num_registros=4, semilla=SEMILLA)
    v.encolar_marca_final(0x00000004, 0x10)
    v.encolar_marca_final(0x99999999, 0xEE)
    print(f"[MODO PROPIO B]  {v.verificar_memoria(datos)}")


if __name__ == "__main__":
    modo_usual()
    modo_propio_a()
    modo_propio_b()