#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo de uso completo del motor K3 33x1.
Muestra:
- Cómo configurar la clave geométrica.
- Cómo colapsar una semilla en fase.
- Cómo recuperar la semilla.
- Cómo aplicar la máscara lila a un mensaje.
"""

from k3_bindings import K3Engine

# 1. Definir la clave geométrica (parámetros 33x1)
clave = {
    "tales": [3, 5, 8, 13, 21],
    "figuras": ["equilatero", "isosceles", "escaleno"],
    "puntos": [6, 12, 18],
    "saltos": [1, 0, 2],
    "primos": [(3, 5), (7, 11)],
    "porcentajes_aportacion": [50, 100],
    "iteraciones_pi": 15
}

motor = K3Engine(clave)

# 2. Semilla original (16 bits)
semilla_original = 43210
print(f"Semilla original: {semilla_original} (0x{semilla_original:04X})")

# 3. Colapsar a fase
hash_fase = motor.encriptar(semilla_original)
print(f"Hash de fase: {hash_fase:.18Lf}")

# 4. Recuperar semilla
semilla_recuperada = motor.desencriptar(hash_fase)
print(f"Semilla recuperada: {semilla_recuperada}")

if semilla_recuperada == semilla_original:
    print("✅ ÉXITO: Reversibilidad perfecta.")
else:
    print("❌ ERROR: No se recuperó la semilla.")

# 5. Aplicar máscara lila a un mensaje
mensaje = b"Proyecto 33x1 - Movilidad Sostenible"
print(f"\nMensaje original: {mensaje}")

# Cifrar
cifrado = motor.aplicar_mascara(mensaje, semilla_original)
print(f"Mensaje cifrado (hex): {cifrado.hex()}")

# Descifrar (usando la misma semilla)
descifrado = motor.aplicar_mascara(cifrado, semilla_original)
print(f"Mensaje descifrado: {descifrado}")

if mensaje == descifrado:
    print("✅ Máscara lila funciona correctamente.")
else:
    print("❌ Error en la máscara.")