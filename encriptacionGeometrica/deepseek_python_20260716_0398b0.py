#!/usr/bin/env python3
from k3 import K3

clave = {
    "tales": [3,5,8,13,21],
    "figuras": ["equilatero","isosceles","escaleno"],
    "puntos": [6,12,18],
    "saltos": [1,0,2],
    "primos": [(3,5),(7,11)],
    "porcentajes_aportacion": [50,100],
    "iteraciones_pi": 15
}

motor = K3(clave)

semilla = 43210
hash_fase = motor.encriptar(semilla)
print(f"Hash: {hash_fase:.18Lf}")
recuperada = motor.desencriptar(hash_fase)
print(f"Recuperada: {recuperada}")

mensaje = b"33x1 funciona"
cifrado = motor.mascara(mensaje, semilla)
descifrado = motor.mascara(cifrado, semilla)
print(f"Original: {mensaje}")
print(f"Descifrado: {descifrado}")