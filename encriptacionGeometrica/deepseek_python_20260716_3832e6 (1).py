from k3_bindings import K3Engine

clave = {
    "tales": [3,5,8,13,21],
    "figuras": ["equilatero", "isosceles", "escaleno"],
    "puntos": [6,12,18],
    "saltos": [1,0,2],
    "primos": [(3,5),(7,11)],
    "porcentajes_aportacion": [50,100],
    "iteraciones_pi": 15
}

motor = K3Engine(clave)
hash_fase = motor.encriptar(43210)
semilla = motor.desencriptar(hash_fase)