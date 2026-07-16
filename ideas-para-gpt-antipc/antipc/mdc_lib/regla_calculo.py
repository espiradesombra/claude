"""Regla de cálculo analógica — port de regla_calculo.hpp (HASHTOOLCODE)."""

from __future__ import annotations

import math


class ReglaCalculo:
    """
    Escala logarítmica + interpolación.
    Uso AntiPC: planificador energético / Battery Mode (no sustituye ALU entera).
    """

    def __init__(self, precision: int = 1000, tamano_escala: int = 5000) -> None:
        self.precision = precision
        self.tamano_escala = tamano_escala
        self.escala_logaritmica = [
            math.pow(10.0, i / precision) for i in range(tamano_escala)
        ]

    def multiplicar(self, pos_a: int, pos_b: int) -> float:
        idx = pos_a + pos_b
        return self.escala_logaritmica[idx] if idx < self.tamano_escala else -1.0

    def posicion_de(self, x: float) -> float:
        if x <= 0.0:
            return -1.0
        return self.precision * math.log10(x)

    def valor_en_posicion(self, pos: float) -> float:
        if pos < 0.0 or pos >= self.tamano_escala - 1:
            return math.pow(10.0, pos / self.precision)
        i0 = int(pos)
        i1 = i0 + 1
        frac = pos - i0
        return self.escala_logaritmica[i0] * (1.0 - frac) + self.escala_logaritmica[i1] * frac

    def multiplicar_valores(self, a: float, b: float) -> float:
        if a <= 0.0 or b <= 0.0:
            return -1.0
        return self.valor_en_posicion(self.posicion_de(a) + self.posicion_de(b))