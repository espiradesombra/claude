"""Máscara de curvatura sagrada, campanas Gauss e intérprete mutante."""
from __future__ import annotations

from dataclasses import dataclass

from nucleo.entropia import MuestraEntropia


@dataclass
class DecisionLila:
    decision: int
    medida: int
    r0: int
    r5: int
    r9: int
    curva: float
    inercia_a: int
    red_x: int


class MascaraLila:
    """(-10^(1/x) + 10^(1/(x+a)) - 1) * medida — organismo Lila v9.9."""

    @staticmethod
    def mascara_lila(x: int, a: int) -> float:
        if x <= 0:
            x = 1
        t1 = -pow(10.0, 1.0 / float(x))
        t2 = pow(10.0, 1.0 / float(x + a))
        return (t1 + t2 - 1.0) * float(x)

    @staticmethod
    def pulso_gauss(valor: int, pico: int) -> int:
        return 1 if abs(valor - pico) <= 1 else 0

    @staticmethod
    def interprete_mutante(bit_externo: int, bit_robado: int, ruido_inercia: int) -> int:
        par = (bit_externo << 1) | (bit_robado & 1)
        rotacion = ruido_inercia % 3
        if par in (0b11, 0b00):
            base = 0
        elif par == 0b01:
            base = 1
        else:
            base = 2
        return (base + rotacion) % 3

    def decidir(self, muestra: MuestraEntropia) -> DecisionLila:
        curva = self.mascara_lila(muestra.red_x, muestra.inercia_a)
        medida = abs(int(curva)) % 10

        b0 = self.pulso_gauss(medida, 0)
        r0 = self.interprete_mutante(b0, muestra.bit_pila, muestra.basura)
        b5 = self.pulso_gauss(medida, 5)
        r5 = self.interprete_mutante(b5, muestra.bit_pila, muestra.basura ^ 1)
        b9 = self.pulso_gauss(medida, 9)
        r9 = self.interprete_mutante(b9, muestra.bit_pila, muestra.basura ^ muestra.nanos)

        decision = (r0 + r5 + r9) % 4
        return DecisionLila(
            decision=decision,
            medida=medida,
            r0=r0,
            r5=r5,
            r9=r9,
            curva=curva,
            inercia_a=muestra.inercia_a,
            red_x=muestra.red_x,
        )