"""MechanicalProvider v0.1 — vector de estado + incremento log (R001)."""

from __future__ import annotations

import math
from dataclasses import dataclass

from mdc_lib.regla_calculo import ReglaCalculo


@dataclass
class MechanicalState:
    position: float
    log_state: float


class MechanicalProvider:
    """Actualiza estado en lugar de reconstruir — regla mecánica toy."""

    def __init__(self) -> None:
        self._regla = ReglaCalculo()
        self._state = MechanicalState(position=0.0, log_state=0.0)

    @property
    def state(self) -> MechanicalState:
        return self._state

    def step(self, delta: float) -> float:
        self._state.position += delta
        if self._state.position > 0:
            self._state.log_state = math.log10(self._state.position)
        return self._state.log_state

    def multiply_increment(self, a: float, b: float) -> float:
        """log(ab) ≈ log(a) + log(b) vía regla de cálculo."""
        return self._regla.multiplicar_valores(a, b)