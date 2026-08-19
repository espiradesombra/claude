"""Acciones del Select Case mutante — criba, saltos 6k±1, resonancia."""
from __future__ import annotations

import time
from dataclasses import dataclass, field

from nucleo.criba import es_primo
from nucleo.mdc_memoria import AleatorovixMemory, diente_frac
from nucleo.mascara_lila import DecisionLila

ACCION_OLVIDO = 0
ACCION_SALTA_6K_POS = 1
ACCION_SALTA_6K_NEG = 2
ACCION_RESONANCIA = 3

ACCION_NOMBRES = {
    ACCION_OLVIDO: "Criba Desmemoriada (Olvido)",
    ACCION_SALTA_6K_POS: "Salto 6k+1 (Salta)",
    ACCION_SALTA_6K_NEG: "Salto 6k-1 (Baila)",
    ACCION_RESONANCIA: "Inercia ZypyZape (Resonancia)",
}


@dataclass
class ResultadoAccion:
    accion: int
    k_semilla: int
    candidato: int | None = None
    es_primo: bool = False
    mensaje: str = ""
    m_mdc: int | None = None
    d_mdc: float | None = None


@dataclass
class EstadoVolatil:
    """Buffer que se borra tras cada ciclo — criba desmemoriada en RAM."""
    codigo_operar: bytearray = field(
        default_factory=lambda: bytearray(
            b"Lila v1.0 - 10^{1/x} + gemelos + pila + mierda + olvido"
        )
    )
    ultimo_resultado: ResultadoAccion | None = None
    traza: list[str] = field(default_factory=list)

    def desmemoriar(self) -> None:
        self.codigo_operar[:] = b"\x00" * len(self.codigo_operar)
        self.ultimo_resultado = None
        self.traza.clear()


class AccionesLila:
    def __init__(self, n_mdc: int | None = None, direccion: str = "rtl") -> None:
        self.n_mdc = n_mdc
        self.direccion = direccion  # rtl = pinza desde arriba, derecha→izquierda
        self.memoria = AleatorovixMemory()
        self._m_actual = 10_000 if direccion == "rtl" else 1

    @staticmethod
    def salto_6k_pos(k: int) -> int:
        return 6 * k + 1

    @staticmethod
    def salto_6k_neg(k: int) -> int:
        v = 6 * k - 1
        return v if v > 1 else 5

    def ejecutar(self, decision: DecisionLila, k_semilla: int) -> ResultadoAccion:
        accion = decision.decision
        if accion == ACCION_OLVIDO:
            return ResultadoAccion(
                accion=accion,
                k_semilla=k_semilla,
                mensaje="Rastro borrado; sector reiniciado sin memoria de ciclo anterior.",
            )
        if accion == ACCION_SALTA_6K_POS:
            cand = self.salto_6k_pos(k_semilla)
            primo = es_primo(cand)
            return ResultadoAccion(
                accion=accion,
                k_semilla=k_semilla,
                candidato=cand,
                es_primo=primo,
                mensaje=f"Candidato 6k+1 = {cand}" + (" [PRIMO]" if primo else ""),
            )
        if accion == ACCION_SALTA_6K_NEG:
            cand = self.salto_6k_neg(k_semilla)
            primo = es_primo(cand)
            return ResultadoAccion(
                accion=accion,
                k_semilla=k_semilla,
                candidato=cand,
                es_primo=primo,
                mensaje=f"Candidato 6k-1 = {cand}" + (" [PRIMO]" if primo else ""),
            )
        # Resonancia ZypyZape
        salto = self.memoria.resonancia(self._m_actual, decision.medida)
        if self.direccion == "rtl":
            self._m_actual = max(1, self._m_actual - salto)
        else:
            self._m_actual += salto
        d_val = None
        if self.n_mdc is not None and self.n_mdc > 1:
            d_val = diente_frac(self.n_mdc, self._m_actual)
        time.sleep((decision.medida % 5) / 5000.0)
        return ResultadoAccion(
            accion=accion,
            k_semilla=k_semilla,
            m_mdc=self._m_actual,
            d_mdc=d_val,
            mensaje=f"Resonancia: m={self._m_actual} salto={salto}"
            + (f" d(m)={d_val:.6f}" if d_val is not None else ""),
        )