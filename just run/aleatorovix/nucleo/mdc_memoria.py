"""MDC v6 — memoria adaptativa y perfil diente de sierra."""
from __future__ import annotations

import math
from fractions import Fraction


def diente_frac(n: int, m: int) -> float:
    """frac(N / (2*(2m+3))) con aritmética exacta."""
    denom = 2 * (2 * m + 3)
    return float(Fraction(n % denom, denom))


def es_candidato_6k(v: int) -> bool:
    return v % 6 in (1, 5)


def siguiente_m_6k(m: int, positivo: bool) -> int:
    """Avanza en L1 = {6k±1} por índice m."""
    if positivo:
        return m + 1
    return max(1, m - 1)


class AleatorovixMemory:
    """Sectores calientes/fríos — evita repetir zonas ya visitadas."""

    def __init__(self) -> None:
        self.heat_map: dict[int, float] = {}
        self.k_inercia: float = 2.0

    def registrar_sector(self, m: int, desfase: float, acel: float) -> None:
        self.heat_map[m] = desfase + abs(acel)

    def es_zona_fria(self, m: int, umbral: float = 0.4) -> bool:
        for m_prev, score in self.heat_map.items():
            if abs(m - m_prev) < 1000 and score > umbral:
                return True
        return False

    def resonancia(self, m: int, medida: int) -> int:
        """Inercia ZypyZape: acelera o frena según calor del sector."""
        desfase = (m % 97) / 97.0
        acel = ((m * 13) % 41) / 41.0 - 0.5
        self.registrar_sector(m, desfase, acel)
        if desfase < 0.1:
            self.k_inercia = 2.0
        else:
            self.k_inercia += 0.5 + (medida / 20.0)
        salto = max(1, int(self.k_inercia * (1 + medida / 10.0)))
        if self.es_zona_fria(m):
            salto *= 2
        return salto

    def explorar_mdc(self, n: int, m_inicio: int | None = None, pasos: int = 8) -> list[tuple[int, float]]:
        """Recorre dientes midiendo d(m) para factorización orientada."""
        if m_inicio is None:
            m_inicio = max(1, math.isqrt(n) // 4)
        trayectoria: list[tuple[int, float]] = []
        m = m_inicio
        for _ in range(pasos):
            d = diente_frac(n, m)
            desfase = abs(d - 0.5)
            acel = d - (trayectoria[-1][1] if trayectoria else d)
            self.registrar_sector(m, desfase, acel)
            trayectoria.append((m, d))
            if self.es_zona_fria(m):
                m += 3
            elif desfase < 0.05:
                m += 1
            else:
                m += max(1, int(desfase * 20))
        return trayectoria