"""
MDC v6 Aleatorovix — memoria adaptativa de sectores (esquema).
Extraído de transcurso.txt — integrar con MDC v5 real cuando esté listo.
"""
from __future__ import annotations

import math


class AleatorovixMemory:
    def __init__(self) -> None:
        self.heat_map: dict[int, float] = {}

    def registrar_sector(self, m: int, desfase: float, acel: float) -> None:
        self.heat_map[m] = desfase + abs(acel)

    def es_zona_fria(self, m: int, umbral: float = 0.4) -> bool:
        for m_prev, score in self.heat_map.items():
            if abs(m - m_prev) < 1000 and score > umbral:
                return True
        return False


def mdc_v6_aleatorovix_demo(n: int, max_pasos: int = 20) -> None:
    memory = AleatorovixMemory()
    k = 2.0
    m = math.isqrt(n) // 2
    print("[*] MDC v6 Aleatorovix (demo memoria adaptativa)")
    pasos = 0
    while m > 1 and pasos < max_pasos:
        if memory.es_zona_fria(m):
            k += 1.0
            m -= max(1, n // int(k * 1_000_000))
            pasos += 1
            continue
        desfase = (m % 97) / 97.0
        acel = ((m * 13) % 41) / 41.0 - 0.5
        memory.registrar_sector(m, desfase, acel)
        if desfase < 0.1:
            k = 2.0
            print(f"  zona caliente m={m} desfase={desfase:.3f}")
        else:
            k += 0.5
            m -= 100
        pasos += 1
    print(f"  sectores visitados: {len(memory.heat_map)}")


if __name__ == "__main__":
    mdc_v6_aleatorovix_demo(10**12 + 39)