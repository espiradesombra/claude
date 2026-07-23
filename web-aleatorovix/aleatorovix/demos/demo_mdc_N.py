#!/usr/bin/env python3
"""Demo: resonancia MDC sobre N — perfil d(m) con memoria adaptativa."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aleatorovix import AleatorovixOrganismo
from nucleo.mdc_memoria import AleatorovixMemory, diente_frac

# 10403 = 101 * 103 (ejemplo del libro)
N = 10_403


def main() -> None:
    print(f"[*] MDC Aleatorovix — N = {N}")
    mem = AleatorovixMemory()
    tray = mem.explorar_mdc(N, pasos=12)
    print("    m       d(m)     |d-0.5|")
    for m, d in tray:
        print(f"    {m:5d}   {d:.6f}   {abs(d - 0.5):.6f}")

    print(f"\n[*] Organismo con --n-mdc {N} (5 ciclos resonancia)")
    org = AleatorovixOrganismo(n_mdc=N, usar_ping=False)
    org.respirar(n=5)

    # Factor exacto si d(m) ~ 0.5
    mejor = min(tray, key=lambda t: abs(t[1] - 0.5))
    m_best, d_best = mejor
    denom = 2 * (2 * m_best + 3)
    if abs(d_best - 0.5) < 1e-9:
        factor = N // (N // denom) if N % denom == 0 else None
        print(f"\n[*] Diente óptimo m={m_best} d={d_best:.6f}")
    else:
        print(f"\n[*] Diente más cercano a 0.5: m={m_best} d={d_best:.6f}")
        print(f"    diente_frac verificación: {diente_frac(N, m_best):.6f}")


if __name__ == "__main__":
    main()