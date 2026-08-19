#!/usr/bin/env python3
"""Demo: organismo + criba híbrida — compara sintonía vs criba clásica."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aleatorovix import AleatorovixOrganismo
from nucleo.criba import criba_desmemoriada, criba_hibrida

LIMITE = 200
CICLOS = 30


def main() -> None:
    print(f"[*] Aleatorovix: {CICLOS} ciclos buscando candidatos 6k±1")
    org = AleatorovixOrganismo(usar_ping=False)
    org.respirar(n=CICLOS, verbose=False)
    sintonizados = sorted(set(org.primos_encontrados))
    print(f"    Primos por organismo: {len(sintonizados)} -> {sintonizados[:15]}{'...' if len(sintonizados) > 15 else ''}")

    clasicos = criba_desmemoriada(LIMITE)
    hibridos = criba_hibrida(LIMITE)
    print(f"[*] Criba desmemoriada hasta {LIMITE}: {len(clasicos)} primos")
    print(f"[*] Criba híbrida hasta {LIMITE}: {len(hibridos)} primos")
    print(f"[*] Coinciden: {clasicos == hibridos}")


if __name__ == "__main__":
    main()