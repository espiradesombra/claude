#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K3 + Aleatorovix completo (CLI + GUI opcional) ❤️🌹
Víctor Manzanares Alberola — uso civil / educativo

Ejecutar:
  python k3_aleatorovix_completo.py              # demo integración
  python k3_aleatorovix_completo.py k3
  python k3_aleatorovix_completo.py aleatorovix
  python k3_aleatorovix_completo.py acordeon     # serie ❤️🌹 (números no estándar)
  python k3_aleatorovix_completo.py ventanas     # ventanas móviles sobre chorro
  python k3_aleatorovix_completo.py industrial   # MUCHOS ❤️🌹 (planta 33×1)
  python k3_aleatorovix_completo.py stats
  python k3_aleatorovix_completo.py gui
  python k3_aleatorovix_completo.py info
  python k3_aleatorovix_completo.py colab        # demos headless (sin GUI)
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from k3_aleatorovix_gui import (  # noqa: E402
    PAIR,
    K_HEART,
    L_ROSE,
    AleatorovixOrganismo,
    AccionesLila,
    bits_a_flores,
    entero_a_flores,
    numero_floral,
    render_acordeon,
)


def demo_k3(n: int = 10) -> None:
    org = AleatorovixOrganismo(usar_ping=False)
    print(f"=== SECUENCIA K3 {PAIR} (floral) ===")
    print("Reglas: 12→+12, 1→+10, 11→+2")
    print(f"{K_HEART}=1  {L_ROSE}=0\n")
    for i, (o, t, b) in enumerate(org.generar_secuencia_k3(n)):
        print(
            f"N={entero_a_flores(i, 4)}: "
            f"{numero_floral(o)} → {numero_floral(t)} → {bits_a_flores(b)}"
        )


def demo_aleatorovix(ciclos: int = 5) -> None:
    org = AleatorovixOrganismo(usar_ping=True)
    print(f"=== ALEATOROVIX {PAIR} ===\n")
    for _ in range(ciclos):
        r = org.ciclo_completo()
        print(f"Ciclo {numero_floral(org.ciclos_completados)} {PAIR}")
        print(f"  a={numero_floral(r['inercia_a'])}  x={numero_floral(r['red_x'])}")
        print(f"  medida={entero_a_flores(r['medida'], 4)}  k={r['k_floral']}")
        print(f"  chorro={r['flores_acordeon']}")
        print(f"  posible={r['flores_posible']}")
        if r["candidato_floral"]:
            print(f"  candidato={r['candidato_floral']}" + (" PRIMO" if r["primo"] else ""))
        print(f"  {r['mensaje']}")
        print(f">>> CRIBA DESMEMORIADA {PAIR}")
        print("-" * 40)


def demo_acordeon(cuantos: int = 12) -> None:
    """
    Modo acordeón: serie de rosas y corazones.
    Cada pliegue muestra un número posible SIN dígitos arábigos estándar.
    """
    org = AleatorovixOrganismo(usar_ping=False)
    print(f"=== MODO ACORDEÓN {PAIR} ===")
    print(f"Leyenda: {K_HEART} = bit 1 (K/msl)   {L_ROSE} = bit 0 (L/lsl)")
    print("Los números se leen solo como flores (binario / nibbles).")
    print("Estirar = espacios entre flores | Plegar = compacto.\n")

    serie = org.serie_acordeon_flores(pliegues=24, cuantos=cuantos)
    for item in serie:
        print(f"── pliegue {entero_a_flores(item['pliegue'], 4)}  "
              f"ancho {entero_a_flores(item['ancho'], 4)} ──")
        print(f"  serie viva:   {item['flores']}")
        print(f"  dígitos-flor: {item['floral_digitos']}")
        print("  pulso:")
        for linea in item["vista"].splitlines():
            print(f"    {linea}")
        print()

    print(f"=== fin acordeón ({entero_a_flores(cuantos)} pliegues) {PAIR} ===")
    from k3_aleatorovix_gui import digito_a_flores

    print(f"\nAlfabeto nibble (dígito → 4 flores):")
    for d in range(10):
        print(f"  [{d}] → {digito_a_flores(d)}")


def demo_ventanas(
    total_bits: int = 128,
    num_ventanas: int = 8,
    min_ancho: int = 8,
    max_ancho: int = 32,
) -> None:
    """Ventanas móviles sobre un chorro largo (idea Colab + DataFrame opcional)."""
    org = AleatorovixOrganismo(usar_ping=False)
    print(f"=== VENTANAS MÓVILES {PAIR} ===")
    print(
        f"chorro={total_bits} bits | ventanas={num_ventanas} | "
        f"ancho∈[{min_ancho},{max_ancho}]\n"
    )
    data = org.acordeon.chorro_ventanas_moviles(
        total_bits=total_bits,
        num_ventanas=num_ventanas,
        min_ancho=min_ancho,
        max_ancho=max_ancho,
    )
    try:
        import pandas as pd

        df = pd.DataFrame(
            [
                {
                    "ventana": w["ventana"],
                    "start": w["start_idx"],
                    "end": w["end_idx"],
                    "ancho": w["ancho"],
                    "flores": w["flores"],
                    "floral": w["floral_digitos"],
                    "valor": w["valor"],
                }
                for w in data
            ]
        )
        print(df.to_string(index=False))
        print()
    except ImportError:
        print("(pandas no instalado — salida texto)\n")

    for w in data:
        print(
            f"ventana {entero_a_flores(w['ventana'], 4)} "
            f"[{w['start_idx']}:{w['end_idx']}] ancho={w['ancho']}"
        )
        print(f"  {w['flores']}")
        print(f"  {w['floral_digitos']}")
        for line in w["vista"].splitlines():
            print(f"    {line}")
        print()
    print(f"=== fin ventanas {PAIR} ===")


def demo_integracion() -> None:
    print("=" * 70)
    print(f"DEMO INTEGRACIÓN: K3 + ALEATOROVIX {PAIR}")
    print("=" * 70)
    org = AleatorovixOrganismo(usar_ping=True)
    print("\n1. Secuencia K3 (floral):")
    print("-" * 70)
    seq = org.generar_secuencia_k3(5)
    for i, (o, t, b) in enumerate(seq):
        print(f"   {entero_a_flores(i, 3)}: {numero_floral(o)} → {numero_floral(t)} → {bits_a_flores(b)}")

    print(f"\n2. Ciclos + acordeón {PAIR}:")
    print("-" * 70)
    for _, t, b in seq:
        r = org.ciclo_completo()
        print(f"\n   K3→ {numero_floral(t)} | k={r['k_floral']}")
        print(f"   chorro {r['flores_acordeon']}")
        print(f"   posible {r['flores_posible']}")
        print(r["acordeon_vista"])
        print(f">>> CRIBA DESMEMORIADA {PAIR}")
        print("----------------------------")


def demo_stats(ciclos: int = 100) -> None:
    from collections import Counter

    org = AleatorovixOrganismo(usar_ping=False)
    for _ in range(ciclos):
        org.ciclo_completo()
    counts = Counter(org.historial_decisiones)
    primos = sorted({c for c in org.primos_encontrados if AccionesLila._es_primo(c)})
    print(f"=== STATS {numero_floral(ciclos)} ciclos {PAIR} ===")
    print(f"Candidatos: {numero_floral(len(org.primos_encontrados))} | "
          f"Primos: {numero_floral(len(primos))} {L_ROSE}")
    for d in range(4):
        print(f"  {AccionesLila.ACCION_NOMBRES[d]}: {numero_floral(counts.get(d, 0))}")
    if primos:
        print("Primos (floral):")
        for p in primos[:15]:
            print(f"  {numero_floral(p)}")
    print(f"{K_HEART} msl  {L_ROSE} lsl")


def demo_info() -> None:
    print(
        f"""
K3 + Aleatorovix {PAIR}
======================
{K_HEART} = K = factor_msl = bit 1 (estirar acordeón)
{L_ROSE} = L = factor_lsl = bit 0 (plegar acordeón)

Números NO estándar:
  - entero → binario → serie de flores
  - dígitos 0-9 → nibble de 4 flores (0000..1001)

Modo acordeón:
  python k3_aleatorovix_completo.py acordeon

Fuentes: organismo_lila_v99.c + gemini-code-1784158392232.py
Uso: civil / educativo (VMA 33×1)
"""
    )


def main() -> None:
    arg = (sys.argv[1] if len(sys.argv) > 1 else "integracion").lower()
    if arg in ("k3",):
        demo_k3()
    elif arg in ("aleatorovix", "lila", "ciclo"):
        demo_aleatorovix()
    elif arg in ("acordeon", "accordion", "flores", "rosas"):
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 12
        demo_acordeon(n)
    elif arg in ("ventanas", "windows", "movil", "móvil", "colab"):
        if arg == "colab":
            demo_aleatorovix(3)
            demo_ventanas()
            demo_acordeon(5)
        else:
            demo_ventanas()
    elif arg in ("industrial", "industria", "planta", "factory"):
        # Delega a la planta floral industrial (miles/millones de ❤️🌹)
        from k3_aleatorovix_industrial import main as industrial_main

        escala = sys.argv[2] if len(sys.argv) > 2 else "industrial"
        industrial_main(["full", escala])
    elif arg in ("stats", "estadisticas"):
        demo_stats()
    elif arg in ("gui", "window"):
        # En Windows no hace falta $DISPLAY; abrir GUI directamente
        from k3_aleatorovix_gui import main as gui_main

        gui_main()
    elif arg in ("info", "help", "-h", "--help"):
        demo_info()
    else:
        demo_integracion()


if __name__ == "__main__":
    main()
