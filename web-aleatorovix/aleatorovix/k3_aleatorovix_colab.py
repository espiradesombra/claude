#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K3 + Aleatorovix — demo headless (Colab / sin pantalla) ❤️🌹

No usa tkinter. pandas es opcional (tablas bonitas si está instalado).

En Colab:
  !python k3_aleatorovix_colab.py
  # o copiar celdas y llamar: main()

En local sin GUI:
  python k3_aleatorovix_colab.py
  python k3_aleatorovix_colab.py ventanas
  python k3_aleatorovix_colab.py ciclo
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

# Import solo del motor (el módulo GUI también define el motor; en Colab
# tkinter a veces no está — fallback embebido si falla el import).
try:
    from k3_aleatorovix_gui import (
        PAIR,
        K_HEART,
        L_ROSE,
        AleatorovixOrganismo,
        bits_a_flores,
        entero_a_flores,
        numero_floral,
        formatear_valor_floral,
        digito_a_flores,
    )
    _IMPORT_OK = True
except Exception as _e:  # noqa: BLE001
    _IMPORT_OK = False
    _IMPORT_ERR = _e


def _try_dataframe(rows: list) -> None:
    """Muestra tabla con pandas si hay; si no, print simple."""
    try:
        import pandas as pd

        df = pd.DataFrame(rows)
        # columnas útiles primero
        cols = [
            c
            for c in (
                "ventana",
                "start_idx",
                "end_idx",
                "ancho",
                "flores",
                "floral_digitos",
                "valor",
            )
            if c in df.columns
        ]
        df = df[cols] if cols else df
        try:
            from IPython.display import display  # type: ignore

            display(df)
        except Exception:
            print(df.to_string(index=False))
    except ImportError:
        for r in rows:
            print(
                f"  v{r.get('ventana')}: "
                f"[{r.get('start_idx')}:{r.get('end_idx')}] "
                f"ancho={r.get('ancho')}  "
                f"{r.get('flores')}  |  {r.get('floral_digitos')}"
            )


def demo_ciclo() -> None:
    print(f"\n=== Ciclo Aleatorovix {PAIR} (sin GUI) ===\n")
    org = AleatorovixOrganismo(usar_ping=False)
    r = org.ciclo_completo()
    print(f"Ciclo {numero_floral(org.ciclos_completados)} {PAIR}")
    print(formatear_valor_floral("inercia(a)", r["inercia_a"]))
    print(formatear_valor_floral("red(x)", r["red_x"]))
    print(f"  decisión: {r['decision_floral']}  {r['mensaje']}")
    if r["candidato_floral"]:
        tag = f" {L_ROSE}PRIMO{L_ROSE}" if r["primo"] else ""
        print(f"  candidato: {r['candidato_floral']}{tag}")
    print(f"  chorro:   {r['flores_acordeon']}")
    print(f"  posible:  {r['flores_posible']}")
    print("\n" + r["acordeon_vista"])
    print("-" * 50)


def demo_ventanas(
    total_bits: int = 128,
    num_ventanas: int = 8,
    min_ancho: int = 8,
    max_ancho: int = 32,
) -> None:
    print(f"\n=== Ventanas móviles (acordeón) {PAIR} ===\n")
    print(
        f"Chorro={total_bits} bits | ventanas={num_ventanas} | "
        f"ancho∈[{min_ancho},{max_ancho}]\n"
        f"Leyenda: {K_HEART}=1  {L_ROSE}=0  (números no arábigos)\n"
    )
    org = AleatorovixOrganismo(usar_ping=False)
    data = org.acordeon.chorro_ventanas_moviles(
        total_bits=total_bits,
        num_ventanas=num_ventanas,
        min_ancho=min_ancho,
        max_ancho=max_ancho,
    )
    _try_dataframe(data)
    print()
    for w in data:
        print(f"── ventana {entero_a_flores(w['ventana'], 4)} ──")
        print(f"  flores: {w['flores']}")
        print(f"  dígitos-flor: {w['floral_digitos']}")
        for line in w["vista"].splitlines():
            print(f"    {line}")
        print()
    print(f"=== fin ventanas {PAIR} ===\n")


def demo_alfabeto() -> None:
    print(f"Alfabeto nibble 0–9 {PAIR}:")
    for d in range(10):
        print(f"  {d} → {digito_a_flores(d)}")


def main() -> None:
    if not _IMPORT_OK:
        print("ERROR: no se pudo importar k3_aleatorovix_gui.py")
        print("Asegúrate de subir ambos archivos al mismo directorio.")
        print(f"Detalle: {_IMPORT_ERR}")
        sys.exit(1)

    arg = (sys.argv[1] if len(sys.argv) > 1 else "todo").lower()
    if arg in ("ciclo", "sim", "demo"):
        demo_ciclo()
    elif arg in ("ventanas", "windows", "movil", "móvil"):
        demo_ventanas()
    elif arg in ("alfa", "alfabeto"):
        demo_alfabeto()
    else:
        demo_ciclo()
        demo_ventanas()
        demo_alfabeto()
        print(
            "\nTip Colab: sube k3_aleatorovix_gui.py + este archivo juntos.\n"
            "GUI local: python k3_aleatorovix_gui.py\n"
        )


if __name__ == "__main__":
    main()
