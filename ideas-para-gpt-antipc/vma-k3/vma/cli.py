"""CLI vma-k3 — verificación K3 desde terminal."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from vma.k3 import K3CrossVerifier


def _parse_mark(s: str) -> tuple[int, int]:
    parts = s.split(":")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("formato marca: VALOR:FLAG (hex o decimal)")
    return int(parts[0], 0), int(parts[1], 0)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="vma-k3",
        description="Auditoría K3 industrial — Theorem K3 v0 (VMA)",
    )
    p.add_argument("input", nargs="?", help="Texto o ruta a archivo binario")
    p.add_argument("--file", "-f", help="Leer datos desde archivo")
    p.add_argument("--bits", type=int, default=32, choices=[8, 16, 32])
    p.add_argument("--desfase", type=int, default=0)
    p.add_argument("--registros", type=int, default=3)
    p.add_argument("--seed", type=lambda x: int(x, 0), default=0x1F2E3D4C)
    p.add_argument("--mark", action="append", type=_parse_mark, default=[], metavar="VAL:FLAG")
    p.add_argument("--modo", choices=["usual", "propio-a", "propio-b"], help="Preset operativo")
    args = p.parse_args(argv)

    if args.modo == "usual":
        args.bits, args.desfase, args.registros = 32, 0, 3
    elif args.modo == "propio-a":
        args.bits, args.desfase, args.registros = 16, 2, 3
    elif args.modo == "propio-b":
        args.bits, args.desfase, args.registros = 32, 0, 4

    if args.file:
        data = Path(args.file).read_bytes()
    elif args.input:
        pth = Path(args.input)
        data = pth.read_bytes() if pth.is_file() else args.input.encode("utf-8")
    else:
        data = sys.stdin.buffer.read()

    v = K3CrossVerifier(
        bits_ancho=args.bits,
        desfase_stride=args.desfase,
        num_registros=args.registros,
        semilla=args.seed,
    )
    for valor, flag in args.mark:
        v.encolar_marca_final(valor, flag)

    print(v.verificar_memoria(data))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())