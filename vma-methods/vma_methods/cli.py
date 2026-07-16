#!/usr/bin/env python3
"""CLI vma-methods — cribas, criva, newton."""

from __future__ import annotations

import argparse
import io
import sys

from vma_methods import cribas, criva, newton, classic


def _cmd_criba(args: argparse.Namespace) -> int:
    limit = int(args.limit)
    mapping = {
        "desmemoriada": cribas.CribaDesmemoriada,
        "6k": cribas.CribaModular6k,
        "hibrida": cribas.CribaHibrida,
    }
    cls = mapping[args.tipo]
    sieve = cls(limit)
    if args.tipo == "hibrida" and args.segmented:
        primes = sieve.segmented_run(seg_size=args.seg_size)
        label = "Híbrida segmentada"
    else:
        primes = sieve.run(verbose=args.verbose)
        label = args.tipo
    print(f"{label} hasta {limit}: {len(primes)} primos")
    if args.show:
        n = min(args.show, len(primes))
        print(primes[:n], "..." if len(primes) > n else "")
    return 0


def _cmd_compare_cribas(args: argparse.Namespace) -> int:
    buf = io.StringIO()
    old = sys.stdout
    sys.stdout = buf
    try:
        cribas.comparar_cribas(limit=args.limit, verbose=True)
    finally:
        sys.stdout = old
    print(buf.getvalue())
    return 0


def _cmd_criva(args: argparse.Namespace) -> int:
    if args.compare:
        xs = [int(x.strip()) for x in args.compare.split(",")]
        criva.compare_criva_vs_pnt(xs)
        return 0
    x = float(args.x)
    est = criva.criva(x) * x
    exact = classic.pi_count(int(x))
    pnt = classic.pnt_estimate(x)
    print(f"x = {x:.0f}")
    print(f"π(x) exacto     : {exact}")
    print(f"Criva π(x) ~     : {est:.1f}")
    print(f"PNT π(x) ~       : {pnt:.1f}")
    print(f"densidad Criva   : {criva.criva(x):.8f}")
    print(f"densidad exacta  : {classic.prime_density(x):.8f}")
    return 0


def _cmd_newton(args: argparse.Namespace) -> int:
    E = float(args.E)
    b = float(args.b)
    if args.familia == "none":
        r = newton.newton_rapido(E, b=b, verbose=args.verbose)
    else:
        r = newton.log_con_oraculo(E, b=b, familia=args.familia, verbose=args.verbose)
    print(f"E = {E}  base b = {b}")
    print(f"log_b(E) exacto : {r['j_exacto']:.12f}")
    print(f"log_b(E) VMA    : {r['j']:.12f}")
    print(f"iteraciones     : {r['iteraciones']}")
    print(f"error           : {r['error']:.3e}")
    if "j0" in r:
        print(f"j0 oráculo      : {r['j0']:.6f}  ({r.get('familia', '')})")
    return 0


def _cmd_pi(args: argparse.Namespace) -> int:
    x = int(args.x)
    n = classic.pi_count(x)
    print(f"π({x}) = {n}")
    print(f"densidad = {n / x:.8f}" if x else "")
    print(f"PNT ~ {classic.pnt_estimate(x):.1f}")
    print(f"Criva ~ {criva.criva(x) * x:.1f}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="vma-methods",
        description="VMA — cribas, Criva, Newton Rápido",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("criba", help="Ejecutar una criba")
    c.add_argument("--limit", type=int, default=5000)
    c.add_argument(
        "--tipo",
        choices=["desmemoriada", "6k", "hibrida"],
        default="6k",
    )
    c.add_argument("--show", type=int, default=0, help="Mostrar N primeros primos")
    c.add_argument("--segmented", action="store_true", help="Solo híbrida: modo segmentado")
    c.add_argument("--seg-size", type=int, default=500)
    c.add_argument("--verbose", action="store_true")
    c.set_defaults(func=_cmd_criba)

    cc = sub.add_parser("compare-cribas", help="Benchmark 3 cribas VMA")
    cc.add_argument("--limit", type=int, default=5000)
    cc.set_defaults(func=_cmd_compare_cribas)

    cv = sub.add_parser("criva", help="Estimar densidad de primos")
    cv.add_argument("--x", type=float, default=10000)
    cv.add_argument("--compare", help="Lista x separada por comas")
    cv.set_defaults(func=_cmd_criva)

    n = sub.add_parser("newton", help="Newton Rápido log_b(E)")
    n.add_argument("--E", type=float, required=True)
    n.add_argument("--b", type=float, default=10.0)
    n.add_argument(
        "--familia",
        default="general",
        choices=["none", "general", "cuadrados", "cubos", "potencia", "kp", "mersenne"],
    )
    n.add_argument("--verbose", action="store_true")
    n.set_defaults(func=_cmd_newton)

    pi = sub.add_parser("pi", help="π(x) exacto (Eratóstenes)")
    pi.add_argument("--x", type=int, default=1000)
    pi.set_defaults(func=_cmd_pi)

    g = sub.add_parser("gui", help="Abrir interfaz gráfica")
    g.set_defaults(func=lambda _: __import__("vma_methods.gui.app", fromlist=["main"]).main())

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())