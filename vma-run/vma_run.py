#!/usr/bin/env python3
"""
VMA RUN — ejecutable unificado (ideas curadas de documentos viejos)
Víctor Manzanares Alberola · 2026

Herramientas en un solo comando:
  k3      auditoría industrial Theorem K3
  mdc     factorización cinemática diofántica
  criva   densidad racional de primos
  criba   criba modular 6k±1
  phase   amplificador de fase K=3 XOR
  hurto   resultados hurto gravitatorio (tabla)
  demo    showcase rápido de todo (<5 s)
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if getattr(sys, "frozen", False):
    # bin\vma-run.exe → just run\
    JUST_RUN = Path(sys.executable).resolve().parent.parent
else:
    JUST_RUN = ROOT.parent
LIB = ROOT / "lib"
VMA_K3 = JUST_RUN / "vma-k3"

for p in (str(LIB), str(VMA_K3)):
    if p not in sys.path:
        sys.path.insert(0, p)


def cmd_k3(args: argparse.Namespace) -> int:
    from vma.k3.cross_verifier import K3CrossVerifier

    data = args.input.encode() if args.text else Path(args.input).read_bytes()
    v = K3CrossVerifier(
        bits_ancho=args.bits,
        num_registros=args.registros,
        desfase_stride=args.desfase,
    )
    h = v.verificar_memoria_int(data)
    print(f"K3 [{args.modo}] bits={args.bits} registros={args.registros}")
    print(f"hash = 0x{h:08X}")
    print(f"bytes = {len(data)}")
    return 0


def cmd_mdc(args: argparse.Namespace) -> int:
    from mdc import mdc_factor

    n = int(args.n)
    f1, f2, info = mdc_factor(n, verbose=args.verbose)
    if f1:
        print(f"{n} = {f1} x {f2}  [{info['reason']}]")
    else:
        print(f"{n} -> candidato primo  [{info['reason']}]")
    return 0


def cmd_criva(args: argparse.Namespace) -> int:
    from criva import criva, compare_criva_vs_pnt

    if args.compare:
        xs = [int(x) for x in args.compare.split(",")]
        compare_criva_vs_pnt(xs)
        return 0
    x = float(args.x)
    est = criva(x) * x
    pnt = x / math.log(x) if x > 1 else 0
    print(f"x = {x:.0f}")
    print(f"Criva pi(x) ~ {est:.1f}")
    print(f"PNT   pi(x) ~ {pnt:.1f}")
    print(f"densidad Criva = {criva(x):.6f}")
    return 0


def cmd_criba(args: argparse.Namespace) -> int:
    from cribas import CribaModular6k

    limit = int(args.limit)
    primes = CribaModular6k(limit).run()
    print(f"Criba modular 6k±1 hasta {limit}: {len(primes)} primos")
    if args.show:
        print(primes[: args.show], "..." if len(primes) > args.show else "")
    return 0


def cmd_phase(args: argparse.Namespace) -> int:
    from k3_phase import run_demo

    print("Amplificador de fase K=3 XOR")
    print("S_A=(0,0,1,0)  S_B=(0,1,1,0)  — diff f1 oculta")
    print(f"{'t':>3}  {'S_A':>16}  {'S_B':>16}  dist")
    print("-" * 48)
    for t, sa, sb, dist in run_demo(args.steps):
        print(f"{t:>3}  {str(sa):>16}  {str(sb):>16}  {dist}")
    print("\nOK dist=2 desde t>=1 -> fase oculta amplificada a v2")
    return 0


def cmd_antipc(args: argparse.Namespace) -> int:
    antipc_dir = JUST_RUN / "antipc" / "src" / "antipc"
    antipc_src = antipc_dir / "udp_benchmark.py"
    if not antipc_src.exists():
        print(f"No encontrado: {antipc_src}")
        return 1
    import subprocess

    env = os.environ.copy()
    env["PYTHONPATH"] = str(antipc_dir) + os.pathsep + env.get("PYTHONPATH", "")
    cmd = [
        sys.executable,
        str(antipc_src),
        "--packets",
        str(args.packets),
        "--hubs",
        str(args.hubs),
    ]
    print(f"AntiPC UDP benchmark: {args.hubs} hubs, {args.packets} paquetes")
    return subprocess.call(cmd, cwd=str(antipc_dir), env=env)


def cmd_hurto(args: argparse.Namespace) -> int:
    csv_path = JUST_RUN / "hurto-gravitatorio" / "quijote_results.csv"
    if not csv_path.exists():
        print(f"No encontrado: {csv_path}")
        return 1
    print("Hurto gravitatorio — NREL 5MW (simulación)")
    print(f"{'N':>3} {'modo':<10} {'P_net':>8} {'eta':>6} {'dP%':>6}")
    print("-" * 40)
    with csv_path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            pnet = float(row["P_net_mean"])
            mark = "+" if pnet > 0 else " "
            print(
                f"{row['N']:>3} {row['mode']:<10} {mark}{pnet:>7.1f}kW "
                f"{float(row['eff_hurto']):>5.1f}x {float(row['millora']):>5.2f}%"
            )
    return 0


def cmd_demo(_: argparse.Namespace) -> int:
    print("=" * 60)
    print("VMA RUN — demo rápida (matemáticas + energía + K3)")
    print("=" * 60)

    print("\n[1/6] K3 auditoría")
    from vma.k3.cross_verifier import K3CrossVerifier

    h = K3CrossVerifier().verificar_memoria_int(b"VMA_JUST_RUN_2026")
    print(f"  hash('VMA_JUST_RUN_2026') = 0x{h:08X}")

    print("\n[2/6] MDC factorización")
    from mdc import mdc_factor

    for n in (10403, 15251, 100621):
        f1, f2, info = mdc_factor(n)
        print(f"  {n} = {f1}x{f2}" if f1 else f"  {n} primo? [{info['reason']}]")

    print("\n[3/6] Criva densidad primos")
    from criva import criva

    for x in (100, 1000, 10000):
        est = criva(x) * x
        print(f"  x={x:>5}  Criva pi(x)~{est:.0f}")

    print("\n[4/6] Criba 6k±1")
    from cribas import CribaModular6k

    from criva import sieve_primes

    criba_n = len(CribaModular6k(500).run())
    exact_n = len(sieve_primes(500))
    print(f"  criba 6k±1 hasta 500: {criba_n}  |  Eratóstenes: {exact_n}")

    print("\n[5/6] Amplificador fase K=3")
    from k3_phase import run_demo

    rows = run_demo(4)
    print(f"  t=1 dist={rows[1][3]}  t=3 dist={rows[3][3]}  (persistente)")

    print("\n[6/6] Hurto gravitatorio")
    cmd_hurto(argparse.Namespace())

    print("\n" + "=" * 60)
    print("Listo. Usa: python vma_run.py <k3|mdc|criva|criba|phase|hurto> --help")
    print("=" * 60)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="vma-run",
        description="VMA toolkit unificado — just run",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("demo", help="showcase rápido de todas las herramientas")
    d.set_defaults(func=cmd_demo)

    k = sub.add_parser("k3", help="hash auditoría industrial K3")
    k.add_argument("input", help="texto o ruta de archivo")
    k.add_argument("--text", action="store_true", help="input es texto literal")
    k.add_argument("--modo", default="usual")
    k.add_argument("--bits", type=int, default=32)
    k.add_argument("--registros", type=int, default=3)
    k.add_argument("--desfase", type=int, default=0)
    k.set_defaults(func=cmd_k3)

    m = sub.add_parser("mdc", help="factorizar N con MDC")
    m.add_argument("n", help="entero a factorizar")
    m.add_argument("-v", "--verbose", action="store_true")
    m.set_defaults(func=cmd_mdc)

    c = sub.add_parser("criva", help="estimar densidad π(x)/x")
    c.add_argument("x", nargs="?", default="1000", help="límite x")
    c.add_argument("--compare", metavar="X,Y,Z", help="tabla comparativa")
    c.set_defaults(func=cmd_criva)

    cr = sub.add_parser("criba", help="criba modular 6k±1")
    cr.add_argument("--limit", default="1000")
    cr.add_argument("--show", type=int, default=0, help="mostrar N primeros")
    cr.set_defaults(func=cmd_criba)

    ph = sub.add_parser("phase", help="teorema amplificador K=3 XOR")
    ph.add_argument("--steps", type=int, default=6)
    ph.set_defaults(func=cmd_phase)

    h = sub.add_parser("hurto", help="tabla hurto gravitatorio")
    h.set_defaults(func=cmd_hurto)

    ap = sub.add_parser("antipc", help="benchmark UDP AntiPC (hubs locales)")
    ap.add_argument("--hubs", type=int, default=4)
    ap.add_argument("--packets", type=int, default=500)
    ap.set_defaults(func=cmd_antipc)

    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())