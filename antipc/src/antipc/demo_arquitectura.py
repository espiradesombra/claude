#!/usr/bin/env python3
"""Demo arquitectura unificada AntiPC + HASHTOOLCODE."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from architecture.antipc_stack import AntiPCStack
from runtime.profile import ExecutionProfile


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(Path(__file__).parent.parent / "sale-it"))
    parser.add_argument("--n", type=int, default=1147, help="semiprimo toy MDC (31×37)")
    args = parser.parse_args()

    stack = AntiPCStack(profile=ExecutionProfile.industrial_audit())

    print()
    print("=" * 72)
    print("  AntiPC — ARQUITECTURA UNIFICADA")
    print("  HASHTOOLCODE (K3+HMAC) + Regla Mecánica (MDC) + Flow Kernel + Red")
    print("=" * 72)

    for layer in stack.architecture_map():
        print(f"\n  [{layer.layer.name}] {layer.name}")
        print(f"    {layer.source}")
        for c in layer.components:
            print(f"      · {c}")

    print("\n" + "-" * 72)
    print("  PERMISOS DE RED (HMAC challenge-response)")
    ok = stack.request_network_permission()
    print(f"    Autorización local: {'OK' if ok else 'DENEGADA'}")

    if ok:
        print("\n" + "-" * 72)
        print("  ENLACE L3 — UDP HUBS + HMAC")
        try:
            net = stack.start_network(hubs=2)
            print(f"    Hubs autenticados: {net['authenticated']}/{net['total']}")
            payloads = [b"AntiPC-L3-test", b"HASHTOOLCODE-linked"]
            digests = stack.dispatch_remote_hash(payloads)
            for seq, d in digests.items():
                print(f"    seq={seq} digest={d:08X} (hub remoto)")
            stack.stop_network()
            print("    Enlace OK — maestro → HMAC → hub → K3 → RESULT")
        except Exception as exc:
            print(f"    Enlace fallido: {exc}")
            stack.stop_network()

    print("\n" + "-" * 72)
    print("  PIPELINE MDC (toy-number)")
    mdc = stack.pipeline_analyze_number(args.n)
    print(f"    N={mdc['n']}  factor={mdc['factor']}  regla={mdc.get('regla_producto')}")

    root = Path(args.root)
    if root.is_dir():
        sample = next(root.rglob("*.py"), None)
        if sample:
            print("\n" + "-" * 72)
            print("  PIPELINE INDUSTRIAL (fichero)")
            audit = stack.pipeline_audit_file(str(sample))
            print(f"    {audit['path']}")
            print(f"    K3={audit['k3_digest']}  MDC_fase={audit.get('mdc_fase', 0):.4f}")

    out = Path(__file__).parent / "output" / "ARQUITECTURA_UNIFICADA.txt"
    stack.export_architecture_report(out)
    print(f"\n  Informe: {out}")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())