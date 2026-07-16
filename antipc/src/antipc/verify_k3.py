#!/usr/bin/env python3
"""Verifica hash HASHTOOLCODE k3hash — Python port vs DLL nativa."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hash_engine import get_backend, k3_hash_buffer, source_origin
from k3_hash import k3_hash_buffer as py_only


def main() -> int:
    samples = [
        b"",
        b"AntiPC",
        b"x^=(B*0x9E3779B1);",
        b"HASHTOOLCODE",
        os.urandom(128),
    ]

    print()
    print("=" * 60)
    print("  Verificacion K3 — HASHTOOLCODE zip")
    print("=" * 60)
    print(f"  Fuente C   : {source_origin()}")
    print(f"  Backend    : {get_backend()}")
    print()

    ok = True
    for sample in samples:
        engine_h = k3_hash_buffer(sample)
        py_h = py_only(sample)
        match = engine_h == py_h
        ok = ok and match
        tag = "OK" if match else "MISMATCH"
        preview = sample[:16].hex() if sample else "(empty)"
        print(f"  [{tag}] engine={engine_h:08X}  python={py_h:08X}  data={preview}")

    print()
    if get_backend().startswith("native"):
        print("  Usando DLL nativa del paquete HASHTOOLCODE.")
    else:
        print("  Sin DLL: port Python del mismo k3hash.c del ZIP.")
        print("  Para compilar DLL: native\\k3hash\\build_k3hash.bat")
    print("=" * 60)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())