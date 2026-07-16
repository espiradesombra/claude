#!/usr/bin/env python3
"""
Demo Microkernel AntiPC v0.1.2-alpha

Flujo: Boot → Registry → IdentityKOP → Ledger → Lookup → Reuse
Inspirado en gptcomputing RunPlan; implementacion propia Python.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from runtime.microkernel import AntiPCMicroKernel


def main() -> None:
    mk = AntiPCMicroKernel()
    mk.boot()

    mk.run_boot_kop()

    payload = b"Hola"
    oid = mk.create_knowledge(payload)
    if oid is None:
        sys.exit(1)

    print(f"Lookup {oid}..........{'OK' if mk.lookup(oid) else 'FAIL'}")

    print("Reuse (2nd IdentityKOP)...")
    oid2 = mk.create_knowledge(payload)
    reused = oid2 == oid or mk.metrics.kop_reused > 0
    print(f"Reuse................{'OK' if reused else 'FAIL'}")

    mk.shutdown()


if __name__ == "__main__":
    main()