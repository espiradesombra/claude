#!/usr/bin/env python3
"""Demo KOP binario — ahorro de tiempo vs lookup JSON (gptcomputing K3 MicroKernel)."""

from __future__ import annotations

import json
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from runtime.kop import (
    KOPId,
    KOP_DEFINITIONS,
    K3MicroKernel,
    build_knowledge_blob,
    lookup_signature_from_blob,
    benchmark_lookup,
)
from runtime.kop_registry import KOPRegistry
from runtime.signature import make_signature


def main() -> None:
    payload = b"test-payload-antiPC-v011"
    sig = make_signature("HASH", [payload.hex()])
    ref_id = "ref-demo-001"

    blob = build_knowledge_blob(ref_id, sig, payload, producer="demo_kop")
    mk = K3MicroKernel.unpack(blob)

    print("=== AntiPC KOP binario (K3 MicroKernel) ===\n")
    print(f"Tamano blob     : {len(blob)} bytes")
    print(f"Magic           : {blob[:4]!r}")
    print(f"Footer          : {blob[-4:]!r}")
    print()

    print("Directorio (acceso directo):")
    for kop_id in KOPId:
        if kop_id in mk.blocks:
            name = KOP_DEFINITIONS[kop_id]["name"]
            data = mk.get_kop(kop_id)
            preview = data[:40] if data else b""
            print(f"  KOP {int(kop_id):03d} {name:12s}  {len(data or b'')} B  {preview!r}")

    print()
    print(f"Signature KOP003: {mk.get_signature_hex()}")
    print(f"Lookup directo  : {lookup_signature_from_blob(blob)}")
    print()
    print("KOP datos     :", ", ".join(KOPRegistry.list_data_kops()))
    print("KOP behavior  :", ", ".join(KOPRegistry.list_behavior_kops()))
    print()

    # Benchmark: JSON parse simulado vs salto directorio
    json_doc = json.dumps(
        {"kop_001": ref_id, "kop_003": sig, "payload": payload.hex()},
        separators=(",", ":"),
    )

    def json_lookup() -> str:
        doc = json.loads(json_doc)
        return doc["kop_003"]

    stats = benchmark_lookup(json_lookup, blob, n=100000)
    print("Benchmark lookup x100000:")
    print(f"  JSON parse     : {stats['json_sec']:.4f} s")
    print(f"  KOP directorio : {stats['binary_sec']:.4f} s")
    print(f"  Aceleracion    : {stats['speedup_x']:.1f}x")
    print(f"  Tiempo ahorrado: {stats['saved_pct']:.1f}%")
    print()
    print("Norma gptcomputing: un KOP nunca modifica otro; Identity inmutable.")
    print("Siguiente: integrar blob en KnowledgeBuffer.publish (v0.1.1).")


if __name__ == "__main__":
    main()