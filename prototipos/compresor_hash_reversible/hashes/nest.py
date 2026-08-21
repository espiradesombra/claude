"""Anidación de familias: cada capa hashea el blob de la anterior."""

from __future__ import annotations

import json
from typing import Callable

from .cinematic import hash_cinematic
from .espejo import hash_espejo
from .k3_core import hash_k3_plain, k3_bytes
from .toffoli import hash_toffoli

FAMILIES: dict[str, Callable[[bytes], dict]] = {
    "cinematic": hash_cinematic,
    "espejo": hash_espejo,
    "toffoli": hash_toffoli,
    "toffoli_k3": hash_toffoli,
    "k3": hash_k3_plain,
}

DEFAULT_PIPELINE = ["cinematic", "espejo", "toffoli", "k3"]


def _layer_blob(layer: dict) -> bytes:
    """Serialización estable para alimentar la siguiente capa."""
    meta = layer.get("meta") or {}
    # compact: family|digest|json_meta
    return (
        f"{layer['family']}|{layer['digest_hex']}|".encode("utf-8")
        + json.dumps(meta, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )


def nest_pipeline(data: bytes, pipeline: list[str] | None = None) -> dict:
    """
    Aplica familias en orden. La capa 0 ve el bloque crudo;
    cada siguiente ve el blob (digest+meta) de la anterior → ANIDACIÓN.
    """
    pipeline = pipeline or list(DEFAULT_PIPELINE)
    layers: list[dict] = []
    current = data
    for name in pipeline:
        fn = FAMILIES.get(name)
        if fn is None:
            raise ValueError(f"familia desconocida: {name}; opciones={list(FAMILIES)}")
        layer = fn(current)
        layers.append(
            {
                "family": layer["family"],
                "reversible": layer.get("reversible", False),
                "digest_hex": layer["digest_hex"],
                "digest": layer["digest"],
                "meta": layer.get("meta", {}),
            }
        )
        current = _layer_blob(layer)

    # digest exterior = K3 de la concatenación de digests de capas
    outer = b"".join(int(L["digest"]).to_bytes(4, "big") for L in layers)
    outer_d = k3_bytes(outer)
    return {
        "pipeline": pipeline,
        "layers": [
            {k: v for k, v in L.items() if k != "digest"}  # digest int no JSON-friendly dup
            for L in layers
        ],
        "outer_digest": outer_d,
        "outer_digest_hex": f"{outer_d:08x}",
        "n_layers": len(layers),
    }
