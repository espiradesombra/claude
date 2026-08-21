"""Familias de hash VMA anidables en el compresor VHC."""

from __future__ import annotations

from .cinematic import hash_cinematic
from .espejo import hash_espejo
from .k3_core import k3_bytes, hash_k3_plain
from .toffoli import hash_toffoli
from .nest import nest_pipeline, FAMILIES

__all__ = [
    "hash_cinematic",
    "hash_espejo",
    "hash_k3_plain",
    "hash_toffoli",
    "k3_bytes",
    "nest_pipeline",
    "FAMILIES",
]
