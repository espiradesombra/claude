"""
Hash engine — HASHTOOLCODE k3hash (fuente oficial del ZIP del usuario).

Origen:
  C:\\Users\\cuent\\Desktop\\HASHTOOLCODE(l)L()(L ) (L).zip
  → k3hash/k3hash/src/k3hash.c

Usa DLL nativa si existe (native/k3hash/build/Release/k3hash.dll).
Si no, port Python bit-a-bit del mismo k3hash.c.
"""

from __future__ import annotations

from pathlib import Path

from native_engine import get_backend as _ne_backend
from native_engine import k3_hash_buffer as _ne_k3_hash_buffer

_ENGINE_DIR = Path(__file__).resolve().parent
_NATIVE_DIR = _ENGINE_DIR / "native" / "k3hash"


def get_backend() -> str:
    return _ne_backend()


def k3_hash_buffer(data: bytes, seed: int | None = None) -> int:
    """Hash principal AntiPC — k3hash del paquete HASHTOOLCODE."""
    return _ne_k3_hash_buffer(data, seed)


def k3_hash(data: bytes, seed: int | None = None) -> int:
    return k3_hash_buffer(data, seed)


def source_origin() -> str:
    candidate = _NATIVE_DIR / "src" / "k3hash.c"
    if candidate.is_file():
        return str(candidate)
    return str(_ENGINE_DIR / "native" / "k3hash" / "src" / "k3hash.c")