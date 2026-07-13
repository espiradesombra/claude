"""
Hash engine — HASHTOOLCODE k3hash (fuente oficial del ZIP del usuario).

Origen:
  C:\\Users\\cuent\\Desktop\\HASHTOOLCODE(l)L()(L ) (L).zip
  → k3hash/k3hash/src/k3hash.c

Usa DLL nativa si existe (native/k3hash/build/Release/k3hash.dll).
Si no, port Python bit-a-bit del mismo k3hash.c.
"""

from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path

_ENGINE_DIR = Path(__file__).resolve().parent
_NATIVE_DIR = _ENGINE_DIR / "native" / "k3hash"
_DLL_CANDIDATES = [
    _NATIVE_DIR / "build" / "Release" / "k3hash.dll",
    _NATIVE_DIR / "build" / "k3hash.dll",
    _ENGINE_DIR / "k3hash.dll",
]

_native_lib: ctypes.CDLL | None = None
_backend: str = "python"


def _load_native() -> ctypes.CDLL | None:
    global _native_lib, _backend
    if _native_lib is not None:
        return _native_lib

    for dll_path in _DLL_CANDIDATES:
        if not dll_path.is_file():
            continue
        try:
            lib = ctypes.CDLL(str(dll_path))

            class K3HashConfig(ctypes.Structure):
                _fields_ = [
                    ("bits_ancho", ctypes.c_int),
                    ("num_registros", ctypes.c_int),
                    ("semilla_inicial", ctypes.c_uint32),
                ]

            lib.k3_config_default.restype = K3HashConfig
            lib.k3_hash_buffer.argtypes = [
                ctypes.c_void_p,
                ctypes.c_size_t,
                ctypes.POINTER(K3HashConfig),
            ]
            lib.k3_hash_buffer.restype = ctypes.c_uint32
            lib._K3HashConfig = K3HashConfig  # type: ignore[attr-defined]

            _native_lib = lib
            _backend = f"native:{dll_path.name}"
            return lib
        except OSError:
            continue
    return None


def get_backend() -> str:
    _load_native()
    return _backend


def k3_hash_buffer(data: bytes, seed: int | None = None) -> int:
    """Hash principal AntiPC — k3hash del paquete HASHTOOLCODE."""
    lib = _load_native()
    if lib is not None:
        cfg = lib.k3_config_default()
        if seed is not None:
            cfg.semilla_inicial = seed & 0xFFFFFFFF
        buf = ctypes.create_string_buffer(data, len(data))
        return int(lib.k3_hash_buffer(buf, len(data), ctypes.byref(cfg)))

    from k3_hash import K3HashConfig, k3_hash_buffer as py_hash

    config = K3HashConfig()
    if seed is not None:
        config.semilla_inicial = seed & 0xFFFFFFFF
    return py_hash(data, config)


def k3_hash(data: bytes, seed: int | None = None) -> int:
    return k3_hash_buffer(data, seed)


def source_origin() -> str:
    return str(_NATIVE_DIR / "src" / "k3hash.c")