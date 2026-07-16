"""
Motor nativo unificado AntiPC — antipc_native.dll (C).

Integra: k3hash, MDC, criba Eratóstenes, Criva, stream K3, convergencia geométrica.
Fuentes en PC: encriptacionGeometrica/, vma-methods/, mdc_lib/, native/k3hash/.
"""

from __future__ import annotations

import ctypes
import sys
import time
from pathlib import Path

_ENGINE_DIR = Path(__file__).resolve().parent
_ANTIPC_ROOT = _ENGINE_DIR.parent.parent

_lib: ctypes.CDLL | None = None
_backend: str = "python"
_loaded_from: Path | None = None


K3HASH_MAX_REGISTROS = 10
K3_REDUNDANT_SEEDS = (0xA5A5A5A5, 0x5A5A5A5A, 0x12345678)


class K3HashConfig(ctypes.Structure):
    _fields_ = [
        ("bits_ancho", ctypes.c_int),
        ("num_registros", ctypes.c_int),
        ("semilla_inicial", ctypes.c_uint32),
    ]


class K3HashCtx(ctypes.Structure):
    _fields_ = [
        ("config", K3HashConfig),
        ("registros", ctypes.c_uint32 * K3HASH_MAX_REGISTROS),
        ("buffer_parcial", ctypes.c_uint8 * 4),
        ("bytes_en_buffer", ctypes.c_int),
        ("ancho_bytes", ctypes.c_int),
    ]


class AntipcMdcCollision(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_int32),
        ("y", ctypes.c_int32),
        ("s", ctypes.c_uint32),
        ("t", ctypes.c_uint32),
        ("k", ctypes.c_uint32),
    ]


class AntipcMdcTrainResult(ctypes.Structure):
    _fields_ = [
        ("steps", ctypes.c_uint32),
        ("n_collisions", ctypes.c_uint32),
        ("hits", AntipcMdcCollision * 64),
    ]


class AntipcGeoClave(ctypes.Structure):
    _fields_ = [
        ("tales_count", ctypes.c_uint32),
        ("tales", ctypes.c_uint32 * 8),
        ("figuras_count", ctypes.c_uint32),
        ("figuras", ctypes.c_int32 * 8),
        ("puntos_count", ctypes.c_uint32),
        ("puntos", ctypes.c_uint32 * 8),
        ("saltos_count", ctypes.c_uint32),
        ("saltos", ctypes.c_int32 * 8),
        ("primos_count", ctypes.c_uint32),
        ("primos_p1", ctypes.c_int32 * 4),
        ("primos_p2", ctypes.c_int32 * 4),
        ("porcentajes", ctypes.c_int32 * 4),
        ("iteraciones_pi", ctypes.c_int32),
    ]


class AntipcNewtonResult(ctypes.Structure):
    _fields_ = [
        ("j", ctypes.c_double),
        ("j_exacto", ctypes.c_double),
        ("error", ctypes.c_double),
        ("iterations", ctypes.c_int),
        ("converged", ctypes.c_int),
    ]


NEWTON_FAMILIA = {
    "general": 0,
    "cuadrados": 1,
    "cubos": 2,
    "potencia": 3,
    "kp": 4,
    "mersenne": 5,
}


def _dll_candidates() -> list[Path]:
    paths = [
        _ENGINE_DIR / "antipc_native.dll",
        _ENGINE_DIR / "native" / "antipc_core" / "build" / "Release" / "antipc_native.dll",
        _ANTIPC_ROOT / "antipc_native.dll",
        _ANTIPC_ROOT / "dist" / "antipc_native.dll",
        _ENGINE_DIR / "native" / "k3hash" / "build" / "Release" / "k3hash.dll",
        _ENGINE_DIR / "k3hash.dll",
        _ANTIPC_ROOT / "k3hash.dll",
    ]
    if getattr(sys, "frozen", False):
        exe_dir = Path(sys.executable).resolve().parent
        paths.insert(0, exe_dir / "antipc_native.dll")
        paths.insert(1, exe_dir / "k3hash.dll")
    return paths


def _bind_antipc_native(lib: ctypes.CDLL) -> None:
    lib.antipc_native_version.restype = ctypes.c_char_p

    lib.antipc_mdc_factor.argtypes = [ctypes.c_uint64]
    lib.antipc_mdc_factor.restype = ctypes.c_uint64

    lib.antipc_sieve_count.argtypes = [ctypes.c_uint32]
    lib.antipc_sieve_count.restype = ctypes.c_uint32

    if hasattr(lib, "antipc_sieve_hibrida_count"):
        lib.antipc_sieve_hibrida_count.argtypes = [ctypes.c_uint32]
        lib.antipc_sieve_hibrida_count.restype = ctypes.c_uint32

    if hasattr(lib, "antipc_sieve_modular6k_count"):
        lib.antipc_sieve_modular6k_count.argtypes = [ctypes.c_uint32]
        lib.antipc_sieve_modular6k_count.restype = ctypes.c_uint32

    if hasattr(lib, "antipc_mdc_scan_trains"):
        lib.antipc_mdc_scan_trains.argtypes = [
            ctypes.c_uint64,
            ctypes.POINTER(AntipcMdcTrainResult),
            ctypes.POINTER(AntipcMdcTrainResult),
        ]
        lib.antipc_mdc_scan_trains.restype = None

    if hasattr(lib, "antipc_newton_oraculo"):
        lib.antipc_newton_oraculo.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_double,
        ]
        lib.antipc_newton_oraculo.restype = ctypes.c_double
        lib.antipc_newton_rapido.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
        ]
        lib.antipc_newton_rapido.restype = AntipcNewtonResult
        lib.antipc_newton_log.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_double,
        ]
        lib.antipc_newton_log.restype = AntipcNewtonResult

    if hasattr(lib, "antipc_geo_clave_default"):
        lib.antipc_geo_clave_default.argtypes = [ctypes.POINTER(AntipcGeoClave)]
        lib.antipc_geo_clave_default.restype = None
        lib.antipc_geo_fase_encrypt.argtypes = [ctypes.c_uint16, ctypes.POINTER(AntipcGeoClave)]
        lib.antipc_geo_fase_encrypt.restype = ctypes.c_double
        lib.antipc_geo_fase_decrypt.argtypes = [
            ctypes.c_double,
            ctypes.POINTER(AntipcGeoClave),
            ctypes.POINTER(ctypes.c_uint16),
        ]
        lib.antipc_geo_fase_decrypt.restype = ctypes.c_int
        lib.antipc_aleatorovix_xor.argtypes = [
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_size_t,
            ctypes.c_uint16,
        ]
        lib.antipc_aleatorovix_xor.restype = None
        lib.antipc_geo_masivo_crypt.argtypes = [
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_size_t,
            ctypes.c_uint16,
            ctypes.POINTER(AntipcGeoClave),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
        ]
        lib.antipc_geo_masivo_crypt.restype = ctypes.c_int

    if hasattr(lib, "antipc_mdc_ksweep_classic"):
        lib.antipc_mdc_ksweep_classic.argtypes = [
            ctypes.c_uint64,
            ctypes.c_uint64,
            ctypes.c_uint64,
        ]
        lib.antipc_mdc_ksweep_classic.restype = ctypes.c_uint64
        lib.antipc_mdc_ksweep_predict.argtypes = [
            ctypes.c_uint64,
            ctypes.c_uint64,
            ctypes.c_uint64,
            ctypes.POINTER(ctypes.c_uint64),
        ]
        lib.antipc_mdc_ksweep_predict.restype = ctypes.c_uint64

    lib.antipc_sieve_fill.argtypes = [
        ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_uint32,
    ]
    lib.antipc_sieve_fill.restype = ctypes.c_uint32

    lib.antipc_criva_estimate.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.c_int]
    lib.antipc_criva_estimate.restype = ctypes.c_double

    lib.antipc_k3_stream_xor.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.c_uint64,
        ctypes.c_uint64,
    ]
    lib.antipc_k3_stream_xor.restype = None

    lib.antipc_geo_converge.argtypes = [
        ctypes.c_char_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_size_t,
    ]
    lib.antipc_geo_converge.restype = ctypes.c_uint64

    lib.k3_config_default.restype = K3HashConfig
    lib.k3_hash_buffer.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(K3HashConfig),
    ]
    lib.k3_hash_buffer.restype = ctypes.c_uint32

    if hasattr(lib, "k3_hash_init"):
        lib.k3_hash_init.argtypes = [ctypes.POINTER(K3HashCtx), ctypes.POINTER(K3HashConfig)]
        lib.k3_hash_init.restype = ctypes.c_int
        lib.k3_hash_update.argtypes = [
            ctypes.POINTER(K3HashCtx),
            ctypes.c_void_p,
            ctypes.c_size_t,
        ]
        lib.k3_hash_update.restype = ctypes.c_int
        lib.k3_hash_final.argtypes = [ctypes.POINTER(K3HashCtx)]
        lib.k3_hash_final.restype = ctypes.c_uint32

    if hasattr(lib, "k3_hash_file"):
        lib.k3_hash_file.argtypes = [
            ctypes.c_char_p,
            ctypes.POINTER(K3HashConfig),
            ctypes.POINTER(ctypes.c_uint32),
        ]
        lib.k3_hash_file.restype = ctypes.c_int

    if hasattr(lib, "antipc_k3_hash_file"):
        lib.antipc_k3_hash_file.argtypes = [
            ctypes.c_char_p,
            ctypes.c_uint32,
            ctypes.POINTER(ctypes.c_uint32),
        ]
        lib.antipc_k3_hash_file.restype = ctypes.c_int
        lib.antipc_k3_fingerprint_file.argtypes = [
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int64),
            ctypes.POINTER(ctypes.c_uint32),
        ]
        lib.antipc_k3_fingerprint_file.restype = ctypes.c_int
        lib.antipc_k3_heavy_hash.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        lib.antipc_k3_heavy_hash.restype = ctypes.c_uint32
        lib.antipc_k3_redundant_hashes.argtypes = [
            ctypes.c_void_p,
            ctypes.c_size_t,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.c_int,
        ]
        lib.antipc_k3_redundant_hashes.restype = ctypes.c_int
        lib.antipc_k3_hamming.argtypes = [ctypes.c_uint32, ctypes.c_uint32]
        lib.antipc_k3_hamming.restype = ctypes.c_int
        lib.antipc_k3_similarity.argtypes = [ctypes.c_uint32, ctypes.c_uint32]
        lib.antipc_k3_similarity.restype = ctypes.c_double

    lib._K3HashConfig = K3HashConfig  # type: ignore[attr-defined]


def _bind_k3hash_only(lib: ctypes.CDLL) -> None:
    lib.k3_config_default.restype = K3HashConfig
    lib.k3_hash_buffer.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(K3HashConfig),
    ]
    lib.k3_hash_buffer.restype = ctypes.c_uint32
    if hasattr(lib, "k3_hash_file"):
        lib.k3_hash_file.argtypes = [
            ctypes.c_char_p,
            ctypes.POINTER(K3HashConfig),
            ctypes.POINTER(ctypes.c_uint32),
        ]
        lib.k3_hash_file.restype = ctypes.c_int
    lib._K3HashConfig = K3HashConfig  # type: ignore[attr-defined]


def load_native(force: bool = False) -> ctypes.CDLL | None:
    global _lib, _backend, _loaded_from
    if _lib is not None and not force:
        return _lib

    for dll_path in _dll_candidates():
        if not dll_path.is_file():
            continue
        try:
            lib = ctypes.CDLL(str(dll_path))
            if hasattr(lib, "antipc_native_version"):
                _bind_antipc_native(lib)
                ver = lib.antipc_native_version().decode("ascii", errors="replace")
                _backend = f"antipc_native:{ver}@{dll_path.name}"
            else:
                _bind_k3hash_only(lib)
                _backend = f"k3hash-only:{dll_path.name}"
            _lib = lib
            _loaded_from = dll_path
            return lib
        except OSError:
            continue
    return None


def get_backend() -> str:
    load_native()
    return _backend


def is_full_native() -> bool:
    load_native()
    return _lib is not None and hasattr(_lib, "antipc_mdc_factor")


def k3_hash_buffer(data: bytes, seed: int | None = None) -> int:
    lib = load_native()
    if lib is not None:
        cfg = lib.k3_config_default()
        if seed is not None:
            cfg.semilla_inicial = seed & 0xFFFFFFFF
        buf = ctypes.create_string_buffer(data, len(data))
        return int(lib.k3_hash_buffer(buf, len(data), ctypes.byref(cfg)))
    from k3_hash import K3HashConfig, k3_hash_buffer as py_hash

    cfg = K3HashConfig()
    if seed is not None:
        cfg.semilla_inicial = seed & 0xFFFFFFFF
    return py_hash(data, cfg)


def k3_hash_stream(data: bytes, seed: int | None = None) -> int:
    """Hash por streaming (K3HashCtx) — equivalente a buffer para datos completos."""
    lib = load_native()
    if lib is not None and hasattr(lib, "k3_hash_init"):
        cfg = lib.k3_config_default()
        if seed is not None:
            cfg.semilla_inicial = seed & 0xFFFFFFFF
        ctx = K3HashCtx()
        lib.k3_hash_init(ctypes.byref(ctx), ctypes.byref(cfg))
        if data:
            buf = ctypes.create_string_buffer(data, len(data))
            lib.k3_hash_update(ctypes.byref(ctx), buf, len(data))
        return int(lib.k3_hash_final(ctypes.byref(ctx)))
    return k3_hash_buffer(data, seed)


def k3_hash_file(path: str | Path, seed: int | None = None) -> int:
    """Hash de fichero por streaming C (sin cargar todo en RAM)."""
    lib = load_native()
    p = str(Path(path).resolve())
    if lib is not None and hasattr(lib, "antipc_k3_hash_file"):
        out = ctypes.c_uint32(0)
        s = 0 if seed is None else (seed & 0xFFFFFFFF)
        if lib.antipc_k3_hash_file(p.encode("utf-8"), s, ctypes.byref(out)) != 0:
            raise OSError(f"No se pudo hashear: {p}")
        return int(out.value)
    if lib is not None and hasattr(lib, "k3_hash_file"):
        cfg = lib.k3_config_default()
        if seed is not None:
            cfg.semilla_inicial = seed & 0xFFFFFFFF
        out = ctypes.c_uint32(0)
        if lib.k3_hash_file(p.encode("utf-8"), ctypes.byref(cfg), ctypes.byref(out)) != 0:
            raise OSError(f"No se pudo hashear: {p}")
        return int(out.value)
    return k3_hash_buffer(Path(p).read_bytes(), seed)


def k3_fingerprint_file(path: str | Path) -> tuple[int, int]:
    """Retorna (tamaño_bytes, hash) — lógica k3dedup."""
    lib = load_native()
    p = str(Path(path).resolve())
    if lib is not None and hasattr(lib, "antipc_k3_fingerprint_file"):
        size = ctypes.c_int64(0)
        digest = ctypes.c_uint32(0)
        if lib.antipc_k3_fingerprint_file(
            p.encode("utf-8"), ctypes.byref(size), ctypes.byref(digest)
        ) != 0:
            raise OSError(f"No se pudo fingerprint: {p}")
        return int(size.value), int(digest.value)
    size = Path(p).stat().st_size
    return size, k3_hash_file(p)


def k3_redundant_hashes(data: bytes, replicas: int = 3) -> list[int]:
    """Hashes redundantes Grafcet (semillas 0xA5A5A5A5, 0x5A5A5A5A, 0x12345678)."""
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_k3_redundant_hashes"):
        arr = (ctypes.c_uint32 * replicas)()
        n = int(
            lib.antipc_k3_redundant_hashes(
                data,
                len(data),
                ctypes.c_int(replicas),
                arr,
                ctypes.c_int(replicas),
            )
        )
        return [int(arr[i]) for i in range(n)]
    return [
        k3_hash_buffer(data, seed=K3_REDUNDANT_SEEDS[i % len(K3_REDUNDANT_SEEDS)])
        for i in range(replicas)
    ]


def k3_heavy_hash(data: bytes) -> int:
    """Cadena 4 rondas Grafcet: hash(digest_le || payload)."""
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_k3_heavy_hash"):
        buf = ctypes.create_string_buffer(data, len(data)) if data else None
        return int(lib.antipc_k3_heavy_hash(buf, len(data)))
    digest = k3_hash_buffer(data)
    for _ in range(4):
        digest = k3_hash_buffer(digest.to_bytes(4, "little") + data)
    return digest


def k3_hamming(a: int, b: int) -> int:
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_k3_hamming"):
        return int(lib.antipc_k3_hamming(ctypes.c_uint32(a), ctypes.c_uint32(b)))
    x = (a ^ b) & 0xFFFFFFFF
    count = 0
    while x:
        count += x & 1
        x >>= 1
    return count


def k3_similarity(a: int, b: int) -> float:
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_k3_similarity"):
        return float(lib.antipc_k3_similarity(ctypes.c_uint32(a), ctypes.c_uint32(b)))
    return 1.0 - k3_hamming(a, b) / 32.0


def k3_dedup_paths(paths: list[str | Path]) -> dict:
    """Agrupa duplicados por (tamaño, hash) — k3dedup.c."""
    entries: list[tuple[str, int, int]] = []
    for raw in paths:
        p = Path(raw)
        if not p.is_file():
            continue
        size, digest = k3_fingerprint_file(p)
        entries.append((str(p.resolve()), size, digest))

    entries.sort(key=lambda e: (e[1], e[2]))

    groups: list[dict] = []
    i = 0
    recoverable = 0
    while i < len(entries):
        j = i + 1
        while (
            j < len(entries)
            and entries[j][1] == entries[i][1]
            and entries[j][2] == entries[i][2]
        ):
            j += 1
        if j - i > 1:
            paths_in_group = [entries[k][0] for k in range(i, j)]
            groups.append(
                {
                    "size": entries[i][1],
                    "hash": f"{entries[i][2]:08X}",
                    "copies": len(paths_in_group),
                    "paths": paths_in_group,
                }
            )
            recoverable += entries[i][1] * (len(paths_in_group) - 1)
        i = j

    return {
        "files_analyzed": len(entries),
        "duplicate_groups": len(groups),
        "bytes_recoverable": recoverable,
        "groups": groups,
    }


def sieve_hibrida_count(limit: int) -> int:
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_sieve_hibrida_count"):
        return int(lib.antipc_sieve_hibrida_count(ctypes.c_uint32(limit)))
    from vma_methods_bridge import sieve_count_python

    return sieve_count_python(limit)


def sieve_modular6k_count(limit: int) -> int:
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_sieve_modular6k_count"):
        return int(lib.antipc_sieve_modular6k_count(ctypes.c_uint32(limit)))
    return sieve_count(limit)


def sieve_desmemoriada_count(limit: int) -> int:
    from vma_methods_bridge import sieve_desmemoriada_count as _py

    return _py(limit)


def mdc_scan_trains(n: int) -> tuple[AntipcMdcTrainResult, AntipcMdcTrainResult]:
    lib = load_native()
    if lib is None or not hasattr(lib, "antipc_mdc_scan_trains"):
        raise RuntimeError("antipc_mdc_scan_trains no disponible")
    tx = AntipcMdcTrainResult()
    ty = AntipcMdcTrainResult()
    lib.antipc_mdc_scan_trains(ctypes.c_uint64(n), ctypes.byref(tx), ctypes.byref(ty))
    return tx, ty


def newton_log(
    E: float,
    b: float = 10.0,
    familia: str = "general",
    *,
    n_exp: int = 2,
    k_known: float = 1.0,
) -> dict:
    lib = load_native()
    if lib is None or not hasattr(lib, "antipc_newton_log"):
        from vma_methods_bridge import newton_log_python

        return newton_log_python(E, b, familia, n_exp=n_exp, k_known=k_known)
    fam = NEWTON_FAMILIA.get(familia.lower(), 0)
    r = lib.antipc_newton_log(
        float(E), float(b), ctypes.c_int(fam), ctypes.c_int(n_exp), float(k_known)
    )
    return {
        "j": r.j,
        "j_exacto": r.j_exacto,
        "iteraciones": r.iterations,
        "error": r.error,
        "converged": bool(r.converged),
        "backend": _backend,
    }


def mdc_ksweep(n: int, *, predict: bool = True) -> tuple[int, int]:
    """K-sweep en rango automático (toy). Retorna (factor_D, evaluaciones)."""
    import math

    lib = load_native()
    if lib is None:
        raise RuntimeError("antipc_native.dll no disponible")

    sq = math.isqrt(n)
    m_max = (sq - 3) // 2
    if m_max < 1:
        return 0, 0
    m_ini = 1
    evals = ctypes.c_uint64(0)

    if predict and hasattr(lib, "antipc_mdc_ksweep_predict"):
        f = int(
            lib.antipc_mdc_ksweep_predict(
                ctypes.c_uint64(n),
                ctypes.c_uint64(m_ini),
                ctypes.c_uint64(m_max),
                ctypes.byref(evals),
            )
        )
        return f, int(evals.value)

    if hasattr(lib, "antipc_mdc_ksweep_classic"):
        f = int(
            lib.antipc_mdc_ksweep_classic(
                ctypes.c_uint64(n),
                ctypes.c_uint64(m_ini),
                ctypes.c_uint64(m_max),
            )
        )
        return f, 0
    return 0, 0


def mdc_factor(n: int) -> int:
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_mdc_factor"):
        return int(lib.antipc_mdc_factor(ctypes.c_uint64(n)))
    from mdc_lib.factoritzacio import factorizar_mdc_toy

    return factorizar_mdc_toy(n)


def sieve_count(limit: int) -> int:
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_sieve_count"):
        return int(lib.antipc_sieve_count(ctypes.c_uint32(limit)))
    from vma_methods_bridge import sieve_count_python

    return sieve_count_python(limit)


def sieve_primes(limit: int, cap: int | None = None) -> list[int]:
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_sieve_fill"):
        max_cap = cap or max(64, limit // 10)
        arr = (ctypes.c_uint32 * max_cap)()
        n = int(lib.antipc_sieve_fill(ctypes.c_uint32(limit), arr, ctypes.c_uint32(max_cap)))
        return [int(arr[i]) for i in range(n)]
    from vma_methods_bridge import sieve_primes_python

    return sieve_primes_python(limit)[: cap or None]


def criva_estimate(x: float, layers: int = 10, iterations: int = 8) -> float:
    lib = load_native()
    if lib is not None and hasattr(lib, "antipc_criva_estimate"):
        return float(lib.antipc_criva_estimate(x, layers, iterations))
    from vma_methods_bridge import criva_python

    return criva_python(x, layers, iterations)


def k3_stream_xor(data: bytearray, base: int = 33, rel: int = 1) -> None:
    lib = load_native()
    if lib is None or not hasattr(lib, "antipc_k3_stream_xor"):
        raise RuntimeError("antipc_native.dll no disponible para k3_stream_xor")
    buf = (ctypes.c_ubyte * len(data)).from_buffer(data)
    lib.antipc_k3_stream_xor(buf, len(data), ctypes.c_uint64(base), ctypes.c_uint64(rel))


def geo_clave_default() -> AntipcGeoClave:
    lib = load_native()
    if lib is None or not hasattr(lib, "antipc_geo_clave_default"):
        raise RuntimeError("antipc_geo_clave_default no disponible")
    clave = AntipcGeoClave()
    lib.antipc_geo_clave_default(ctypes.byref(clave))
    return clave


def geo_masivo_crypt(
    data: bytes,
    semilla: int = 43210,
    *,
    decrypt: bool = False,
    hash_fase: float | None = None,
    clave: AntipcGeoClave | None = None,
) -> tuple[bytes, float]:
    """Cifrado/descifrado masivo Aleatorovix + fase geométrica."""
    lib = load_native()
    if lib is None or not hasattr(lib, "antipc_geo_masivo_crypt"):
        raise RuntimeError("antipc_geo_masivo_crypt no disponible")

    if clave is None:
        clave = geo_clave_default()

    out = ctypes.create_string_buffer(len(data))
    hash_out = ctypes.c_double(hash_fase or 0.0)
    hash_ptr = ctypes.byref(hash_out) if decrypt or hash_fase is None else ctypes.byref(hash_out)

    ok = lib.antipc_geo_masivo_crypt(
        data,
        out,
        len(data),
        ctypes.c_uint16(semilla & 0xFFFF),
        ctypes.byref(clave),
        hash_ptr,
        1 if decrypt else 0,
    )
    if not ok:
        raise RuntimeError("geo_masivo_crypt falló (backtracking o parámetros)")
    return out.raw, float(hash_out.value)


def geo_converge(bits: str, tales: list[int], puntos: list[int]) -> int:
    lib = load_native()
    if lib is None or not hasattr(lib, "antipc_geo_converge"):
        raise RuntimeError("antipc_native.dll no disponible para geo_converge")
    t_arr = (ctypes.c_uint32 * len(tales))(*tales)
    p_arr = (ctypes.c_uint32 * len(puntos))(*puntos)
    raw = bits.encode("ascii")
    return int(
        lib.antipc_geo_converge(
            raw,
            len(raw),
            t_arr,
            len(tales),
            p_arr,
            len(puntos),
        )
    )


def status_report() -> str:
    load_native()
    lines = [
        f"Backend     : {_backend}",
        f"DLL path    : {_loaded_from or '—'}",
        f"Full native : {is_full_native()}",
    ]
    if _lib is not None and hasattr(_lib, "antipc_native_version"):
        lines.append(f"Version C   : {_lib.antipc_native_version().decode()}")
    return "\n".join(lines)


def bench_native(limit: int = 50_000, mdc_n: int = 1_047_029) -> dict:
    """Compara tiempos Python vs C donde aplique."""
    load_native()
    results: dict = {"backend": _backend, "limit": limit, "mdc_n": mdc_n}

    t0 = time.perf_counter()
    c_count = sieve_count(limit)
    results["sieve_c_ms"] = int((time.perf_counter() - t0) * 1000)
    results["sieve_count"] = c_count

    from vma_methods_bridge import sieve_count_python

    t0 = time.perf_counter()
    py_count = sieve_count_python(limit)
    results["sieve_py_ms"] = int((time.perf_counter() - t0) * 1000)
    results["sieve_py_count"] = py_count

    payload = b"AntiPC-native-bench-" + str(limit).encode()
    t0 = time.perf_counter()
    h1 = k3_hash_buffer(payload)
    results["hash_ms"] = int((time.perf_counter() - t0) * 1000)
    results["hash_hex"] = f"{h1:08X}"

    t0 = time.perf_counter()
    f = mdc_factor(mdc_n)
    results["mdc_ms"] = int((time.perf_counter() - t0) * 1000)
    results["mdc_factor"] = f

    bits = "0101101011000010"
    tales = [3, 5, 8, 13, 21]
    puntos = [6, 12, 18]
    if is_full_native():
        t0 = time.perf_counter()
        geo = geo_converge(bits, tales, puntos)
        results["geo_ms"] = int((time.perf_counter() - t0) * 1000)
        results["geo_perimeter"] = geo

    results["criva_10k"] = criva_estimate(10_000.0)

    if hasattr(_lib, "antipc_mdc_scan_trains"):
        t0 = time.perf_counter()
        mdc_scan_trains(mdc_n if mdc_n < 10_000_000 else 1147)
        results["mdc_trains_ms"] = int((time.perf_counter() - t0) * 1000)

    if hasattr(_lib, "antipc_sieve_hibrida_count"):
        t0 = time.perf_counter()
        results["sieve_hibrida_count"] = sieve_hibrida_count(limit)
        results["sieve_hibrida_ms"] = int((time.perf_counter() - t0) * 1000)

    if hasattr(_lib, "antipc_newton_log"):
        t0 = time.perf_counter()
        nr = newton_log(121.0, familia="cuadrados")
        results["newton_ms"] = int((time.perf_counter() - t0) * 1000)
        results["newton_j"] = nr["j"]
        results["newton_ok"] = nr["converged"]

    if hasattr(_lib, "antipc_mdc_ksweep_predict"):
        t0 = time.perf_counter()
        f, ev = mdc_ksweep(100003 * 100019, predict=True)
        results["ksweep_ms"] = int((time.perf_counter() - t0) * 1000)
        results["ksweep_factor"] = f
        results["ksweep_evals"] = ev

    if hasattr(_lib, "antipc_geo_masivo_crypt"):
        msg = b"AntiPC-Aleatorovix-bench"
        t0 = time.perf_counter()
        enc, hf = geo_masivo_crypt(msg, 43210)
        dec, _ = geo_masivo_crypt(enc, 43210, decrypt=True, hash_fase=hf)
        results["aleatorovix_ms"] = int((time.perf_counter() - t0) * 1000)
        results["aleatorovix_ok"] = dec == msg

    if hasattr(_lib, "antipc_k3_heavy_hash"):
        t0 = time.perf_counter()
        results["k3_heavy"] = f"{k3_heavy_hash(payload):08X}"
        results["k3_heavy_ms"] = int((time.perf_counter() - t0) * 1000)
        results["k3_redundant"] = [f"{h:08X}" for h in k3_redundant_hashes(payload, 3)]

    return results