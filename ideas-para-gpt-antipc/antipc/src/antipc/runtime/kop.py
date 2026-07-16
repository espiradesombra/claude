"""
Knowledge Operations (KOP) — formato binario K3 MicroKernel.

Segun gptcomputing.txt (~9728): directorio con offsets, acceso directo sin
recorrer el objeto. Ahorra tiempo vs JSON/texto en lookup de Identity,
Signature, Trust, etc.

Layout:
  [Header 16B][KOP blobs...][Directory N*20B][Footer 4B]
"""

from __future__ import annotations

import hashlib
import struct
import time
from dataclasses import dataclass, field
from enum import IntEnum, IntFlag
from typing import Callable

MAGIC = b"K3MK"
FOOTER = b"K3ND"
FILE_VERSION = 1

FILE_HEADER = struct.Struct("<4s I H H I")  # magic, ver, kop_count, reserved, dir_offset
DIR_ENTRY = struct.Struct("<H H I I Q")     # id, version, offset, size, flags
KOP_HEADER = struct.Struct("<H H I Q")      # id, version, size, flags


class KOPId(IntEnum):
    IDENTITY = 1
    DNA = 2
    SIGNATURE = 3
    TRUST = 4
    HISTORY = 5
    COMPRESSION = 6
    ENCRYPTION = 7
    NETWORK = 8
    HASH = 9


class KOPKind(IntEnum):
    DATA = 0
    BEHAVIOR = 1


class KOPFlags(IntFlag):
    NONE = 0
    IMMUTABLE = 1 << 0
    CACHEABLE = 1 << 1
    DISTRIBUTABLE = 1 << 2
    DETERMINISTIC = 1 << 3
    STREAMABLE = 1 << 4


# Metadatos por KOP (definiciones del txt)
KOP_DEFINITIONS: dict[KOPId, dict] = {
    KOPId.IDENTITY: {
        "name": "Identity",
        "kind": KOPKind.DATA,
        "flags": KOPFlags.IMMUTABLE | KOPFlags.CACHEABLE,
        "desc": "Identificador estable del objeto de conocimiento",
    },
    KOPId.DNA: {
        "name": "DNA",
        "kind": KOPKind.DATA,
        "flags": KOPFlags.IMMUTABLE,
        "desc": "Producer, policy, parents — como nacio el conocimiento",
    },
    KOPId.SIGNATURE: {
        "name": "Signature",
        "kind": KOPKind.DATA,
        "flags": KOPFlags.IMMUTABLE | KOPFlags.CACHEABLE | KOPFlags.DETERMINISTIC,
        "desc": "S = f(op, inputs) — clave de lookup en KnowledgeBuffer",
    },
    KOPId.TRUST: {
        "name": "Trust",
        "kind": KOPKind.DATA,
        "flags": KOPFlags.CACHEABLE,
        "desc": "Consenso evaluadores (Axioma 17)",
    },
    KOPId.HISTORY: {
        "name": "History",
        "kind": KOPKind.DATA,
        "flags": KOPFlags.IMMUTABLE,
        "desc": "Append-only ledger pointer",
    },
    KOPId.COMPRESSION: {
        "name": "Compression",
        "kind": KOPKind.BEHAVIOR,
        "flags": KOPFlags.DETERMINISTIC,
        "desc": "Transforma payload; genera nueva version",
    },
    KOPId.ENCRYPTION: {
        "name": "Encryption",
        "kind": KOPKind.BEHAVIOR,
        "flags": KOPFlags.DETERMINISTIC,
        "desc": "Cifra; genera nueva version",
    },
    KOPId.NETWORK: {
        "name": "Network",
        "kind": KOPKind.BEHAVIOR,
        "flags": KOPFlags.DISTRIBUTABLE | KOPFlags.STREAMABLE,
        "desc": "Replicacion / exchange en fabric UDP",
    },
    KOPId.HASH: {
        "name": "Hash",
        "kind": KOPKind.BEHAVIOR,
        "flags": KOPFlags.CACHEABLE | KOPFlags.DETERMINISTIC | KOPFlags.STREAMABLE,
        "desc": "K3 hash — transformacion cacheable con reuse en KnowledgeBuffer",
    },
}


@dataclass(frozen=True)
class KOPBlock:
    kop_id: KOPId
    version: int
    flags: int
    data: bytes

    @property
    def header_bytes(self) -> bytes:
        return KOP_HEADER.pack(
            int(self.kop_id), self.version, len(self.data), int(self.flags)
        )

    def pack(self) -> bytes:
        return self.header_bytes + self.data


@dataclass
class K3MicroKernel:
    """Objeto binario con directorio — acceso O(1) por KOP id."""

    blocks: dict[KOPId, KOPBlock] = field(default_factory=dict)
    _directory_cache: dict[KOPId, tuple[int, int]] = field(
        default_factory=dict, repr=False
    )
    _raw: bytes | None = field(default=None, repr=False)

    def set_kop(
        self,
        kop_id: KOPId,
        data: bytes,
        version: int = 1,
        flags: int | None = None,
    ) -> None:
        if flags is None:
            flags = int(KOP_DEFINITIONS[kop_id]["flags"])
        self.blocks[kop_id] = KOPBlock(kop_id, version, flags, data)
        self._raw = None

    def get_kop(self, kop_id: KOPId) -> bytes | None:
        block = self.blocks.get(kop_id)
        if block is not None:
            return block.data
        if self._raw is not None and kop_id in self._directory_cache:
            off, size = self._directory_cache[kop_id]
            start = off + KOP_HEADER.size
            return self._raw[start : start + size]
        return None

    def get_signature_hex(self) -> str | None:
        raw = self.get_kop(KOPId.SIGNATURE)
        if raw is None:
            return None
        return raw.decode("ascii", errors="replace").strip("\x00")

    def pack(self) -> bytes:
        if not self.blocks:
            raise ValueError("K3MicroKernel vacio")

        parts: list[bytes] = []
        dir_entries: list[tuple[KOPId, int, int, KOPBlock]] = []
        offset = FILE_HEADER.size

        for kop_id in sorted(self.blocks.keys(), key=int):
            block = self.blocks[kop_id]
            blob = block.pack()
            dir_entries.append((kop_id, offset, len(block.data), block))
            parts.append(blob)
            offset += len(blob)

        dir_offset = offset
        directory = bytearray()
        for kop_id, off, data_len, block in dir_entries:
            directory.extend(
                DIR_ENTRY.pack(
                    int(kop_id),
                    block.version,
                    off,
                    data_len,
                    int(block.flags),
                )
            )

        header = FILE_HEADER.pack(
            MAGIC,
            FILE_VERSION,
            len(dir_entries),
            0,
            dir_offset,
        )
        footer = FOOTER
        raw = header + b"".join(parts) + bytes(directory) + footer
        self._raw = raw
        self._directory_cache = {
            kop_id: (off, size) for kop_id, off, size, _ in dir_entries
        }
        return raw

    @classmethod
    def unpack(cls, raw: bytes) -> K3MicroKernel:
        if len(raw) < FILE_HEADER.size + len(FOOTER):
            raise ValueError("blob K3 demasiado corto")
        magic, ver, count, _reserved, dir_offset = FILE_HEADER.unpack_from(raw)
        if magic != MAGIC:
            raise ValueError(f"magic invalido: {magic!r}")
        if raw[-4:] != FOOTER:
            raise ValueError("footer K3ND no encontrado")

        mk = cls()
        mk._raw = raw
        for i in range(count):
            base = dir_offset + i * DIR_ENTRY.size
            kid, kver, off, size, flags = DIR_ENTRY.unpack_from(raw, base)
            kop_id = KOPId(kid)
            data_start = off + KOP_HEADER.size
            data = raw[data_start : data_start + size]
            mk.blocks[kop_id] = KOPBlock(kop_id, kver, flags, data)
            mk._directory_cache[kop_id] = (off, size)
        return mk


def build_knowledge_blob(
    ref_id: str,
    signature_hex: str,
    payload: bytes,
    producer: str = "antipc",
    trust: float = 1.0,
) -> bytes:
    """Construye K3MK con KOP 001-005 + payload en DNA extension."""
    mk = K3MicroKernel()
    mk.set_kop(KOPId.IDENTITY, ref_id.encode("utf-8")[:64])
    dna = f"producer={producer}|parents=0".encode("utf-8")
    mk.set_kop(KOPId.DNA, dna)
    sig = signature_hex.encode("ascii")[:32]
    mk.set_kop(KOPId.SIGNATURE, sig.ljust(32, b"\x00")[:32])
    trust_b = struct.pack("<d", max(0.0, min(1.0, trust)))
    mk.set_kop(KOPId.TRUST, trust_b)
    mk.set_kop(KOPId.HISTORY, b"v1")
    if payload:
        mk.set_kop(KOPId.COMPRESSION, payload, flags=int(KOPFlags.STREAMABLE))
    return mk.pack()


def lookup_signature_from_blob(raw: bytes) -> str | None:
    """Acceso directo via directorio — sin parsear JSON."""
    if raw[:4] != MAGIC:
        return None
    _magic, _ver, count, _res, dir_offset = FILE_HEADER.unpack_from(raw)
    for i in range(count):
        base = dir_offset + i * DIR_ENTRY.size
        kid, _kver, off, size, _flags = DIR_ENTRY.unpack_from(raw, base)
        if kid == int(KOPId.SIGNATURE):
            start = off + KOP_HEADER.size
            return raw[start : start + size].decode("ascii", errors="replace").strip("\x00")
    return None


def benchmark_lookup(json_fn: Callable[[], str], blob: bytes, n: int = 50000) -> dict:
    """Compara lookup JSON simulado vs salto directo KOP 003."""
    t0 = time.perf_counter()
    for _ in range(n):
        json_fn()
    t_json = time.perf_counter() - t0

    t0 = time.perf_counter()
    for _ in range(n):
        lookup_signature_from_blob(blob)
    t_bin = time.perf_counter() - t0

    return {
        "iterations": n,
        "json_sec": t_json,
        "binary_sec": t_bin,
        "speedup_x": t_json / t_bin if t_bin > 0 else 0,
        "saved_pct": (1 - t_bin / t_json) * 100 if t_json > 0 else 0,
    }