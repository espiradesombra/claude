#!/usr/bin/env python3
"""
Compresor por enlazamiento de hashes reversibles (prototipo VMA)
===============================================================
Idea (hash cinemático + metadatos):
  - El "hash" no es solo digest ciego: lleva metadatos para REPRODUCIR / VERIFICAR.
  - Compresión = deduplicar bloques que colisionan en (digest, meta-estructura)
    y encadenar H_i = K3(H_{i-1} || bloque_i) de forma verificable.

Modos:
  pack   archivo → .vhc (VMA Hash Chain)
  unpack .vhc → archivo
  verify .vhc (recalcula cadena)

No pretende batir zstd en datos aleatorios: brilla en corpus con bloques repetidos
(docs VMA, dumps, packs). Extensible a regeneración total cuando el bloque
sea generable solo con params (modo "espejo").
"""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import sys
import zipfile
from pathlib import Path

# --- K3 mínimo (32-bit), alineado con hashtool / antipc ---
OFFSETS = (5, 7, 9, 11, 13, 17, 19, 23)
MAGICO = 0x9E3779B1
FINAL_MUL_A = 0x85EBCA6B
FINAL_MUL_B = 0xC2B2AE35
SEED = 0x1F2E3D4C
N_REGS = 4


def _rotl32(v: int, n: int) -> int:
    n &= 31
    v &= 0xFFFFFFFF
    if n == 0:
        return v
    return ((v << n) | (v >> (32 - n))) & 0xFFFFFFFF


def _compress(estado: list[int], b: int) -> None:
    n = N_REGS
    b &= 0xFFFFFFFF
    x = (estado[0] ^ (b & estado[1 % n])) & 0xFFFFFFFF
    for j in range(2, n):
        x ^= _rotl32(estado[j], OFFSETS[j % 8])
    x ^= _rotl32(estado[1 % n], OFFSETS[4])
    x = (x + _rotl32(estado[0], OFFSETS[5])) & 0xFFFFFFFF
    x ^= (b * MAGICO) & 0xFFFFFFFF
    x ^= b
    x = (x + _rotl32(b, OFFSETS[6])) & 0xFFFFFFFF
    x ^= _rotl32(b, OFFSETS[7])
    x ^= _rotl32(x, OFFSETS[2])
    x = (x * MAGICO) & 0xFFFFFFFF
    x ^= x >> 15
    prev0 = estado[0]
    estado[0] = x & 0xFFFFFFFF
    temp_prev = prev0
    for i in range(1, n):
        temp_actual = estado[i]
        estado[i] = (_rotl32(estado[i], OFFSETS[0] + i) ^ temp_prev) & 0xFFFFFFFF
        temp_prev = temp_actual
        if i == n - 1:
            estado[i] = (estado[i] ^ b) & 0xFFFFFFFF


def _finalize(regs: list[int]) -> int:
    acc = regs[0]
    for i in range(1, N_REGS):
        acc ^= _rotl32(regs[i], 5 + i)
        acc &= 0xFFFFFFFF
    acc ^= acc >> 16
    acc = (acc * FINAL_MUL_A) & 0xFFFFFFFF
    acc ^= acc >> 13
    acc = (acc * FINAL_MUL_B) & 0xFFFFFFFF
    acc ^= acc >> 16
    return acc & 0xFFFFFFFF


def k3_bytes(data: bytes, seed: int = SEED) -> int:
    regs = [(seed ^ (i * MAGICO)) & 0xFFFFFFFF for i in range(N_REGS)]
    # pad to 4
    if len(data) % 4:
        data = data + b"\x00" * (4 - len(data) % 4)
    for i in range(0, len(data), 4):
        bloque = (data[i] << 24) | (data[i + 1] << 16) | (data[i + 2] << 8) | data[i + 3]
        _compress(regs, bloque)
    return _finalize(regs)


def cinematic_walk(bits: list[int]) -> dict:
    """
    Hash cinemático (versión compacta de la idea VMA):
      bit 1 → acelerar
      bit 0 → decelerar 1/3 (solo el primero de una racha), luego la racha
              de ceros consecutivos acelera al contar longitud
    Devuelve digest + metadatos REVERSIBLES (bastan para re-simular).
    """
    v = 1.0
    t = 0.0
    choques = 0
    rachas_cero: list[int] = []
    eventos: list[dict] = []
    i = 0
    while i < len(bits):
        if bits[i] == 1:
            v *= 2.0
            t += 1.0 / max(v, 1e-9)
            eventos.append({"i": i, "bit": 1, "v": v, "op": "accel"})
            i += 1
        else:
            # primer 0 de la racha: decel 1/3
            v *= 2.0 / 3.0
            # longitud de ceros consecutivos → aceleración por conteo
            j = i
            while j < len(bits) and bits[j] == 0:
                j += 1
            run = j - i
            rachas_cero.append(run)
            v *= 1.0 + 0.15 * run
            t += run / max(v, 1e-9)
            choques += 1
            eventos.append({"i": i, "bit": 0, "run": run, "v": v, "op": "zero_run"})
            i = j

    meta = {
        "pasadas": 1,
        "factor1": 2.0,
        "factor0": 2.0 / 3.0,
        "modo0": "primer_0_decel_luego_racha_acelera",
        "v_final": v,
        "t_final": t,
        "numChoques": choques,
        "rachas_cero": rachas_cero,
        "n_bits": len(bits),
    }
    # digest mezclado con K3 de la trayectoria
    traj = struct.pack("<ddI", float(v), float(t), choques & 0xFFFFFFFF)
    traj += bytes(min(b, 255) for b in rachas_cero[:64])
    digest = k3_bytes(traj)
    return {"digest": digest, "digest_hex": f"{digest:08x}", "meta": meta}


def bytes_to_bits(data: bytes) -> list[int]:
    bits: list[int] = []
    for b in data:
        for k in range(7, -1, -1):
            bits.append((b >> k) & 1)
    return bits


def pack_file(src: Path, dst: Path, block_size: int = 4096) -> dict:
    raw = src.read_bytes()
    blocks = [raw[i : i + block_size] for i in range(0, len(raw), block_size)] or [b""]

    dictionary: dict[str, bytes] = {}  # sha256 → payload
    links: list[dict] = []
    prev_h = 0
    saved = 0

    for idx, blk in enumerate(blocks):
        content_key = hashlib.sha256(blk).hexdigest()
        cin = cinematic_walk(bytes_to_bits(blk))
        # enlace de cadena reversible: H_i = K3( H_{i-1} || blk || digest_cin )
        chain_payload = prev_h.to_bytes(4, "big") + blk + cin["digest"].to_bytes(4, "big")
        h_i = k3_bytes(chain_payload)

        if content_key in dictionary:
            links.append(
                {
                    "i": idx,
                    "ref": content_key,
                    "cin": cin["digest_hex"],
                    "meta": cin["meta"],
                    "H": f"{h_i:08x}",
                    "H_prev": f"{prev_h:08x}",
                    "len": len(blk),
                }
            )
            saved += len(blk)
        else:
            dictionary[content_key] = blk
            links.append(
                {
                    "i": idx,
                    "store": content_key,
                    "cin": cin["digest_hex"],
                    "meta": cin["meta"],
                    "H": f"{h_i:08x}",
                    "H_prev": f"{prev_h:08x}",
                    "len": len(blk),
                }
            )
        prev_h = h_i

    manifest = {
        "format": "VHC1",
        "source_name": src.name,
        "source_size": len(raw),
        "block_size": block_size,
        "n_blocks": len(blocks),
        "links": links,
        "dict_keys": list(dictionary.keys()),
        "note": (
            "Reversible: meta cinemática + H_i=K3(H_prev||block||cin). "
            "Contenedor ZIP: manifest.json + dict/<sha256>."
        ),
    }
    if dst.suffix.lower() != ".vhc":
        dst = dst.with_suffix(dst.suffix + ".vhc") if dst.suffix else Path(str(dst) + ".vhc")
    with zipfile.ZipFile(dst, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("manifest.json", json.dumps(manifest, indent=2))
        for key, payload in dictionary.items():
            zf.writestr(f"dict/{key}", payload)
    packed = dst.stat().st_size
    return {
        "source_size": len(raw),
        "packed_size": packed,
        "ratio": packed / max(len(raw), 1),
        "blocks": len(blocks),
        "unique_blocks": len(dictionary),
        "bytes_deduped": saved,
        "out": str(dst),
    }


def _load_vhc(src: Path) -> tuple[dict, dict[str, bytes]]:
    with zipfile.ZipFile(src, "r") as zf:
        container = json.loads(zf.read("manifest.json").decode("utf-8"))
        dictionary: dict[str, bytes] = {}
        for name in zf.namelist():
            if name.startswith("dict/") and not name.endswith("/"):
                dictionary[name.split("/", 1)[1]] = zf.read(name)
    assert container.get("format") == "VHC1"
    return container, dictionary


def unpack_file(src: Path, dst: Path) -> dict:
    container, dictionary = _load_vhc(src)
    out = bytearray()
    prev_h = 0
    for link in container["links"]:
        key = link.get("ref") or link.get("store")
        blk = dictionary[key]
        if len(blk) != link["len"]:
            raise ValueError(f"len mismatch block {link['i']}")
        cin = cinematic_walk(bytes_to_bits(blk))
        if cin["digest_hex"] != link["cin"]:
            raise ValueError(f"cinematic digest mismatch block {link['i']}")
        chain_payload = prev_h.to_bytes(4, "big") + blk + cin["digest"].to_bytes(4, "big")
        h_i = k3_bytes(chain_payload)
        if f"{h_i:08x}" != link["H"]:
            raise ValueError(f"chain H mismatch block {link['i']}")
        out.extend(blk)
        prev_h = h_i
    dst.write_bytes(bytes(out))
    if len(out) != container["source_size"]:
        raise ValueError("size mismatch after unpack")
    return {"ok": True, "size": len(out), "out": str(dst)}


def verify_file(src: Path) -> dict:
    tmp = src.with_suffix(".vhc.verify.bin")
    try:
        unpack_file(src, tmp)
        return {"ok": True, "message": "cadena + digests cinemáticos OK"}
    finally:
        if tmp.exists():
            tmp.unlink()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Compresor VHC — enlaces de hash reversibles")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_pack = sub.add_parser("pack", help="comprimir → .vhc (zip)")
    p_pack.add_argument("src")
    p_pack.add_argument("-o", "--out", default=None)
    p_pack.add_argument("-b", "--block-size", type=int, default=4096)

    p_un = sub.add_parser("unpack", help="descomprimir")
    p_un.add_argument("src")
    p_un.add_argument("-o", "--out", required=True)

    p_ver = sub.add_parser("verify", help="verificar cadena")
    p_ver.add_argument("src")

    args = ap.parse_args(argv)
    if args.cmd == "pack":
        src = Path(args.src)
        out = Path(args.out) if args.out else src.with_suffix(src.suffix + ".vhc")
        stats = pack_file(src, out, block_size=args.block_size)
        print(json.dumps(stats, indent=2))
        return 0
    if args.cmd == "unpack":
        print(json.dumps(unpack_file(Path(args.src), Path(args.out)), indent=2))
        return 0
    if args.cmd == "verify":
        print(json.dumps(verify_file(Path(args.src)), indent=2))
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
