#!/usr/bin/env python3
"""
Compresor VHC2 — enlaces + ANIDACIÓN de hashes reversibles
==========================================================
Capas (pipeline por defecto):
  cinematic → espejo → toffoli → k3

Cada capa hashea el blob (digest+meta) de la anterior.
Cadena de bloques: H_i = K3(H_{i-1} || bloque || outer_digest)

Ver MAPA_HASHES_REVERSIBLES.md
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

from hashes import FAMILIES, nest_pipeline
from hashes.k3_core import k3_bytes
from hashes.nest import DEFAULT_PIPELINE


def pack_file(
    src: Path,
    dst: Path,
    block_size: int = 4096,
    pipeline: list[str] | None = None,
) -> dict:
    pipeline = pipeline or list(DEFAULT_PIPELINE)
    raw = src.read_bytes()
    blocks = [raw[i : i + block_size] for i in range(0, len(raw), block_size)] or [b""]

    dictionary: dict[str, bytes] = {}
    links: list[dict] = []
    prev_h = 0
    saved = 0

    for idx, blk in enumerate(blocks):
        content_key = hashlib.sha256(blk).hexdigest()
        nested = nest_pipeline(blk, pipeline)
        outer = nested["outer_digest"]
        chain_payload = prev_h.to_bytes(4, "big") + blk + outer.to_bytes(4, "big")
        h_i = k3_bytes(chain_payload)

        link = {
            "i": idx,
            "pipeline": pipeline,
            "layers": nested["layers"],
            "outer": nested["outer_digest_hex"],
            "H": f"{h_i:08x}",
            "H_prev": f"{prev_h:08x}",
            "len": len(blk),
        }
        if content_key in dictionary:
            link["ref"] = content_key
            saved += len(blk)
        else:
            dictionary[content_key] = blk
            link["store"] = content_key
        links.append(link)
        prev_h = h_i

    manifest = {
        "format": "VHC2",
        "source_name": src.name,
        "source_size": len(raw),
        "block_size": block_size,
        "n_blocks": len(blocks),
        "pipeline": pipeline,
        "links": links,
        "dict_keys": list(dictionary.keys()),
        "families_available": list(FAMILIES.keys()),
        "note": "VHC2: capas anidadas reversibles + dedup + cadena K3 entre bloques.",
    }
    if not str(dst).endswith(".vhc"):
        dst = Path(str(dst) + ("" if str(dst).endswith(".vhc") else ".vhc"))
    with zipfile.ZipFile(dst, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("manifest.json", json.dumps(manifest, indent=2))
        for key, payload in dictionary.items():
            zf.writestr(f"dict/{key}", payload)

    packed = dst.stat().st_size
    return {
        "format": "VHC2",
        "pipeline": pipeline,
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
        dictionary = {
            name.split("/", 1)[1]: zf.read(name)
            for name in zf.namelist()
            if name.startswith("dict/") and not name.endswith("/")
        }
    return container, dictionary


def unpack_file(src: Path, dst: Path) -> dict:
    container, dictionary = _load_vhc(src)
    fmt = container.get("format")
    if fmt not in ("VHC1", "VHC2"):
        raise ValueError(f"formato desconocido: {fmt}")
    out = bytearray()
    prev_h = 0
    for link in container["links"]:
        key = link.get("ref") or link.get("store")
        blk = dictionary[key]
        if len(blk) != link["len"]:
            raise ValueError(f"len mismatch block {link['i']}")

        if fmt == "VHC2":
            pipe = link.get("pipeline") or container.get("pipeline") or DEFAULT_PIPELINE
            nested = nest_pipeline(blk, pipe)
            if nested["outer_digest_hex"] != link["outer"]:
                raise ValueError(f"outer nest digest mismatch block {link['i']}")
            # verificar cada capa
            for a, b in zip(nested["layers"], link["layers"]):
                if a["family"] != b["family"] or a["digest_hex"] != b["digest_hex"]:
                    raise ValueError(f"layer mismatch block {link['i']}: {a['family']}")
            outer = nested["outer_digest"]
        else:
            # VHC1 legacy: solo cinematic outer in 'cin'
            from hashes.cinematic import hash_cinematic

            cin = hash_cinematic(blk)
            if cin["digest_hex"] != link.get("cin"):
                raise ValueError(f"cinematic mismatch block {link['i']}")
            outer = cin["digest"]

        chain_payload = prev_h.to_bytes(4, "big") + blk + outer.to_bytes(4, "big")
        h_i = k3_bytes(chain_payload)
        if f"{h_i:08x}" != link["H"]:
            raise ValueError(f"chain H mismatch block {link['i']}")
        out.extend(blk)
        prev_h = h_i

    dst.write_bytes(bytes(out))
    if len(out) != container["source_size"]:
        raise ValueError("size mismatch")
    return {"ok": True, "size": len(out), "format": fmt, "out": str(dst)}


def verify_file(src: Path) -> dict:
    tmp = src.with_suffix(".verify.bin")
    try:
        r = unpack_file(src, tmp)
        return {"ok": True, "format": r["format"], "message": "capas anidadas + cadena OK"}
    finally:
        if tmp.exists():
            tmp.unlink()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Compresor VHC2 — hashes reversibles anidados")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_pack = sub.add_parser("pack")
    p_pack.add_argument("src")
    p_pack.add_argument("-o", "--out", default=None)
    p_pack.add_argument("-b", "--block-size", type=int, default=4096)
    p_pack.add_argument(
        "-p",
        "--pipeline",
        default=",".join(DEFAULT_PIPELINE),
        help=f"familias separadas por coma; disponibles: {','.join(FAMILIES)}",
    )

    p_un = sub.add_parser("unpack")
    p_un.add_argument("src")
    p_un.add_argument("-o", "--out", required=True)

    p_ver = sub.add_parser("verify")
    p_ver.add_argument("src")

    p_list = sub.add_parser("families", help="listar familias de hash")

    args = ap.parse_args(argv)
    if args.cmd == "families":
        print(json.dumps({"families": list(FAMILIES), "default": DEFAULT_PIPELINE}, indent=2))
        return 0
    if args.cmd == "pack":
        src = Path(args.src)
        out = Path(args.out) if args.out else Path(str(src) + ".vhc")
        pipe = [x.strip() for x in args.pipeline.split(",") if x.strip()]
        print(json.dumps(pack_file(src, out, args.block_size, pipe), indent=2))
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
