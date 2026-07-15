#!/usr/bin/env python3
"""Regenera MANIFIESTO_HASHES.json para el paquete 33x1."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

DESKTOP = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent

# Rutas relativas al Desktop (ajusta si mueves carpetas)
CANDIDATES = [
    "repo/vma-methods/vma_methods/cribas.py",
    "repo/vma-methods/vma_methods/criva.py",
    "repo/vma-methods/vma_methods/newton.py",
    "repo/graficas y explicaciones/00_INDICE.txt",
    "repo/teoremas/INDICE_MAESTRO.txt",
    "repo/teoremas/29_metodo_pitagorico_visual_divisores.txt",
    "grok/antipc/dist/k3hash.dll",
    "grok/antipc/src/antipc/mdc_lib/mdc_analyze.py",
    "repo/graficas y explicaciones/13_ZypyZape_A_vs_B.png",
    "33x1/02_USO_CIVIL.txt",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    files = []
    for rel in CANDIDATES:
        p = DESKTOP / rel.replace("/", "\\")
        if not p.is_file():
            files.append({"path": rel, "status": "missing"})
            continue
        files.append({
            "path": rel,
            "bytes": p.stat().st_size,
            "sha256": sha256_file(p),
            "status": "ok",
        })
    ok = sum(1 for f in files if f.get("status") == "ok")
    payload = {
        "generated": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        "author": "Víctor Manzanares Alberola",
        "plan": "33x1",
        "files_ok": ok,
        "files_total": len(files),
        "files": files,
    }
    out_path = OUT / "MANIFIESTO_HASHES.json"
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Manifiesto: {out_path}")
    print(f"Firmados: {ok}/{len(files)}")


if __name__ == "__main__":
    main()