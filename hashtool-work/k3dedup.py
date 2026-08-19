#!/usr/bin/env python3
"""
k3dedup.py — portable: detecta duplicados exactos (tamaño + hash K3).

Igual idea que el ejemplo C `k3dedup.c` del pack HASHTOOL / k3hash.
No borra nada: solo lista grupos. Verifica a mano antes de eliminar.

Uso (Windows PowerShell):
  python k3dedup.py C:\\ruta\\a\\carpeta
  python k3dedup.py C:\\ruta --ext .py,.docx
  Get-ChildItem -Recurse -File C:\\ruta | ForEach-Object FullName | python k3dedup.py -

Uso (Linux/macOS):
  python3 k3dedup.py /ruta
  find /ruta -type f | python3 k3dedup.py -
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

# --- K3 mínimo (mismo algoritmo que antipc/k3_hash.py / k3hash.c) ---
OFFSETS = (5, 7, 9, 11, 13, 17, 19, 23)
MAGICO = 0x9E3779B1
FINAL_MUL_A = 0x85EBCA6B
FINAL_MUL_B = 0xC2B2AE35
SEED = 0x1F2E3D4C
N_REGS = 4
ANCHO = 4  # 32-bit


def _rotl32(value: int, positions: int) -> int:
    positions &= 31
    value &= 0xFFFFFFFF
    if positions == 0:
        return value
    return ((value << positions) | (value >> (32 - positions))) & 0xFFFFFFFF


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


def k3_hash_file(path: Path, chunk: int = 1 << 20) -> int:
    regs = [(SEED ^ (i * MAGICO)) & 0xFFFFFFFF for i in range(N_REGS)]
    partial = bytearray()
    with path.open("rb") as f:
        while True:
            data = f.read(chunk)
            if not data:
                break
            buf = partial + data
            i = 0
            while i + ANCHO <= len(buf):
                bloque = 0
                for b in range(ANCHO):
                    bloque = (bloque << 8) | buf[i + b]
                _compress(regs, bloque)
                i += ANCHO
            partial = bytearray(buf[i:])
    if partial:
        bloque = 0
        for b in range(ANCHO):
            byte = partial[b] if b < len(partial) else 0
            bloque = (bloque << 8) | byte
        _compress(regs, bloque)
    return _finalize(regs)


def collect_paths(root: Path, exts: set[str] | None) -> list[Path]:
    out: list[Path] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        # skip junk / vcs
        base = Path(dirpath)
        if any(p in {".git", "node_modules", "__pycache__", ".venv"} for p in base.parts):
            continue
        for name in filenames:
            p = base / name
            if exts and p.suffix.lower() not in exts:
                continue
            out.append(p)
    return out


def paths_from_stdin() -> list[Path]:
    out: list[Path] = []
    for line in sys.stdin:
        s = line.strip().strip('"')
        if s:
            out.append(Path(s))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Detecta duplicados exactos con hash K3 (no borra nada)."
    )
    ap.add_argument(
        "ruta",
        help="Carpeta a escanear, o '-' para leer rutas por stdin",
    )
    ap.add_argument(
        "--ext",
        default="",
        help="Filtro de extensiones separadas por coma, ej: .py,.docx,.pdf",
    )
    ap.add_argument(
        "--delete-report",
        action="store_true",
        help="Escribe sugerencias de borrado (conserva 1ª ruta de cada grupo) a stderr",
    )
    args = ap.parse_args()

    exts = {e.strip().lower() if e.strip().startswith(".") else f".{e.strip().lower()}"
            for e in args.ext.split(",") if e.strip()} or None

    if args.ruta == "-":
        files = paths_from_stdin()
    else:
        root = Path(args.ruta)
        if not root.is_dir():
            print(f"[error] no es carpeta: {root}", file=sys.stderr)
            return 2
        files = collect_paths(root, exts)

    # Agrupa primero por tamaño (rápido), luego hash solo si hay colisión de tamaño
    by_size: dict[int, list[Path]] = defaultdict(list)
    skipped = 0
    for p in files:
        try:
            by_size[p.stat().st_size].append(p)
        except OSError as e:
            print(f"[aviso] no se pudo stat: {p} ({e})", file=sys.stderr)
            skipped += 1

    groups_out: list[tuple[int, int, list[Path]]] = []
    hashed = 0
    for size, paths in by_size.items():
        if len(paths) < 2:
            continue
        by_hash: dict[int, list[Path]] = defaultdict(list)
        for p in paths:
            try:
                h = k3_hash_file(p)
                hashed += 1
                by_hash[h].append(p)
            except OSError as e:
                print(f"[aviso] no se pudo leer: {p} ({e})", file=sys.stderr)
                skipped += 1
        for h, same in by_hash.items():
            if len(same) > 1:
                groups_out.append((size, h, same))

    groups_out.sort(key=lambda g: (-g[0] * (len(g[2]) - 1), -g[0]))

    recoverable = 0
    for i, (size, h, paths) in enumerate(groups_out, 1):
        recoverable += size * (len(paths) - 1)
        print(f"=== Grupo {i} (tam={size} bytes, hash=0x{h:08X}, {len(paths)} copias) ===")
        for p in paths:
            print(f"  {p}")
        if args.delete_report:
            keep = paths[0]
            for p in paths[1:]:
                print(f"# sugerido borrar (mantener {keep}): {p}", file=sys.stderr)

    print(
        f"\n[k3dedup] {len(files)} rutas, {hashed} hasheados (solo tamaños repetidos), "
        f"{len(groups_out)} grupos, ~{recoverable} bytes recuperables, avisos={skipped}.",
        file=sys.stderr,
    )
    print(
        "[k3dedup] Coincidencia tamaño+hash 32-bit. Verifica byte a byte antes de borrar críticos.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
