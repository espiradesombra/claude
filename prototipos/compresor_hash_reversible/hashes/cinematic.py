"""Hash cinemático 1/0 — reversible con metadatos de trayectoria."""

from __future__ import annotations

import struct

from .k3_core import k3_bytes


def hash_cinematic(data: bytes) -> dict:
    bits: list[int] = []
    for b in data:
        for k in range(7, -1, -1):
            bits.append((b >> k) & 1)

    v = 1.0
    t = 0.0
    choques = 0
    rachas_cero: list[int] = []
    i = 0
    while i < len(bits):
        if bits[i] == 1:
            v *= 2.0
            t += 1.0 / max(v, 1e-9)
            i += 1
        else:
            v *= 2.0 / 3.0
            j = i
            while j < len(bits) and bits[j] == 0:
                j += 1
            run = j - i
            rachas_cero.append(run)
            v *= 1.0 + 0.15 * run
            t += run / max(v, 1e-9)
            choques += 1
            i = j

    meta = {
        "pasadas": 1,
        "factor1": 2.0,
        "factor0": 2.0 / 3.0,
        "modo0": "primer_0_decel_luego_racha_acelera",
        "v_final": v,
        "t_final": t,
        "numChoques": choques,
        "rachas_cero": rachas_cero[:128],
        "n_bits": len(bits),
    }
    traj = struct.pack("<ddI", float(v), float(t), choques & 0xFFFFFFFF)
    traj += bytes(min(b, 255) for b in rachas_cero[:64])
    digest = k3_bytes(traj)
    return {
        "family": "cinematic",
        "reversible": True,
        "digest": digest,
        "digest_hex": f"{digest:08x}",
        "meta": meta,
    }
