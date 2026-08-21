"""Hash espejo degradante — R = R0·exp(-k·I·t·cosθ); reversible con params."""

from __future__ import annotations

import math

from .k3_core import k3_bytes


def hash_espejo(
    data: bytes,
    *,
    R0: float = 0.95,
    k: float = 0.12,
    intensidad: float = 1.0,
    angulo_deg: float = 30.0,
    tiempo: float = 1.0,
    clave: bytes = b"VMA-ESPEJO",
) -> dict:
    seed = k3_bytes(data) / 0xFFFFFFFF
    R0_d = R0 * (0.9 + seed * 0.2)
    k_d = k * (0.5 + seed * 1.0)
    I_d = intensidad * (0.5 + seed * 1.0)
    th_d = angulo_deg * (0.5 + seed * 1.0)
    t_d = tiempo * (0.5 + seed * 1.0)
    cos_th = math.cos(th_d * math.pi / 180.0)
    R = R0_d * math.exp(-k_d * I_d * t_d * cos_th)

    blob = (
        f"R={R:.12f}|R0={R0_d:.12f}|k={k_d:.12f}|I={I_d:.12f}|"
        f"th={th_d:.12f}|t={t_d:.12f}|clave={clave!r}|len={len(data)}"
    ).encode("utf-8")
    digest = k3_bytes(blob + data[:64])
    return {
        "family": "espejo",
        "reversible": True,
        "digest": digest,
        "digest_hex": f"{digest:08x}",
        "meta": {
            "R0": R0_d,
            "k": k_d,
            "intensidad": I_d,
            "angulo": th_d,
            "tiempo": t_d,
            "reflectividad": R,
            "clave": clave.decode("utf-8", errors="replace"),
        },
    }
