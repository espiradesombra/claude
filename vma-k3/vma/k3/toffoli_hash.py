"""
Theorem K3 — variante Toffoli + signal detector (banderitas intermedias).

Hash streaming con acumulación reversible y puntos de guardado.
NO criptográfico — integridad / telemetría industrial.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List


def toffoli_gate(state: int, block: int, prev: int) -> int:
    return (state ^ (block & prev)) & 0xFFFFFFFF


@dataclass
class Banderita:
    checksum: int
    signal_phase: int


@dataclass
class K3ToffoliResult:
    final_hash: int
    remaining_bytes: int
    timestamp: float
    puntos_control: Dict[int, Banderita] = field(default_factory=dict)


def theorem_k3_hash(
    data: bytes,
    block_bits: int = 32,
    banderitas: List[int] | None = None,
    anchor_state: int = 0xAAAA5555,
    anchor_prev: int = 0x5555AAAA,
) -> K3ToffoliResult:
    banderitas = banderitas or []
    byte_step = max(1, block_bits // 8)
    hash_state = anchor_state & 0xFFFFFFFF
    prev_state = anchor_prev & 0xFFFFFFFF
    puntos: Dict[int, Banderita] = {}
    ts_base = int(time.time())

    for i in range(0, len(data), byte_step):
        block_val = int.from_bytes(data[i : i + byte_step], "big")
        next_state = toffoli_gate(hash_state, block_val, prev_state)
        prev_state = hash_state
        hash_state = next_state & 0xFFFFFFFF

        if i in banderitas:
            signal = (hash_state ^ ts_base) & 0xFF
            puntos[i] = Banderita(hash_state, signal)

    return K3ToffoliResult(
        final_hash=hash_state,
        remaining_bytes=len(data) % byte_step,
        timestamp=time.time(),
        puntos_control=puntos,
    )