"""K3-Toffoli — acumulación reversible por puerta; banderitas opcionales."""

from __future__ import annotations

from .k3_core import k3_bytes


def toffoli_gate(state: int, block: int, prev: int) -> int:
    return (state ^ (block & prev)) & 0xFFFFFFFF


def hash_toffoli(data: bytes, block_bits: int = 32) -> dict:
    step = max(1, block_bits // 8)
    state = 0xAAAA5555
    prev = 0x5555AAAA
    checkpoints: list[dict] = []
    for i in range(0, len(data), step):
        block_val = int.from_bytes(data[i : i + step], "big")
        nxt = toffoli_gate(state, block_val, prev)
        prev = state
        state = nxt
        if i % (step * 8) == 0:
            checkpoints.append({"i": i, "state": state})

    # mezcla final con K3 del estado
    digest = k3_bytes(state.to_bytes(4, "big") + prev.to_bytes(4, "big") + data[:32])
    return {
        "family": "toffoli_k3",
        "reversible": True,  # puerta; la cadena completa necesita data o checkpoints
        "digest": digest,
        "digest_hex": f"{digest:08x}",
        "meta": {
            "final_state": state,
            "prev_state": prev,
            "block_bits": block_bits,
            "checkpoints": checkpoints[:16],
        },
    }
