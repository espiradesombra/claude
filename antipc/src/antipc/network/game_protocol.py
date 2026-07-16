"""
Protocolo UDP fijo 64B para MMO — encaja en SlotRing (recvinto directo).

Paquete de estado de jugador: posición, velocidad, tick, shard, fingerprint K3.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from enum import IntEnum

MAGIC = 0x4D4D  # "MM" — Massive Multiplayer AntiPC
VERSION = 1
GAME_PKT_SIZE = 64

_GAME_STRUCT = struct.Struct("!HBBHIIfff fffIIHH14s")


class GameMsgType(IntEnum):
    MOVE = 1
    ACTION = 2
    HEARTBEAT = 3
    SPAWN = 4


@dataclass(frozen=True)
class GameStatePacket:
    entity_id: int
    tick: int
    shard_id: int
    px: float
    py: float
    pz: float
    vx: float
    vy: float
    vz: float
    seq: int
    k3_tag: int
    yaw: int
    flags: int
    msg_type: int = GameMsgType.MOVE


def _k3_tag_for(fields_without_tag: tuple) -> int:
    body = struct.pack("!HBBHIIfff fffIIHH14s", *fields_without_tag, b"\x00" * 14)
    try:
        from native_engine import k3_hash_buffer

        return k3_hash_buffer(body) & 0xFFFFFFFF
    except Exception:
        from hash_engine import k3_hash_buffer

        return k3_hash_buffer(body) & 0xFFFFFFFF


def pack_state(
    *,
    entity_id: int,
    tick: int,
    shard_id: int = 0,
    px: float = 0.0,
    py: float = 0.0,
    pz: float = 0.0,
    vx: float = 0.0,
    vy: float = 0.0,
    vz: float = 0.0,
    seq: int = 0,
    yaw: int = 0,
    flags: int = 0,
    msg_type: int = GameMsgType.MOVE,
    k3_tag: int | None = None,
) -> bytes:
    """Serializa exactamente 64 bytes listos para recvinto(slot)."""
    base = (
        MAGIC,
        VERSION,
        int(msg_type) & 0xFF,
        int(shard_id) & 0xFFFF,
        int(entity_id) & 0xFFFFFFFF,
        int(tick) & 0xFFFFFFFF,
        float(px),
        float(py),
        float(pz),
        float(vx),
        float(vy),
        float(vz),
        int(seq) & 0xFFFFFFFF,
        0,
        int(yaw) & 0xFFFF,
        int(flags) & 0xFFFF,
    )
    if k3_tag is None:
        k3_tag = _k3_tag_for(base)
    return struct.pack(
        "!HBBHIIfff fffIIHH14s",
        MAGIC,
        VERSION,
        int(msg_type) & 0xFF,
        int(shard_id) & 0xFFFF,
        int(entity_id) & 0xFFFFFFFF,
        int(tick) & 0xFFFFFFFF,
        float(px),
        float(py),
        float(pz),
        float(vx),
        float(vy),
        float(vz),
        int(seq) & 0xFFFFFFFF,
        int(k3_tag) & 0xFFFFFFFF,
        int(yaw) & 0xFFFF,
        int(flags) & 0xFFFF,
        b"\x00" * 14,
    )


def unpack_state(data: bytes | memoryview) -> GameStatePacket | None:
    raw = bytes(data[:GAME_PKT_SIZE])
    if len(raw) < _GAME_STRUCT.size:
        return None
    (
        magic,
        ver,
        msg_type,
        shard_id,
        entity_id,
        tick,
        px,
        py,
        pz,
        vx,
        vy,
        vz,
        seq,
        k3_tag,
        yaw,
        flags,
        _pad,
    ) = _GAME_STRUCT.unpack(raw)
    if magic != MAGIC or ver != VERSION:
        return None
    return GameStatePacket(
        entity_id=entity_id,
        tick=tick,
        shard_id=shard_id,
        px=px,
        py=py,
        pz=pz,
        vx=vx,
        vy=vy,
        vz=vz,
        seq=seq,
        k3_tag=k3_tag,
        yaw=yaw,
        flags=flags,
        msg_type=msg_type,
    )


def verify_k3_tag(data: bytes | memoryview) -> bool:
    pkt = unpack_state(data)
    if pkt is None:
        return False
    expected = pack_state(
        entity_id=pkt.entity_id,
        tick=pkt.tick,
        shard_id=pkt.shard_id,
        px=pkt.px,
        py=pkt.py,
        pz=pkt.pz,
        vx=pkt.vx,
        vy=pkt.vy,
        vz=pkt.vz,
        seq=pkt.seq,
        yaw=pkt.yaw,
        flags=pkt.flags,
        msg_type=pkt.msg_type,
    )
    return bytes(data[:GAME_PKT_SIZE]) == expected