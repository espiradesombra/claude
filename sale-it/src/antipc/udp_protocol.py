"""Binary UDP protocol — WORK/RESULT + autenticación HMAC (HASHTOOLCODE)."""

from __future__ import annotations

import struct
from enum import IntEnum

HEADER = struct.Struct("!B I H")
RESULT = struct.Struct("!B I I")
DONE = struct.Struct("!B")
AUTH_REQ = struct.Struct("!B 32s")       # user_id utf-8 padded
AUTH_CHALLENGE = struct.Struct("!B 16s")  # nonce
AUTH_RESPONSE = struct.Struct("!B 32s")  # hmac digest
AUTH_STATUS = struct.Struct("!B b")       # ok=1 fail=0


class MsgType(IntEnum):
    WORK = 0x01
    RESULT = 0x02
    DONE = 0x03
    PING = 0x04
    PONG = 0x05
    AUTH_REQUEST = 0x06
    AUTH_CHALLENGE = 0x07
    AUTH_RESPONSE = 0x08
    AUTH_OK = 0x09
    AUTH_FAIL = 0x0A


def _pad32(text: str) -> bytes:
    b = text.encode("utf-8")[:32]
    return b + b"\x00" * (32 - len(b))


def pack_work(seq: int, payload: bytes) -> bytes:
    return HEADER.pack(MsgType.WORK, seq, len(payload)) + payload


def unpack_work(data: bytes) -> tuple[int, bytes]:
    msg_type, seq, length = HEADER.unpack_from(data)
    if msg_type != MsgType.WORK:
        raise ValueError(f"expected WORK, got {msg_type}")
    start = HEADER.size
    return seq, data[start : start + length]


def pack_result(seq: int, digest: int) -> bytes:
    return RESULT.pack(MsgType.RESULT, seq, digest & 0xFFFFFFFF)


def unpack_result(data: bytes) -> tuple[int, int]:
    msg_type, seq, digest = RESULT.unpack_from(data)
    if msg_type != MsgType.RESULT:
        raise ValueError(f"expected RESULT, got {msg_type}")
    return seq, digest


def pack_done() -> bytes:
    return DONE.pack(MsgType.DONE)


def pack_ping() -> bytes:
    return DONE.pack(MsgType.PING)


def pack_pong() -> bytes:
    return DONE.pack(MsgType.PONG)


def pack_auth_request(user_id: str) -> bytes:
    return AUTH_REQ.pack(MsgType.AUTH_REQUEST, _pad32(user_id))


def unpack_auth_request(data: bytes) -> str:
    msg_type, uid = AUTH_REQ.unpack_from(data)
    if msg_type != MsgType.AUTH_REQUEST:
        raise ValueError("expected AUTH_REQUEST")
    return uid.rstrip(b"\x00").decode("utf-8")


def pack_auth_challenge(nonce: bytes) -> bytes:
    return AUTH_CHALLENGE.pack(MsgType.AUTH_CHALLENGE, nonce[:16].ljust(16, b"\x00"))


def unpack_auth_challenge(data: bytes) -> bytes:
    msg_type, nonce = AUTH_CHALLENGE.unpack_from(data)
    if msg_type != MsgType.AUTH_CHALLENGE:
        raise ValueError("expected AUTH_CHALLENGE")
    return nonce


def pack_auth_response(digest: bytes) -> bytes:
    return AUTH_RESPONSE.pack(MsgType.AUTH_RESPONSE, digest[:32].ljust(32, b"\x00"))


def unpack_auth_response(data: bytes) -> bytes:
    msg_type, digest = AUTH_RESPONSE.unpack_from(data)
    if msg_type != MsgType.AUTH_RESPONSE:
        raise ValueError("expected AUTH_RESPONSE")
    return digest


def pack_auth_ok() -> bytes:
    return AUTH_STATUS.pack(MsgType.AUTH_OK, 1)


def pack_auth_fail() -> bytes:
    return AUTH_STATUS.pack(MsgType.AUTH_FAIL, 0)