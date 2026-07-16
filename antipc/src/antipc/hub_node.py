#!/usr/bin/env python3
"""AntiPC UDP hub — lock-free + K3 + autenticación HMAC + lotes WORK/RESULT."""

from __future__ import annotations

import argparse
import os
import socket
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from architecture.network_auth import NetworkAuthGate
from network.network_compute import hub_drain_ring_burst, hub_process_work, tune_socket
from ring_buffer import SPSCRingBuffer
from udp_protocol import (
    MsgType,
    pack_auth_challenge,
    pack_auth_fail,
    pack_auth_ok,
    pack_pong,
    unpack_auth_request,
    unpack_auth_response,
    unpack_work,
    HEADER,
    DONE,
    AUTH_CHALLENGE,
)


def run_hub(
    port: int,
    master_addr: tuple[str, int],
    auth: NetworkAuthGate | None = None,
    require_auth: bool = False,
    verbose: bool = False,
) -> None:
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    tune_socket(sock)
    sock.bind(("0.0.0.0", port))
    sock.settimeout(0.001)

    ring: SPSCRingBuffer[tuple[int, bytes, tuple[str, int]]] = SPSCRingBuffer(8192)
    cache: dict[int, int] = {}
    cache_hits = [0]
    authenticated: set[tuple[str, int]] = set()
    pending_auth: dict[tuple[str, int], str] = {}
    processed = 0
    running = True

    if verbose:
        print(f"[hub:{port}] auth={'ON' if require_auth else 'OFF'} batch=ON", flush=True)

    while running:
        try:
            data, addr = sock.recvfrom(65535)
        except (BlockingIOError, TimeoutError, OSError):
            data = b""

        if data:
            if len(data) == DONE.size and data[0] == MsgType.DONE:
                running = False
                continue
            if len(data) == DONE.size and data[0] == MsgType.PING:
                sock.sendto(pack_pong(), addr)
                continue

            if require_auth and auth and len(data) >= AUTH_CHALLENGE.size:
                if data[0] == MsgType.AUTH_REQUEST:
                    user_id = unpack_auth_request(data)
                    _, nonce = auth.issue_challenge(user_id)
                    pending_auth[addr] = user_id
                    sock.sendto(pack_auth_challenge(nonce), addr)
                    if verbose:
                        print(f"[hub:{port}] AUTH_CHALLENGE → {addr}", flush=True)
                    continue
                if data[0] == MsgType.AUTH_RESPONSE and addr in pending_auth:
                    user_id = pending_auth.pop(addr)
                    digest = unpack_auth_response(data)
                    if auth.verify_response(user_id, digest):
                        authenticated.add(addr)
                        sock.sendto(pack_auth_ok(), addr)
                        if verbose:
                            print(f"[hub:{port}] AUTH_OK {user_id} @ {addr}", flush=True)
                    else:
                        sock.sendto(pack_auth_fail(), addr)
                    continue

            if data[0] == MsgType.WORK_BATCH:
                if require_auth and addr not in authenticated:
                    continue
                for pkt, reply_addr in hub_process_work(
                    data, addr, cache, cache_hits=cache_hits
                ):
                    sock.sendto(pkt, reply_addr)
                    processed += 1
                continue

            if len(data) >= HEADER.size and data[0] == MsgType.WORK:
                if require_auth and addr not in authenticated:
                    continue
                seq, payload = unpack_work(data)
                while not ring.push((seq, payload, addr)):
                    pass

        for pkt, reply_addr in hub_drain_ring_burst(
            ring, cache, cache_hits, max_burst=32, min_batch=4
        ):
            sock.sendto(pkt, reply_addr)
            processed += 1

    sock.close()
    if verbose:
        print(
            f"[hub:{port}] done processed={processed} hits={cache_hits[0]}",
            flush=True,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="AntiPC UDP hub")
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--master-host", default="127.0.0.1")
    parser.add_argument("--master-port", type=int, default=19700)
    parser.add_argument("--master-key", default="", help="hex clave HMAC compartida")
    parser.add_argument("--require-auth", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    auth = None
    if args.require_auth and args.master_key:
        gate = NetworkAuthGate()
        gate.master_key = bytes.fromhex(args.master_key)
        auth = gate

    run_hub(
        args.port,
        (args.master_host, args.master_port),
        auth=auth,
        require_auth=args.require_auth and auth is not None,
        verbose=args.verbose,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())