"""L3 — UDP Fabric enlazado a AntiPCStack + HMAC."""

from __future__ import annotations

import os
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, field

from architecture.network_auth import NetworkAuthGate
from udp_protocol import (
    MsgType,
    pack_auth_challenge,
    pack_auth_ok,
    pack_auth_request,
    pack_auth_response,
    pack_done,
    pack_ping,
    pack_result,
    pack_work,
    unpack_auth_challenge,
    unpack_result,
    AUTH_CHALLENGE,
    RESULT,
)

BASE_HUB_PORT = 19701
MASTER_PORT = 19700
AUTH_TIMEOUT = 2.0


@dataclass
class HubProcess:
    port: int
    process: subprocess.Popen[bytes]
    authenticated: bool = False


@dataclass
class UdpFabric:
    """Maestro UDP: spawn hubs, HMAC, despacho WORK/RESULT."""

    auth: NetworkAuthGate
    user_id: str
    master_host: str = "127.0.0.1"
    master_port: int = MASTER_PORT
    base_hub_port: int = BASE_HUB_PORT
    hubs: list[HubProcess] = field(default_factory=list)
    _sock: socket.socket | None = None
    _script_dir: str = field(default_factory=lambda: os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    def _master_key_hex(self) -> str:
        return self.auth.master_key.hex()

    def _open_socket(self) -> socket.socket:
        if self._sock is not None:
            return self._sock
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        sock.bind((self.master_host, self.master_port))
        sock.settimeout(0.01)
        self._sock = sock
        return sock

    def start_hubs(self, count: int) -> None:
        hub_script = os.path.join(self._script_dir, "hub_node.py")
        key_hex = self._master_key_hex()
        for i in range(count):
            port = self.base_hub_port + i
            proc = subprocess.Popen(
                [
                    sys.executable, hub_script,
                    "--port", str(port),
                    "--master-host", self.master_host,
                    "--master-port", str(self.master_port),
                    "--master-key", key_hex,
                    "--require-auth",
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            self.hubs.append(HubProcess(port, proc))
        time.sleep(0.4)
        self._open_socket()

    def authenticate_hub(self, hub_port: int) -> bool:
        """Handshake HMAC UDP con un hub (cliente → servidor hub)."""
        sock = self._open_socket()
        hub_addr = (self.master_host, hub_port)
        deadline = time.perf_counter() + AUTH_TIMEOUT

        sock.sendto(pack_auth_request(self.user_id), hub_addr)

        nonce: bytes | None = None
        while time.perf_counter() < deadline:
            try:
                data, addr = sock.recvfrom(1024)
            except (TimeoutError, BlockingIOError, OSError):
                continue
            if addr[1] != hub_port or len(data) < AUTH_CHALLENGE.size:
                continue
            if data[0] == MsgType.AUTH_CHALLENGE:
                nonce = unpack_auth_challenge(data)
                break

        if nonce is None:
            return False

        response = self.auth.client_response(self.user_id, nonce)
        sock.sendto(pack_auth_response(response), hub_addr)

        while time.perf_counter() < deadline:
            try:
                data, addr = sock.recvfrom(1024)
            except (TimeoutError, BlockingIOError, OSError):
                continue
            if addr[1] == hub_port and data[0] == MsgType.AUTH_OK:
                return True
            if addr[1] == hub_port and data[0] == MsgType.AUTH_FAIL:
                return False
        return False

    def authenticate_all(self) -> int:
        ok = 0
        for hub in self.hubs:
            hub.authenticated = self.authenticate_hub(hub.port)
            if hub.authenticated:
                ok += 1
        return ok

    def wait_pong(self, hub_port: int) -> bool:
        sock = self._open_socket()
        sock.sendto(pack_ping(), (self.master_host, hub_port))
        deadline = time.perf_counter() + 1.0
        while time.perf_counter() < deadline:
            try:
                data, _ = sock.recvfrom(256)
                if data and data[0] == MsgType.PONG:
                    return True
            except (TimeoutError, BlockingIOError, OSError):
                pass
        return False

    def dispatch_hashes(self, payloads: list[bytes]) -> dict[int, int]:
        """Envía WORK a hubs autenticados; devuelve {seq: digest}."""
        if not any(h.authenticated for h in self.hubs):
            raise RuntimeError("no authenticated hubs")

        sock = self._open_socket()
        hub_addrs = [
            (self.master_host, h.port) for h in self.hubs if h.authenticated
        ]
        received: dict[int, int] = {}
        pending = len(payloads)
        next_send = 0

        while pending > 0 or next_send < len(payloads):
            while next_send < len(payloads):
                addr = hub_addrs[next_send % len(hub_addrs)]
                sock.sendto(pack_work(next_send, payloads[next_send]), addr)
                next_send += 1
                if next_send % 32 == 0:
                    break

            try:
                data, _ = sock.recvfrom(4096)
            except (TimeoutError, BlockingIOError, OSError):
                continue

            if len(data) >= RESULT.size and data[0] == MsgType.RESULT:
                seq, digest = unpack_result(data)
                if seq not in received:
                    received[seq] = digest
                    pending -= 1

        return received

    def shutdown(self) -> None:
        if self._sock:
            for hub in self.hubs:
                try:
                    self._sock.sendto(pack_done(), (self.master_host, hub.port))
                except OSError:
                    pass
            time.sleep(0.2)
            self._sock.close()
            self._sock = None
        for hub in self.hubs:
            hub.process.terminate()
            try:
                hub.process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                hub.process.kill()
        self.hubs.clear()