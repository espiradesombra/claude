"""
Cluster MMO L3 — validación heavy K3 en hubs WORK/RESULT por shard.

Maestro: ingest UDP (B) + existencia Grafcet local (D parcial).
Hubs:   heavy hash remoto (misma ALU que hub_node.py / AntiPCStack L3).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from architecture.antipc_stack import AntiPCStack
from architecture.udp_fabric import UdpFabric
from network.network_compute import PipelinedDispatcher, tune_socket


def _shard_for_payload(payload: bytes, num_shards: int, fallback: int) -> int:
    from .game_protocol import unpack_state

    pkt = unpack_state(payload)
    if pkt is None:
        return fallback % max(1, num_shards)
    return pkt.entity_id % max(1, num_shards)


@dataclass
class GameClusterFabric:
    """Arranca hubs L3, autentica HMAC y despacha WORK por shard de entidad."""

    num_shards: int = 4
    hubs: int = 4
    hub_hosts: list[str] | None = None
    hub_base_port: int = 19701
    stack: AntiPCStack = field(default_factory=AntiPCStack)
    fabric: UdpFabric | None = None
    authenticated: int = 0
    validations_sent: int = 0
    validations_ok: int = 0
    validations_fallback: int = 0
    _timeout_s: float = 8.0

    def start(self) -> dict:
        if not self.stack.request_network_permission():
            raise PermissionError("HMAC local falló — red no autorizada")
        self.fabric = self.stack.fabric
        if self.fabric is None:
            raise RuntimeError("UdpFabric no inicializado")
        if self.hub_hosts:
            n = len(self.hub_hosts)
            self.fabric.register_remote_hubs(
                self.hub_hosts, base_port=self.hub_base_port
            )
            hub_count = n
        else:
            hub_count = max(self.hubs, self.num_shards)
            self.fabric.start_hubs(hub_count)
        self.authenticated = self.fabric.authenticate_all()
        if self.authenticated == 0:
            raise RuntimeError("ningún hub autenticado")
        return {
            "hubs": hub_count,
            "authenticated": self.authenticated,
            "shards": self.num_shards,
            "remote": bool(self.hub_hosts),
            "endpoints": [
                {"host": h.host, "port": h.port, "remote": h.remote}
                for h in (self.fabric.hubs if self.fabric else [])
            ],
        }

    def shutdown(self) -> None:
        if self.stack:
            self.stack.stop_network()
        self.fabric = None

    def _hub_addrs(self) -> list[tuple[str, int]]:
        if self.fabric is None:
            return []
        return [
            self.fabric.hub_address(h)
            for h in self.fabric.hubs
            if h.authenticated
        ]

    def validate_heavy_batch(
        self,
        payloads: list[bytes],
        *,
        num_shards: int | None = None,
    ) -> list[int]:
        """
        Envía cada payload al hub del shard (entity_id % shards).
        Retorna lista de digest heavy en el mismo orden.
        """
        if not payloads:
            return []

        hub_addrs = self._hub_addrs()
        if not hub_addrs:
            return self._fallback_heavy(payloads)

        shards = max(1, num_shards or self.num_shards)
        fabric = self.fabric
        assert fabric is not None
        sock = fabric._open_socket()
        tune_socket(sock)
        n = len(payloads)
        items = list(enumerate(payloads))

        def route_fn(seq: int, payload: bytes) -> int:
            return _shard_for_payload(payload, shards, seq)

        dispatcher = PipelinedDispatcher(
            window=fabric.dispatch_window,
            batch_size=fabric.dispatch_batch,
            use_batch=True,
            timeout_s=self._timeout_s,
        )
        received = dispatcher.dispatch(sock, items, hub_addrs, route_fn=route_fn)
        self.validations_sent += n

        digests: list[int] = []
        for seq in range(n):
            if seq in received:
                digests.append(received[seq])
                self.validations_ok += 1
            else:
                digests.append(self._heavy_local(payloads[seq]))
                self.validations_fallback += 1
        return digests

    @staticmethod
    def _heavy_local(payload: bytes) -> int:
        try:
            from native_engine import k3_heavy_hash

            return int(k3_heavy_hash(payload))
        except Exception:
            from grafcet import _heavy_hash

            return _heavy_hash(payload)

    def _fallback_heavy(self, payloads: list[bytes]) -> list[int]:
        self.validations_fallback += len(payloads)
        return [self._heavy_local(p) for p in payloads]