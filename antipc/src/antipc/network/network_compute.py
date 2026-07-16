"""
Cómputo en red AntiPC — pipelining UDP, lotes WORK/RESULT, balanceo por hub.
"""

from __future__ import annotations

import socket
import time
from dataclasses import dataclass, field

from udp_protocol import (
    MsgType,
    pack_result,
    pack_result_batch,
    pack_work,
    pack_work_batch,
    unpack_result,
    unpack_result_batch,
    unpack_work,
    unpack_work_batch,
    HEADER,
    RESULT,
)


def tune_socket(sock: socket.socket, *, rcvbuf: int = 1 << 20, sndbuf: int = 1 << 20) -> None:
    try:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, rcvbuf)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, sndbuf)
    except OSError:
        pass
    try:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    except OSError:
        pass


def heavy_digest(payload: bytes) -> int:
    try:
        from native_engine import k3_heavy_hash

        return int(k3_heavy_hash(payload))
    except Exception:
        from hash_engine import k3_hash_buffer

        digest = k3_hash_buffer(payload)
        for _ in range(4):
            digest = k3_hash_buffer(digest.to_bytes(4, "little") + payload)
        return digest


@dataclass
class DispatchMetrics:
    items: int = 0
    received: int = 0
    batches_sent: int = 0
    batches_recv: int = 0
    singles_sent: int = 0
    fallback_local: int = 0
    elapsed_s: float = 0.0
    hub_rounds: int = 0

    @property
    def throughput(self) -> float:
        return self.received / self.elapsed_s if self.elapsed_s else 0.0


@dataclass
class PipelinedDispatcher:
    """
    Despacho WORK con ventana en vuelo + lotes UDP.
    Maximiza utilización de hubs (E(N) en P_util).
    """

    window: int = 128
    batch_size: int = 24
    use_batch: bool = True
    timeout_s: float = 12.0
    metrics: DispatchMetrics = field(default_factory=DispatchMetrics)

    def dispatch(
        self,
        sock: socket.socket,
        items: list[tuple[int, bytes]],
        hub_addrs: list[tuple[str, int]],
        *,
        route_fn: callable | None = None,
    ) -> dict[int, int]:
        if not items or not hub_addrs:
            return {}

        t0 = time.perf_counter()
        n = len(items)
        received: dict[int, int] = {}
        pending = n
        next_idx = 0
        inflight = 0
        self.metrics = DispatchMetrics(items=n)

        route = route_fn or (lambda seq, _payload: seq % len(hub_addrs))

        while pending > 0 or next_idx < n:
            while next_idx < n and inflight < self.window:
                remaining = n - next_idx
                if self.use_batch and remaining >= 4:
                    chunk_end = min(next_idx + self.batch_size, n)
                    chunk = items[next_idx:chunk_end]
                    hub_i = route(chunk[0][0], chunk[0][1])
                    addr = hub_addrs[hub_i % len(hub_addrs)]
                    sock.sendto(pack_work_batch(chunk), addr)
                    self.metrics.batches_sent += 1
                    inflight += len(chunk)
                    next_idx = chunk_end
                else:
                    seq, payload = items[next_idx]
                    hub_i = route(seq, payload)
                    addr = hub_addrs[hub_i % len(hub_addrs)]
                    sock.sendto(pack_work(seq, payload), addr)
                    self.metrics.singles_sent += 1
                    inflight += 1
                    next_idx += 1

            if pending == 0 and next_idx >= n:
                break

            try:
                data, _ = sock.recvfrom(65535)
            except (TimeoutError, BlockingIOError, OSError):
                if time.perf_counter() - t0 > self.timeout_s:
                    break
                continue

            if not data:
                continue

            if data[0] == MsgType.RESULT_BATCH:
                for seq, digest in unpack_result_batch(data):
                    if seq not in received:
                        received[seq] = digest
                        pending -= 1
                        inflight = max(0, inflight - 1)
                self.metrics.batches_recv += 1
                self.metrics.received = len(received)
                continue

            if len(data) >= RESULT.size and data[0] == MsgType.RESULT:
                seq, digest = unpack_result(data)
                if seq not in received:
                    received[seq] = digest
                    pending -= 1
                    inflight = max(0, inflight - 1)
                    self.metrics.received = len(received)

            if time.perf_counter() - t0 > self.timeout_s:
                break

        self.metrics.elapsed_s = time.perf_counter() - t0
        self.metrics.hub_rounds = self.metrics.batches_sent + self.metrics.singles_sent
        return received


def hub_process_work(
    data: bytes,
    addr: tuple[str, int],
    cache: dict[int, int],
    *,
    cache_hits: list[int],
) -> list[tuple[bytes, tuple[str, int]]]:
    """Procesa WORK o WORK_BATCH; retorna paquetes de respuesta a enviar."""
    out: list[tuple[bytes, tuple[str, int]]] = []

    if data[0] == MsgType.WORK_BATCH:
        batch = unpack_work_batch(data)
        results: list[tuple[int, int]] = []
        for seq, payload in batch:
            from hash_engine import k3_hash_buffer

            key = k3_hash_buffer(payload)
            if key in cache:
                digest = cache[key]
                cache_hits[0] += 1
            else:
                digest = heavy_digest(payload)
                cache[key] = digest
            results.append((seq, digest))
        if results:
            out.append((pack_result_batch(results), addr))
        return out

    if len(data) >= HEADER.size and data[0] == MsgType.WORK:
        seq, payload = unpack_work(data)
        from hash_engine import k3_hash_buffer

        key = k3_hash_buffer(payload)
        if key in cache:
            digest = cache[key]
            cache_hits[0] += 1
        else:
            digest = heavy_digest(payload)
            cache[key] = digest
        out.append((pack_result(seq, digest), addr))
    return out


def hub_drain_ring_burst(
    ring,
    cache: dict[int, int],
    cache_hits: list[int],
    *,
    max_burst: int = 32,
    min_batch: int = 4,
) -> list[tuple[bytes, tuple[str, int]]]:
    """Vacía el ring del hub agrupando RESULT_BATCH cuando hay varios ítems."""
    pending: list[tuple[int, bytes, tuple[str, int]]] = []
    for _ in range(max_burst):
        item = ring.pop()
        if item is None:
            break
        pending.append(item)

    if not pending:
        return []

    by_addr: dict[tuple[str, int], list[tuple[int, int]]] = {}
    from hash_engine import k3_hash_buffer

    for seq, payload, reply_addr in pending:
        key = k3_hash_buffer(payload)
        if key in cache:
            digest = cache[key]
            cache_hits[0] += 1
        else:
            digest = heavy_digest(payload)
            cache[key] = digest
        by_addr.setdefault(reply_addr, []).append((seq, digest))

    replies: list[tuple[bytes, tuple[str, int]]] = []
    for addr, results in by_addr.items():
        if len(results) >= min_batch:
            replies.append((pack_result_batch(results), addr))
        else:
            for seq, digest in results:
                replies.append((pack_result(seq, digest), addr))
    return replies


def classic_dispatch(
    sock: socket.socket,
    items: list[tuple[int, bytes]],
    hub_addrs: list[tuple[str, int]],
    *,
    timeout_s: float = 12.0,
    burst_send: int = 32,
) -> tuple[dict[int, int], DispatchMetrics]:
    """Despacho secuencial legacy (pausa cada burst_send WORK)."""
    t0 = time.perf_counter()
    n = len(items)
    received: dict[int, int] = {}
    pending = n
    next_send = 0
    metrics = DispatchMetrics(items=n)

    while pending > 0 or next_send < n:
        sent_burst = 0
        while next_send < n and sent_burst < burst_send:
            seq, payload = items[next_send]
            addr = hub_addrs[next_send % len(hub_addrs)]
            sock.sendto(pack_work(seq, payload), addr)
            metrics.singles_sent += 1
            next_send += 1
            sent_burst += 1

        try:
            data, _ = sock.recvfrom(65535)
        except (TimeoutError, BlockingIOError, OSError):
            if time.perf_counter() - t0 > timeout_s:
                break
            continue

        if len(data) >= RESULT.size and data[0] == MsgType.RESULT:
            seq, digest = unpack_result(data)
            if seq not in received:
                received[seq] = digest
                pending -= 1
                metrics.received = len(received)

        if time.perf_counter() - t0 > timeout_s:
            break

    metrics.elapsed_s = time.perf_counter() - t0
    return received, metrics


@dataclass
class NetworkBenchResult:
    items: int
    hubs: int
    classic: DispatchMetrics
    pipelined: DispatchMetrics
    parity_ok: bool

    @property
    def speedup(self) -> float:
        c = self.classic.throughput
        p = self.pipelined.throughput
        return p / c if c else 0.0


def run_network_bench(
    fabric,
    *,
    count: int = 512,
    hubs: int = 4,
    window: int = 128,
    batch_size: int = 24,
    payload_size: int = 64,
    timeout_s: float = 15.0,
) -> NetworkBenchResult:
    """Compara despacho clásico vs pipelined sobre hubs locales autenticados."""
    if not any(h.authenticated for h in fabric.hubs):
        fabric.start_hubs(hubs)
        fabric.authenticate_all()

    hub_addrs = [
        fabric.hub_address(h) for h in fabric.hubs if h.authenticated
    ]
    if not hub_addrs:
        raise RuntimeError("ningún hub autenticado para bench")

    sock = fabric._open_socket()
    tune_socket(sock)
    sock.settimeout(0.002)

    payloads = [
        f"antipc-bench-{i}".encode() + bytes([i % 256] * payload_size)
        for i in range(count)
    ]
    items = list(enumerate(payloads))

    classic_results, classic_metrics = classic_dispatch(
        sock, items, hub_addrs, timeout_s=timeout_s
    )

    dispatcher = PipelinedDispatcher(
        window=window,
        batch_size=batch_size,
        use_batch=True,
        timeout_s=timeout_s,
    )
    pipelined_results = dispatcher.dispatch(sock, items, hub_addrs)
    pipelined_metrics = dispatcher.metrics

    parity = classic_results == pipelined_results

    return NetworkBenchResult(
        items=count,
        hubs=len(hub_addrs),
        classic=classic_metrics,
        pipelined=pipelined_metrics,
        parity_ok=parity,
    )