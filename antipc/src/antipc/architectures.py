"""Four AntiPC architectures for reproducible CPU comparison."""

from __future__ import annotations

import queue
import threading
import time
from dataclasses import dataclass, field
from typing import Callable

try:
    from .grafcet import GrafcetEngine, Packet
    from .hash_engine import k3_hash
    from .ring_buffer import SPSCRingBuffer
except ImportError:
    from grafcet import GrafcetEngine, Packet
    from hash_engine import k3_hash
    from ring_buffer import SPSCRingBuffer


@dataclass
class BenchmarkResult:
    name: str
    packets: int
    elapsed_s: float
    processed: int
    cache_hits: int = 0
    hubs: int = 1

    @property
    def throughput(self) -> float:
        return self.processed / self.elapsed_s if self.elapsed_s else 0.0

    @property
    def cpu_load_proxy(self) -> float:
        """Relative ALU pressure: lower is better (inverse throughput per copy)."""
        return (self.elapsed_s / max(self.processed, 1)) * 1000.0


def _make_packets(count: int, payload_size: int, repeat_ratio: float = 0.15) -> list[Packet]:
    """Synthetic UDP-like payloads with controlled duplication for cache tests."""
    import random

    rng = random.Random(42)
    base_pool = [rng.randbytes(payload_size) for _ in range(max(16, count // 10))]
    packets: list[Packet] = []
    for i in range(count):
        if rng.random() < repeat_ratio:
            payload = base_pool[rng.randrange(len(base_pool))]
        else:
            payload = rng.randbytes(payload_size)
        packets.append(Packet(seq=i, payload=payload))
    return packets


def _heavy_hash(payload: bytes) -> int:
    """Simulates ALU-bound HASH→INDEX→QUERY pipeline."""
    digest = k3_hash(payload)
    for _ in range(4):
        digest = k3_hash(digest.to_bytes(4, "little") + payload)
    return digest


# --- Architecture A: conventional (memcpy + locked queue) ---

def run_conventional(packets: list[Packet]) -> BenchmarkResult:
    q: queue.Queue[bytes] = queue.Queue()
    processed = 0

    def worker() -> None:
        nonlocal processed
        while True:
            try:
                data = q.get(timeout=0.05)
            except queue.Empty:
                break
            _heavy_hash(data)
            processed += 1
            q.task_done()

    t0 = time.perf_counter()
    thread = threading.Thread(target=worker, daemon=True)
    thread.start()
    for pkt in packets:
        copy = bytes(pkt.payload)  # explicit memcpy simulation
        q.put(copy)
    q.join()
    thread.join(timeout=1.0)
    elapsed = time.perf_counter() - t0
    return BenchmarkResult("A — Convencional (memcpy + cola)", len(packets), elapsed, processed)


# --- Architecture B: lock-free ring buffer ---

def run_lockfree(packets: list[Packet]) -> BenchmarkResult:
    ring: SPSCRingBuffer[bytes] = SPSCRingBuffer(4096)
    processed = 0
    done = threading.Event()

    def worker() -> None:
        nonlocal processed
        while not done.is_set() or len(ring) > 0:
            item = ring.pop()
            if item is None:
                time.sleep(0)
                continue
            _heavy_hash(item)
            processed += 1

    t0 = time.perf_counter()
    thread = threading.Thread(target=worker, daemon=True)
    thread.start()
    for pkt in packets:
        data = pkt.payload
        while not ring.push(data):
            time.sleep(0)
    done.set()
    thread.join(timeout=2.0)
    elapsed = time.perf_counter() - t0
    return BenchmarkResult("B — Lock-free ring buffer", len(packets), elapsed, processed)


# --- Architecture C: distributed hubs ---

def run_distributed(packets: list[Packet], hubs: int = 4) -> BenchmarkResult:
    """Partitions packets across hub workers; master only aggregates."""
    buckets: list[list[Packet]] = [[] for _ in range(hubs)]
    for i, pkt in enumerate(packets):
        buckets[i % hubs].append(pkt)

    processed = 0
    cache: dict[int, int] = {}
    lock = threading.Lock()

    def hub_worker(batch: list[Packet]) -> None:
        nonlocal processed
        local_cache: dict[int, int] = {}
        local_count = 0
        for pkt in batch:
            key = k3_hash(pkt.payload)
            if key in local_cache:
                digest = local_cache[key]
            else:
                digest = _heavy_hash(pkt.payload)
                local_cache[key] = digest
            local_count += 1
        with lock:
            processed += local_count
            cache.update(local_cache)

    t0 = time.perf_counter()
    threads = [
        threading.Thread(target=hub_worker, args=(bucket,), daemon=True)
        for bucket in buckets if bucket
    ]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    elapsed = time.perf_counter() - t0
    return BenchmarkResult(
        f"C — Distribuido ({hubs} hubs)",
        len(packets),
        elapsed,
        processed,
        cache_hits=len(packets) - len(cache),
        hubs=hubs,
    )


# --- Architecture D: Grafcet + existence matrix + redundancy ---

def run_grafcet(packets: list[Packet], batch_size: int = 32) -> BenchmarkResult:
    engine = GrafcetEngine(batch_size=batch_size)
    t0 = time.perf_counter()
    for pkt in packets:
        engine.ingest(pkt)
    engine.finalize()
    elapsed = time.perf_counter() - t0
    return BenchmarkResult(
        "D — Grafcet + existencia + redundancia",
        len(packets),
        elapsed,
        engine.processed,
        cache_hits=engine.cache_hits,
    )