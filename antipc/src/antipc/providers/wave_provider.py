"""WaveProvider v0.1 — inferencia desde latencia red/WiFi (R003)."""

from __future__ import annotations

import socket
import time
from dataclasses import dataclass

from hash_engine import get_backend
from runtime.plugin import Capability, Cost

from .compute_provider import ComputeProvider, ProviderResult


@dataclass
class WaveSample:
    host: str
    latency_us: int | None
    jitter_us: int
    inference: float
    backend: str


class WaveProvider(ComputeProvider):
    """
    No sustituye la ALU entera: mide el medio físico (RTT, jitter)
    y devuelve un escalar para decisiones del runtime.
    """

    def __init__(self, host: str = "8.8.8.8", port: int = 53, timeout_ms: int = 120) -> None:
        self.host = host
        self.port = port
        self.timeout_ms = timeout_ms
        self._last_us: int | None = None

    def measure_rtt_us(self) -> int | None:
        try:
            t0 = time.perf_counter_ns()
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(self.timeout_ms / 1000.0)
            sock.connect((self.host, self.port))
            sock.close()
            return (time.perf_counter_ns() - t0) // 1000
        except OSError:
            return None

    def sample(self) -> WaveSample:
        rtt = self.measure_rtt_us()
        jitter = 0
        if rtt is not None and self._last_us is not None:
            jitter = abs(rtt - self._last_us)
        self._last_us = rtt

        if rtt is None:
            inference = 0.0
        else:
            # Normalizado 0..1 — onda lenta = valor alto (más "inercia de red")
            inference = min(1.0, rtt / 200_000.0)

        return WaveSample(
            host=self.host,
            latency_us=rtt,
            jitter_us=jitter,
            inference=inference,
            backend=get_backend(),
        )

    def as_entropy_byte(self) -> int:
        s = self.sample()
        raw = (s.latency_us or time.perf_counter_ns()) ^ s.jitter_us
        return int(raw) & 0xFF

    def capability(self) -> Capability:
        return Capability(
            deterministic=False,
            distributable=False,
            cacheable=False,
            requires_network=True,
            requires_permission=True,
        )

    def estimate(self) -> Cost:
        return Cost(alu_units=0, memory_bytes=64)

    def execute(self) -> ProviderResult:
        s = self.sample()
        return ProviderResult(
            payload=f"{s.inference:.6f}".encode("ascii"),
            metadata={
                "host": s.host,
                "latency_us": s.latency_us or -1,
                "jitter_us": s.jitter_us,
                "inference": s.inference,
                "backend": s.backend,
            },
        )