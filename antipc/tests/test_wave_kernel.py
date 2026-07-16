"""Tests WaveMode + WaveProvider en FlowKernel."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "antipc"))

from plugins.k3_plugin import K3HashPlugin
from providers.wave_provider import WaveProvider, WaveSample
from runtime.kernel import FlowKernel
from runtime.modes import WaveMode


class _FakeWave(WaveProvider):
    def __init__(self, inference: float) -> None:
        super().__init__(host="127.0.0.1")
        self._inference = inference

    def sample(self) -> WaveSample:
        return WaveSample(
            host=self.host,
            latency_us=100_000,
            jitter_us=0,
            inference=self._inference,
            backend="test",
        )


def test_wave_mode_reuse_hit() -> None:
    k = FlowKernel(mode=WaveMode(), wave_provider=_FakeWave(0.8))
    k.network_permission = True
    k.register_plugin(K3HashPlugin())
    ref = k.ingest_raw(b"wave-test", label="t")
    k.submit_operation("K3_HASH", [ref])
    k.run_until_idle()
    k.submit_operation("K3_HASH", [ref])
    k.run_until_idle()
    st = k.stats()
    assert st["mode"] == "WAVE"
    assert st["knowledge_hits"] == 1
    assert st["wave_inference"] == 0.8


def test_wave_mode_defers_non_cacheable_alu() -> None:
    mode = WaveMode(defer_threshold=0.5, alu_defer_min=4)
    assert mode.should_defer(0.9, 16, cacheable=False) is True
    assert mode.should_defer(0.9, 16, cacheable=True) is False
    assert mode.should_defer(0.2, 16, cacheable=False) is False