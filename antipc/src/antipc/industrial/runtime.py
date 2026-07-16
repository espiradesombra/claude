"""IndustrialRuntime — perfil + resolver + workers + telemetría."""

from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from time import perf_counter

from runtime.event_bus import Event, EventType
from runtime.kernel import FlowKernel
from runtime.modes import BenchmarkMode, LocalMode
from runtime.plugin import PluginContext
from runtime.plugin_manager import PluginManager
from runtime.profile import ExecutionProfile
from runtime.reference import Reference, ReferenceRecord, ReferenceState
from runtime.resolver import KnowledgeResolver, ResolutionStatus
from runtime.scheduler import PendingOperation
from runtime.telemetry import IndustrialTelemetry


@dataclass
class IndustrialRuntime:
    profile: ExecutionProfile = field(default_factory=ExecutionProfile.server)
    kernel: FlowKernel = field(default_factory=FlowKernel)
    manager: PluginManager = field(default_factory=PluginManager)
    telemetry: IndustrialTelemetry = field(default_factory=IndustrialTelemetry)
    resolver: KnowledgeResolver | None = None
    _lock: threading.Lock = field(default_factory=threading.Lock)

    def __post_init__(self) -> None:
        if self.profile.benchmark_enabled:
            self.kernel.mode = BenchmarkMode()
        else:
            self.kernel.mode = LocalMode()
        self.kernel.network_permission = self.profile.network_permission
        self.resolver = KnowledgeResolver(
            self.kernel.knowledge,
            self.profile,
        )

    def bootstrap(self, *plugins) -> None:
        for p in plugins:
            self.manager.register(p)
            self.kernel.register_plugin(p)

    def fingerprint_file(self, path: str) -> Reference | None:
        """Pipeline industrial: ruta → K3_FILE con resolver."""
        path_ref = self.kernel.ingest_raw(path.encode("utf-8"), label="PATH")
        return self._run_plugin("K3_FILE", [path_ref], path.encode("utf-8"))

    def _run_plugin(self, name: str, input_refs: list[Reference],
                    payload: bytes) -> Reference | None:
        plugin = self.manager.get(name)
        ctx = PluginContext(input_refs, [payload] if payload else
                            [r.payload for r in input_refs if r.payload])
        if not plugin.validate(ctx):
            self.telemetry.log("validate_fail", name, "", "FAILED")
            return None

        sig = plugin.signature(ctx)
        cost = plugin.estimate(ctx).alu_units
        resolution = self.resolver.resolve(sig, cost)

        if resolution.status in (ResolutionStatus.FOUND_LOCAL, ResolutionStatus.FOUND_SHARED):
            self.manager.record_reuse(name)
            self.kernel.scheduler.cancelled += 1
            self.telemetry.log("reuse", name, sig, resolution.status.name, alu_saved=True)
            if resolution.reference:
                return resolution.reference
            known = self.kernel.knowledge.lookup(sig)
            return known

        if resolution.status == ResolutionStatus.WAIT and not self.profile.distributed_execution:
            resolution = resolution  # fall through to execute locally

        t0 = perf_counter()
        ref, _ = self.manager.execute_timed(name, ctx)
        elapsed = (perf_counter() - t0) * 1000
        self.kernel.references[ref.ref_id] = ReferenceRecord(ref, ReferenceState.PUBLISHED)
        self.kernel.knowledge.publish(sig, ref)
        self.kernel.scheduler.executed += 1
        self.kernel.bus.emit(Event(EventType.REFERENCE_PUBLISHED, ref.ref_id, sig))
        self.telemetry.log("execute", name, sig, "EXECUTE", elapsed_ms=elapsed)
        return ref

    def scan_files(self, paths: list[str], workers: int = 4) -> list[Reference]:
        """Escaneo distribuido — simula hubs industriales con thread pool."""
        results: list[Reference | None] = [None] * len(paths)

        if not self.profile.distributed_execution or workers <= 1:
            out = []
            for p in paths:
                r = self.fingerprint_file(p)
                if r:
                    out.append(r)
            return out

        def task(idx: int, path: str) -> tuple[int, Reference | None]:
            return idx, self.fingerprint_file(path)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(task, i, p) for i, p in enumerate(paths)]
            for fut in as_completed(futures):
                idx, ref = fut.result()
                results[idx] = ref
        return [r for r in results if r is not None]

    def run_dedup_report(self, file_refs: list[Reference]) -> Reference | None:
        if not file_refs:
            return None
        ctx = PluginContext(file_refs, [])
        plugin = self.manager.get("K3_DEDUP")
        if not plugin.validate(ctx):
            return None
        ref, _ = self.manager.execute_timed("K3_DEDUP", ctx)
        self.telemetry.log("dedup", "K3_DEDUP", ref.signature, "REPORT")
        return ref

    def report(self) -> dict:
        k = self.kernel.stats()
        k["profile"] = {
            "distributed": self.profile.distributed_execution,
            "network": self.profile.network_enabled,
            "research": self.profile.research_mode,
            "battery": self.profile.battery_optimizer,
        }
        k["plugin_stats"] = {
            name: {
                "executions": s.executions,
                "reuse_hits": s.reuse_hits,
                "avg_ms": round(s.average_time_ms, 2),
                "failures": s.failures,
            }
            for name, s in self.manager.stats.items()
        }
        return k