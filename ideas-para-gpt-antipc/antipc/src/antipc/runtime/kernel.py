"""Flow Kernel — event-driven runtime v0.1.0."""

from __future__ import annotations

from dataclasses import dataclass, field

from .event_bus import Event, EventBus, EventType
from .knowledge import KnowledgeBuffer
from .modes import RuntimeMode, LocalMode, RuntimeModeName, WaveMode
from .plugin import IPlugin, PluginContext
from .reference import Reference, ReferenceRecord, ReferenceState
from .scheduler import PendingOperation, Scheduler


@dataclass
class FlowKernel:
    """
    AntiPC Flow Kernel.

    Moves knowledge; never knows what HASH/OCR/AI mean — only plugins.
    """

    mode: RuntimeMode = field(default_factory=LocalMode)
    knowledge: KnowledgeBuffer = field(default_factory=KnowledgeBuffer)
    bus: EventBus = field(default_factory=EventBus)
    scheduler: Scheduler = field(default_factory=Scheduler)
    plugins: dict[str, IPlugin] = field(default_factory=dict)
    references: dict[str, ReferenceRecord] = field(default_factory=dict)
    network_permission: bool = False
    wave_provider: object | None = None
    _last_wave_inference: float = field(default=0.0, repr=False)
    _wave_deferred: int = field(default=0, repr=False)
    _wave_defer_streak: int = field(default=0, repr=False)
    _wave_defer_streak_max: int = 3

    def register_plugin(self, plugin: IPlugin) -> None:
        info = plugin.info()
        self.plugins[info.name] = plugin

    def _refresh_wave_inference(self) -> float:
        if self.mode.name != RuntimeModeName.WAVE or self.wave_provider is None:
            return 0.0
        sample = self.wave_provider.sample()
        self._last_wave_inference = sample.inference
        return self._last_wave_inference

    def ingest_raw(self, payload: bytes, label: str = "RAW") -> Reference:
        """Acquire step: wrap bytes as initial reference."""
        from .signature import make_signature

        content_sig = make_signature("ACQUIRE", [payload.hex()])
        ref = Reference.create(
            signature=content_sig,
            payload=payload,
            metadata={"label": label},
        )
        self.references[ref.ref_id] = ReferenceRecord(ref, ReferenceState.READY)
        self.bus.emit(Event(EventType.REFERENCE_CREATED, ref.ref_id, ref.signature))
        return ref

    def submit_operation(self, plugin_name: str, input_refs: list[Reference]) -> str:
        plugin = self.plugins[plugin_name]
        payloads = [r.payload for r in input_refs if r.payload is not None]
        ctx = PluginContext(input_refs, payloads)
        if not plugin.validate(ctx):
            raise ValueError(f"invalid input for plugin {plugin_name}")

        sig = plugin.signature(ctx)
        cache = self.mode.cache_policy()

        if cache.reuse_knowledge:
            known = self.knowledge.lookup(sig)
            if known is not None:
                self.scheduler.cancel_if_reusable(sig)
                self.references[known.ref_id] = ReferenceRecord(
                    known, ReferenceState.PUBLISHED
                )
                self.bus.emit(
                    Event(EventType.REFERENCE_PUBLISHED, known.ref_id, sig),
                    merge_equal=True,
                )
                return known.ref_id

        if self.mode.name == RuntimeModeName.WAVE:
            self._refresh_wave_inference()

        est = plugin.estimate(ctx)
        priority = est.alu_units
        if hasattr(self.mode, "priority_boost"):
            priority += self.mode.priority_boost(
                plugin.capability().cacheable, False
            )
            if (
                isinstance(self.mode, WaveMode)
                and plugin.capability().cacheable
                and self._last_wave_inference > 0.0
            ):
                priority += self._last_wave_inference * 2.0

        op = PendingOperation(
            plugin=plugin,
            ctx=ctx,
            signature=sig,
            required=len(input_refs),
            available=len(input_refs),
            priority=priority,
            alu_units=est.alu_units,
        )
        self.scheduler.submit(op)
        self.bus.emit(Event(EventType.OPERATION_READY, sig, sig))
        return sig

    def run_once(self) -> Reference | None:
        """Process one READY operation (worker pool unit)."""
        net = self.mode.network_policy()
        if net.enabled and net.require_permission and not self.network_permission:
            return None

        op = self.scheduler.peek_ready()
        if op is None:
            return None

        if isinstance(self.mode, WaveMode) and self.wave_provider is not None:
            inference = self._refresh_wave_inference()
            if (
                self._wave_defer_streak < self._wave_defer_streak_max
                and self.mode.should_defer(
                    inference,
                    op.alu_units,
                    op.plugin.capability().cacheable,
                )
            ):
                self._wave_deferred += 1
                self._wave_defer_streak += 1
                return None

        op = self.scheduler.pop_ready()
        if op is None:
            return None

        self._wave_defer_streak = 0
        ref = op.plugin.execute(op.ctx)
        self.references[ref.ref_id] = ReferenceRecord(ref, ReferenceState.PUBLISHED)
        self.knowledge.publish(ref.signature, ref)
        self.scheduler.executed += 1
        self.bus.emit(Event(EventType.REFERENCE_PUBLISHED, ref.ref_id, ref.signature))
        return ref

    def run_until_idle(self, max_steps: int = 10_000) -> int:
        steps = 0
        while steps < max_steps:
            result = self.run_once()
            if result is None and not self.scheduler.ready_ops():
                break
            if result is not None:
                steps += 1
        return steps

    def stats(self) -> dict[str, int | float | None]:
        out: dict[str, int | float | None] = {
            "mode": self.mode.name.name,
            "plugins": len(self.plugins),
            "references": len(self.references),
            "executed": self.scheduler.executed,
            "cancelled": self.scheduler.cancelled,
            "knowledge_queries": self.knowledge.queries,
            "knowledge_hits": self.knowledge.hits,
            "events_coalesced": self.bus.coalesced,
            "wave_inference": None,
            "wave_deferred": 0,
        }
        if self.mode.name == RuntimeModeName.WAVE:
            out["wave_inference"] = self._last_wave_inference
            out["wave_deferred"] = self._wave_deferred
        return out