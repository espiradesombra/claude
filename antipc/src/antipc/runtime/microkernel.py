"""AntiPC Microkernel v0.1.2-alpha — Boot → Registry → KOP → Ledger → Reuse."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .event_bus import EventBus
from .knowledge import KnowledgeBuffer
from .kop import KOPId
from .kop_executor import KOPExecutor
from .kops.boot_kop import BootKOP
from .kops.hash_kop import HashKOP
from .kops.identity_kop import IdentityKOP
from .ledger import KnowledgeLedger
from .metrics import MetricsCollector
from .policy import PolicyEngine
from .registry import Registry
from .runtime_context import KOPStatus, RuntimeContext
from .scheduler import Scheduler


VERSION = "0.1.2-alpha"


@dataclass
class AntiPCMicroKernel:
    registry: Registry = field(default_factory=Registry)
    events: EventBus = field(default_factory=EventBus)
    scheduler: Scheduler = field(default_factory=Scheduler)
    knowledge: KnowledgeBuffer = field(default_factory=KnowledgeBuffer)
    ledger: KnowledgeLedger = field(default_factory=KnowledgeLedger)
    policy: PolicyEngine = field(default_factory=PolicyEngine)
    metrics: MetricsCollector = field(default_factory=MetricsCollector)
    executor: KOPExecutor = field(default_factory=KOPExecutor)
    _booted: bool = False
    _ctx: RuntimeContext | None = None

    def _context(self) -> RuntimeContext:
        if self._ctx is None:
            self._ctx = RuntimeContext(
                registry=self.registry,
                scheduler=self.scheduler,
                events=self.events,
                knowledge=self.knowledge,
                ledger=self.ledger,
                policy=self.policy,
                metrics=self.metrics,
            )
        return self._ctx

    def boot(self, verbose: bool = True) -> None:
        if verbose:
            self._line("======================================")
            self._line(f"AntiPC Runtime {VERSION}")
            self._line("======================================")
            self._line("")
            self._line("Booting Kernel...")
        self.executor.register_handler(IdentityKOP())
        self.executor.register_handler(HashKOP())
        checks = [
            ("Registry", lambda: self.registry is not None),
            ("EventBus", lambda: self.events is not None),
            ("Scheduler", lambda: self.scheduler is not None),
            ("KOP Registry", lambda: len(self.executor.describe(KOPId.IDENTITY)) > 0),
            ("HashKOP", lambda: self.executor.describe(KOPId.HASH)["registered"]),
            ("Ledger", lambda: self.ledger is not None),
            ("Policy", lambda: self.policy is not None),
            ("Metrics", lambda: self.metrics is not None),
        ]
        for name, fn in checks:
            ok = fn()
            if verbose:
                self._line(f"{name:.<22}{'OK' if ok else 'FAIL'}")
            if not ok:
                raise RuntimeError(f"boot failed: {name}")
        self._booted = True
        if verbose:
            self._line("")
            self._line("Kernel Ready")
            self._line("")

    def run_boot_kop(self, verbose: bool = True) -> bool:
        if not self._booted:
            self.boot(verbose=False)
        if verbose:
            self._line("Executing BootKOP...")
        result = BootKOP().execute(self._context(), None)
        ok = result.status == KOPStatus.OK
        if verbose:
            self._line(f"BootKOP{'':.<16}{'OK' if ok else 'FAIL'}")
        return ok

    def create_knowledge(self, payload: bytes, verbose: bool = True) -> str | None:
        ctx = self._context()
        result = self.executor.execute(ctx, KOPId.IDENTITY, payload)
        if result.status not in (KOPStatus.OK, KOPStatus.SKIPPED):
            if verbose:
                self._line(f"IdentityKOP FAIL: {result.message}")
            return None
        if verbose:
            self._line(f"IdentityKOP{'':.<12}{'OK' if result.status == KOPStatus.OK else 'REUSE'}")
            self._line(f"Knowledge Objects.....{self.registry.count()}")
        if result.data is not None and hasattr(result.data, "id"):
            return result.data.id()
        return None

    def lookup(self, object_id: str) -> bool:
        obj = self.registry.find(object_id)
        self.metrics.record_query(obj is not None)
        return obj is not None

    def hash_payload(self, payload: bytes, verbose: bool = True) -> str | None:
        if not self._booted:
            self.boot(verbose=False)
        ctx = self._context()
        result = self.executor.execute(ctx, KOPId.HASH, payload)
        if result.status not in (KOPStatus.OK, KOPStatus.SKIPPED):
            if verbose:
                self._line(f"HashKOP FAIL: {result.message}")
            return None
        ref = result.data
        digest = ""
        if ref is not None and hasattr(ref, "metadata"):
            digest = str(ref.metadata.get("digest_hex", ""))
        if not digest and result.status == KOPStatus.SKIPPED:
            handler = self.executor.get_handler(KOPId.HASH)
            if handler is not None:
                sig = handler.knowledge_signature(payload)
                if sig is not None:
                    known = ctx.knowledge.lookup(sig)
                    if known is not None:
                        digest = str(known.metadata.get("digest_hex", ""))
        if verbose:
            tag = "REUSE" if result.status == KOPStatus.SKIPPED else "OK"
            self._line(f"HashKOP{'':.<16}{tag}  {digest}")
        return digest or None

    def export_metrics(self, path: Path | str) -> Path:
        return self.metrics.export(path, version=VERSION)

    def shutdown(self, verbose: bool = True, metrics_path: Path | str | None = None) -> None:
        if verbose:
            self._line("")
            self._line(f"Ledger entries........{self.ledger.count()}")
            self._line(f"P_util................{self.metrics.p_util:.2%}")
            self._line(f"ALU saved.............{self.metrics.alu_saved_pct:.1f}%")
            self._line("")
            self._line("Shutdown...")
            if metrics_path is not None:
                dest = self.export_metrics(metrics_path)
                self._line(f"Metrics JSON........{dest}")
            self._line("Done.")
        elif metrics_path is not None:
            self.export_metrics(metrics_path)
        self._booted = False

    @staticmethod
    def _line(msg: str) -> None:
        print(msg)