"""Scheduler: ¿Qué está listo? — Ready Counter + Execution Queue."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .plugin import IPlugin, PluginContext


@dataclass
class PendingOperation:
    plugin: IPlugin
    ctx: PluginContext
    signature: str
    required: int
    available: int = 0
    priority: float = 0.0
    cancelled: bool = False

    @property
    def ready(self) -> bool:
        return not self.cancelled and self.available >= self.required


@dataclass
class Scheduler:
    pending: list[PendingOperation] = field(default_factory=list)
    executed: int = 0
    cancelled: int = 0

    def submit(self, op: PendingOperation) -> None:
        self.pending.append(op)
        self.pending.sort(key=lambda o: o.priority, reverse=True)

    def ready_ops(self) -> list[PendingOperation]:
        return [op for op in self.pending if op.ready]

    def pop_ready(self) -> PendingOperation | None:
        for i, op in enumerate(self.pending):
            if op.ready:
                return self.pending.pop(i)
        return None

    def cancel_if_reusable(self, signature: str) -> bool:
        for op in self.pending:
            if op.signature == signature and not op.ready:
                op.cancelled = True
                self.cancelled += 1
                return True
        return False