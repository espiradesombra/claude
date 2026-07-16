"""Tipos compartidos del microkernel — sin imports circulares."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any


class KOPStatus(Enum):
    OK = auto()
    SKIPPED = auto()
    DENIED = auto()
    ERROR = auto()


@dataclass(frozen=True)
class KOPCost:
    cpu_units: float = 1.0
    ram_bytes: int = 0
    network_units: float = 0.0
    storage_bytes: int = 0
    trust_gain: float = 0.0


@dataclass
class KOPResult:
    status: KOPStatus
    kop_id: int
    message: str = ""
    data: Any = None
    cost: KOPCost = field(default_factory=KOPCost)


@dataclass
class RuntimeContext:
    registry: Any
    scheduler: Any
    events: Any
    knowledge: Any
    ledger: Any
    policy: Any
    metrics: Any