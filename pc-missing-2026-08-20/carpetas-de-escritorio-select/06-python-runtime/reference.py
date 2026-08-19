"""Immutable references and lifecycle states."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any
import uuid


class ReferenceState(Enum):
    WAITING = auto()
    READY = auto()
    PUBLISHED = auto()
    EXPIRED = auto()


@dataclass(frozen=True)
class Reference:
    """Immutable knowledge identifier r ∈ R."""

    ref_id: str
    signature: str
    metadata: dict[str, Any] = field(default_factory=dict)
    payload: bytes | None = None
    parents: tuple[str, ...] = ()

    @staticmethod
    def create(
        signature: str,
        payload: bytes | None = None,
        metadata: dict[str, Any] | None = None,
        parents: tuple[str, ...] = (),
    ) -> Reference:
        return Reference(
            ref_id=str(uuid.uuid4()),
            signature=signature,
            metadata=dict(metadata or {}),
            payload=payload,
            parents=parents,
        )


@dataclass
class ReferenceRecord:
    reference: Reference
    state: ReferenceState = ReferenceState.WAITING