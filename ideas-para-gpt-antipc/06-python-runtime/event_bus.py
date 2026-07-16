"""Event bus with coalescing — Ley de Activación (no polling)."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Callable


class EventType(Enum):
    REFERENCE_CREATED = auto()
    REFERENCE_PUBLISHED = auto()
    OPERATION_READY = auto()
    OPERATION_CANCELLED = auto()


@dataclass
class Event:
    event_type: EventType
    ref_id: str
    signature: str = ""
    epoch: int = 0
    tick: int = 0


@dataclass
class EventBus:
    _queue: deque[Event] = field(default_factory=deque)
    _coalesce: dict[str, Event] = field(default_factory=dict)
    _epoch: int = 1
    _tick: int = 0
    coalesced: int = 0

    def emit(self, event: Event, merge_equal: bool = True) -> None:
        self._tick += 1
        event.epoch = self._epoch
        event.tick = self._tick
        key = f"{event.event_type.name}:{event.ref_id}:{event.signature}"
        if merge_equal and key in self._coalesce:
            self.coalesced += 1
            return
        self._coalesce[key] = event
        self._queue.append(event)

    def drain(self) -> list[Event]:
        events = list(self._queue)
        self._queue.clear()
        self._coalesce.clear()
        return events

    def subscribe_loop(self, handler: Callable[[Event], None], max_events: int = 10_000) -> int:
        processed = 0
        while self._queue and processed < max_events:
            handler(self._queue.popleft())
            processed += 1
        return processed