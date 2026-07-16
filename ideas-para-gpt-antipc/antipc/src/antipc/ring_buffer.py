"""Single-producer / single-consumer lock-free ring buffer (SPSC)."""

from __future__ import annotations

from typing import Generic, Optional, TypeVar

T = TypeVar("T")


class SPSCRingBuffer(Generic[T]):
    """
    Lock-free circular buffer for one writer and one reader.

    Uses memory barriers via Python's GIL-safe atomic int operations on indices.
    Suitable for decoupling network acquisition from ALU-bound workers.
    """

    __slots__ = ("_capacity", "_mask", "_slots", "_head", "_tail")

    def __init__(self, capacity: int) -> None:
        if capacity <= 0 or capacity & (capacity - 1) != 0:
            raise ValueError("capacity must be a positive power of two")
        self._capacity = capacity
        self._mask = capacity - 1
        self._slots: list[Optional[T]] = [None] * capacity
        self._head = 0
        self._tail = 0

    def push(self, item: T) -> bool:
        next_head = (self._head + 1) & self._mask
        if next_head == self._tail:
            return False
        self._slots[self._head] = item
        self._head = next_head
        return True

    def pop(self) -> Optional[T]:
        if self._tail == self._head:
            return None
        item = self._slots[self._tail]
        self._slots[self._tail] = None
        self._tail = (self._tail + 1) & self._mask
        return item

    def __len__(self) -> int:
        return (self._head - self._tail) & self._mask