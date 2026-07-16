"""Búfer circular de slots fijos 64B — destino directo de recvinto (AntiPC/ZypyZape)."""

from __future__ import annotations

import threading
from dataclasses import dataclass, field


DATA_SIZE = 64


@dataclass
class SlotRing:
    """
    SPSC ring de slots preasignados (power-of-two capacity).
    El hilo de red escribe con recvinto(slot.payload); el worker solo lee por índice.
    """

    capacity: int = 1024
    data_size: int = DATA_SIZE
    _mask: int = field(init=False)
    _slots: list[bytearray] = field(init=False)
    _head: int = 0
    _tail: int = 0
    _lock: threading.Lock = field(default_factory=threading.Lock)

    def __post_init__(self) -> None:
        if self.capacity <= 0 or self.capacity & (self.capacity - 1) != 0:
            raise ValueError("capacity must be power of two")
        self._mask = self.capacity - 1
        self._slots = [bytearray(self.data_size) for _ in range(self.capacity)]

    def write_index(self) -> int | None:
        """Índice del slot libre para recvinto, o None si lleno."""
        with self._lock:
            next_head = (self._head + 1) & self._mask
            if next_head == self._tail:
                return None
            return self._head

    def commit_write(self) -> None:
        with self._lock:
            self._head = (self._head + 1) & self._mask

    def read_index(self) -> int | None:
        with self._lock:
            if self._tail == self._head:
                return None
            return self._tail

    def commit_read(self) -> None:
        with self._lock:
            self._tail = (self._tail + 1) & self._mask

    def slot_view(self, index: int) -> memoryview:
        return memoryview(self._slots[index])

    def slot_bytes(self, index: int) -> bytes:
        return bytes(self._slots[index])

    def __len__(self) -> int:
        with self._lock:
            return (self._head - self._tail) & self._mask