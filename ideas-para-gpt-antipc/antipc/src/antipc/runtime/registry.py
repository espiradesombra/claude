"""Registry — todo objeto registrado implementa IRegistryObject."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, TypeVar, runtime_checkable

ObjectID = str


@dataclass(frozen=True)
class Descriptor:
    kind: str
    name: str
    version: str = "1"
    capabilities: tuple[str, ...] = ()


@runtime_checkable
class IRegistryObject(Protocol):
    def id(self) -> ObjectID: ...
    def descriptor(self) -> Descriptor: ...


T = TypeVar("T", bound=IRegistryObject)


class Registry:
    """No almacena conocimiento — solo sabe qué existe y dónde."""

    def __init__(self) -> None:
        self._objects: dict[ObjectID, IRegistryObject] = {}

    def register(self, obj: IRegistryObject) -> ObjectID:
        oid = obj.id()
        self._objects[oid] = obj
        return oid

    def find(self, oid: ObjectID) -> IRegistryObject | None:
        return self._objects.get(oid)

    def find_typed(self, oid: ObjectID, kind: str) -> IRegistryObject | None:
        obj = self._objects.get(oid)
        if obj is None or obj.descriptor().kind != kind:
            return None
        return obj

    def find_all(self, kind: str | None = None) -> list[IRegistryObject]:
        items = list(self._objects.values())
        if kind is None:
            return items
        return [o for o in items if o.descriptor().kind == kind]

    def count(self) -> int:
        return len(self._objects)