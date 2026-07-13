"""Universal Plugin System (UPS) — IPlugin contract."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto

from .reference import Reference


class PluginType(Enum):
    ACQUIRE = auto()
    TRANSFORM = auto()
    INFER = auto()
    VERIFY = auto()
    PUBLISH = auto()
    EFFECT = auto()


@dataclass
class PluginInfo:
    name: str
    author: str
    version: int
    api_version: int
    plugin_type: PluginType


@dataclass
class Capability:
    deterministic: bool = True
    distributable: bool = False
    cacheable: bool = True
    reversible: bool = False
    streamable: bool = False
    requires_network: bool = False
    requires_permission: bool = False


@dataclass
class Cost:
    alu_units: int = 1
    memory_bytes: int = 0


@dataclass
class PluginContext:
    input_refs: list[Reference]
    input_payloads: list[bytes]


class IPlugin(ABC):
    @abstractmethod
    def info(self) -> PluginInfo: ...

    @abstractmethod
    def signature(self, ctx: PluginContext) -> str: ...

    @abstractmethod
    def capability(self) -> Capability: ...

    @abstractmethod
    def estimate(self, ctx: PluginContext) -> Cost: ...

    @abstractmethod
    def validate(self, ctx: PluginContext) -> bool: ...

    @abstractmethod
    def execute(self, ctx: PluginContext) -> Reference: ...