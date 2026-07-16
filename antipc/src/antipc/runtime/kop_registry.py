"""Registro global de KOP — ids oficiales y validacion de versiones."""

from __future__ import annotations

from .kop import KOPId, KOP_DEFINITIONS, KOPKind


class KOPRegistry:
    """Runtime no necesita conocer todos los KOP — solo los registrados."""

    _handlers: dict[KOPId, int] = {}  # kop_id -> max_version soportada

    @classmethod
    def register(cls, kop_id: KOPId, max_version: int = 1) -> None:
        cls._handlers[kop_id] = max_version

    @classmethod
    def is_supported(cls, kop_id: KOPId, version: int) -> bool:
        max_v = cls._handlers.get(kop_id)
        if max_v is None:
            return kop_id in KOP_DEFINITIONS
        return version <= max_v

    @classmethod
    def list_data_kops(cls) -> list[str]:
        return [
            KOP_DEFINITIONS[k]["name"]
            for k in KOPId
            if KOP_DEFINITIONS[k]["kind"] == KOPKind.DATA
        ]

    @classmethod
    def list_behavior_kops(cls) -> list[str]:
        return [
            KOP_DEFINITIONS[k]["name"]
            for k in KOPId
            if KOP_DEFINITIONS[k]["kind"] == KOPKind.BEHAVIOR
        ]


# KOP oficiales v0.1.1
for kid in KOPId:
    KOPRegistry.register(kid, max_version=1)