"""Knowledge Resolver v0.0.064 — resolve(signature), no exists()."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto

from .knowledge import KnowledgeBuffer
from .profile import ExecutionProfile
from .reference import Reference


class ResolutionStatus(Enum):
    FOUND_LOCAL = auto()
    FOUND_SHARED = auto()
    FOUND_REMOTE = auto()
    EXECUTE = auto()
    WAIT = auto()
    FAILED = auto()


@dataclass
class Resolution:
    status: ResolutionStatus
    reference: Reference | None
    source: str
    confidence: float
    estimated_cost: int


@dataclass
class KnowledgeResolver:
    knowledge: KnowledgeBuffer
    profile: ExecutionProfile
    shared_cache: dict[str, Reference] | None = None

    def resolve(self, signature: str, estimated_cost: int = 1) -> Resolution:
        if self.profile.research_mode:
            return Resolution(
                ResolutionStatus.EXECUTE,
                None,
                "research:force_execute",
                1.0,
                estimated_cost,
            )

        local = self.knowledge.lookup(signature)
        if local is not None:
            return Resolution(
                ResolutionStatus.FOUND_LOCAL,
                local,
                "knowledge_buffer:local",
                1.0,
                0,
            )

        if self.shared_cache and signature in self.shared_cache:
            ref = self.shared_cache[signature]
            return Resolution(
                ResolutionStatus.FOUND_SHARED,
                ref,
                "knowledge_buffer:shared",
                0.95,
                0,
            )

        if self.profile.network_enabled and self.profile.network_permission:
            return Resolution(
                ResolutionStatus.WAIT,
                None,
                "network:pending_remote",
                0.5,
                estimated_cost,
            )

        return Resolution(
            ResolutionStatus.EXECUTE,
            None,
            "alu:execute",
            1.0,
            estimated_cost,
        )