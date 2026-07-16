"""ExecutionProfile — capacidades independientes (gptcomputing v0.0.06)."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ExecutionProfile:
    network_enabled: bool = False
    network_permission: bool = False
    battery_optimizer: bool = False
    benchmark_enabled: bool = False
    research_mode: bool = False
    distributed_execution: bool = False
    publish_knowledge: bool = True

    @staticmethod
    def portable() -> ExecutionProfile:
        return ExecutionProfile(
            battery_optimizer=True,
            publish_knowledge=True,
        )

    @staticmethod
    def server() -> ExecutionProfile:
        return ExecutionProfile(
            distributed_execution=True,
            publish_knowledge=True,
        )

    @staticmethod
    def cluster() -> ExecutionProfile:
        return ExecutionProfile(
            network_enabled=True,
            network_permission=True,
            distributed_execution=True,
            publish_knowledge=True,
        )

    @staticmethod
    def laboratory() -> ExecutionProfile:
        return ExecutionProfile(
            research_mode=True,
            benchmark_enabled=True,
            publish_knowledge=False,
        )

    @staticmethod
    def industrial_audit() -> ExecutionProfile:
        return ExecutionProfile(
            distributed_execution=True,
            benchmark_enabled=True,
            publish_knowledge=True,
        )