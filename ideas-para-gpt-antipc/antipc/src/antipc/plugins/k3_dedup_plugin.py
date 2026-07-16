"""K3_DEDUP — agrupa huellas por (tamaño, hash) — lógica industrial k3dedup."""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from runtime.plugin import Capability, Cost, IPlugin, PluginContext, PluginInfo, PluginType
from runtime.reference import Reference
from runtime.signature import make_signature


class K3DedupPlugin(IPlugin):
    PLUGIN_ID = "K3_DEDUP"

    def info(self) -> PluginInfo:
        return PluginInfo(
            name=self.PLUGIN_ID,
            author="Víctor Manzanares Alberola",
            version=1,
            api_version=1,
            plugin_type=PluginType.VERIFY,
        )

    def capability(self) -> Capability:
        return Capability(deterministic=True, cacheable=True)

    def signature(self, ctx: PluginContext) -> str:
        keys = []
        for ref in ctx.input_refs:
            keys.append(f"{ref.metadata.get('size')}:{ref.metadata.get('digest_hex')}")
        return make_signature("K3_DEDUP", sorted(keys))

    def estimate(self, ctx: PluginContext) -> Cost:
        return Cost(alu_units=len(ctx.input_refs))

    def validate(self, ctx: PluginContext) -> bool:
        return len(ctx.input_refs) >= 1 and all(
            "digest_hex" in r.metadata for r in ctx.input_refs
        )

    def execute(self, ctx: PluginContext) -> Reference:
        groups: dict[str, list[str]] = {}
        for ref in ctx.input_refs:
            key = f"{ref.metadata['size']}:{ref.metadata['digest_hex']}"
            groups.setdefault(key, []).append(ref.metadata.get("path", "?"))

        duplicates = {k: v for k, v in groups.items() if len(v) > 1}
        recoverable = sum(
            int(k.split(":")[0]) * (len(v) - 1)
            for k, v in duplicates.items()
        )

        report = json.dumps({
            "groups": len(groups),
            "duplicate_groups": len(duplicates),
            "bytes_recoverable": recoverable,
            "duplicates": duplicates,
        }, ensure_ascii=False).encode("utf-8")

        return Reference.create(
            signature=self.signature(ctx),
            payload=report,
            metadata={
                "plugin": self.PLUGIN_ID,
                "duplicate_groups": len(duplicates),
                "bytes_recoverable": recoverable,
            },
            parents=tuple(r.ref_id for r in ctx.input_refs),
        )