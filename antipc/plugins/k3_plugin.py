"""K3 HASH plugin — first UPS module."""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hash_engine import k3_hash_buffer
from runtime.plugin import (
    Capability,
    Cost,
    IPlugin,
    PluginContext,
    PluginInfo,
    PluginType,
)
from runtime.reference import Reference
from runtime.signature import make_signature


class K3HashPlugin(IPlugin):
    def info(self) -> PluginInfo:
        return PluginInfo(
            name="K3_HASH",
            author="Víctor Manzanares Alberola",
            version=1,
            api_version=1,
            plugin_type=PluginType.TRANSFORM,
        )

    def capability(self) -> Capability:
        return Capability(
            deterministic=True,
            distributable=True,
            cacheable=True,
            streamable=True,
        )

    def signature(self, ctx: PluginContext) -> str:
        payload_key = ctx.input_payloads[0].hex() if ctx.input_payloads else ""
        return make_signature("K3_HASH", [payload_key])

    def estimate(self, ctx: PluginContext) -> Cost:
        size = sum(len(p) for p in ctx.input_payloads)
        return Cost(alu_units=max(1, size // 32), memory_bytes=size)

    def validate(self, ctx: PluginContext) -> bool:
        return len(ctx.input_payloads) == 1 and len(ctx.input_payloads[0]) > 0

    def execute(self, ctx: PluginContext) -> Reference:
        payload = ctx.input_payloads[0]
        digest = k3_hash_buffer(payload)
        sig = self.signature(ctx)
        parents = tuple(r.ref_id for r in ctx.input_refs)
        return Reference.create(
            signature=sig,
            payload=digest.to_bytes(4, "little"),
            metadata={"plugin": "K3_HASH", "digest_hex": f"{digest:08X}"},
            parents=parents,
        )