"""K3_FILE — huella industrial de fichero (tamaño + hash K3). Equivalente k3dedup."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hash_engine import k3_hash_buffer
from runtime.plugin import Capability, Cost, IPlugin, PluginContext, PluginInfo, PluginType
from runtime.reference import Reference
from runtime.signature import make_signature


class K3FilePlugin(IPlugin):
    PLUGIN_ID = "K3_FILE"

    def info(self) -> PluginInfo:
        return PluginInfo(
            name=self.PLUGIN_ID,
            author="Víctor Manzanares Alberola",
            version=1,
            api_version=1,
            plugin_type=PluginType.ACQUIRE,
        )

    def capability(self) -> Capability:
        return Capability(deterministic=True, distributable=True, cacheable=True, streamable=True)

    def signature(self, ctx: PluginContext) -> str:
        path = ctx.input_payloads[0].decode("utf-8", errors="replace")
        try:
            size = os.path.getsize(path)
        except OSError:
            size = -1
        return make_signature("K3_FILE", [path, str(size)])

    def estimate(self, ctx: PluginContext) -> Cost:
        path = ctx.input_payloads[0].decode("utf-8", errors="replace")
        try:
            size = os.path.getsize(path)
        except OSError:
            size = 1024
        return Cost(alu_units=max(1, size // 65536), memory_bytes=size)

    def validate(self, ctx: PluginContext) -> bool:
        if len(ctx.input_payloads) != 1:
            return False
        path = ctx.input_payloads[0].decode("utf-8", errors="replace")
        return os.path.isfile(path)

    def execute(self, ctx: PluginContext) -> Reference:
        path = ctx.input_payloads[0].decode("utf-8", errors="replace")
        size = os.path.getsize(path)
        with open(path, "rb") as f:
            data = f.read()
        digest = k3_hash_buffer(data)
        sig = self.signature(ctx)
        return Reference.create(
            signature=sig,
            payload=digest.to_bytes(4, "little"),
            metadata={
                "plugin": self.PLUGIN_ID,
                "path": path,
                "size": size,
                "digest_hex": f"{digest:08X}",
            },
            parents=tuple(r.ref_id for r in ctx.input_refs),
        )