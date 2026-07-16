"""MDC_FACTOR — factorización toy (solo números ≤10 dígitos)."""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mdc_lib.factoritzacio import MAX_TOY_N, factorizar_mdc_toy
from runtime.plugin import Capability, Cost, IPlugin, PluginContext, PluginInfo, PluginType
from runtime.reference import Reference
from runtime.signature import make_signature


class MdcFactorPlugin(IPlugin):
    PLUGIN_ID = "MDC_FACTOR"

    def info(self) -> PluginInfo:
        return PluginInfo(
            name=self.PLUGIN_ID,
            author="Víctor Manzanares Alberola",
            version=1,
            api_version=1,
            plugin_type=PluginType.VERIFY,
        )

    def capability(self) -> Capability:
        return Capability(deterministic=True, cacheable=True, distributable=True)

    def signature(self, ctx: PluginContext) -> str:
        n = int.from_bytes(ctx.input_payloads[0][:8], "little")
        return make_signature("MDC_FACTOR", [str(n)])

    def estimate(self, ctx: PluginContext) -> Cost:
        n = int.from_bytes(ctx.input_payloads[0][:8], "little")
        return Cost(alu_units=max(1, n.bit_length()))

    def validate(self, ctx: PluginContext) -> bool:
        if len(ctx.input_payloads) != 1 or len(ctx.input_payloads[0]) < 8:
            return False
        n = int.from_bytes(ctx.input_payloads[0][:8], "little")
        return 4 <= n <= MAX_TOY_N

    def execute(self, ctx: PluginContext) -> Reference:
        n = int.from_bytes(ctx.input_payloads[0][:8], "little")
        factor = factorizar_mdc_toy(n)
        result = {"n": n, "factor": factor, "toy_limit": MAX_TOY_N}
        return Reference.create(
            signature=self.signature(ctx),
            payload=json.dumps(result).encode(),
            metadata={"plugin": self.PLUGIN_ID, **result},
        )