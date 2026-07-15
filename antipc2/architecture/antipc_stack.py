"""
AntiPC Stack — arquitectura unificada HASHTOOLCODE + MDC + Runtime + Red.

Capas:
  L0  Fundación     K3 hash, HMAC auth (HASHTOOLCODE)
  L1  Dominio MDC   Regla Mecánica, fase, factorización toy
  L2  Flow Kernel   Referencias, resolver, event bus, plugins
  L3  Red           UDP hubs, permisos
  L4  Industrial    Telemetría, perfiles, auditoría
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path

from architecture.network_auth import NetworkAuthGate
from architecture.udp_fabric import UdpFabric
from hash_engine import get_backend, source_origin
from industrial.runtime import IndustrialRuntime
from plugins.k3_dedup_plugin import K3DedupPlugin
from plugins.k3_file_plugin import K3FilePlugin
from plugins.k3_plugin import K3HashPlugin
from plugins.mdc_factor_plugin import MdcFactorPlugin
from plugins.mdc_phase_plugin import MdcPhasePlugin
from plugins.mdc_regla_plugin import MdcReglaPlugin
from runtime.profile import ExecutionProfile
from runtime.reference import Reference


class Layer(Enum):
    L0_FOUNDATION = auto()
    L1_MDC = auto()
    L2_RUNTIME = auto()
    L3_NETWORK = auto()
    L4_INDUSTRIAL = auto()


@dataclass
class LayerInfo:
    layer: Layer
    name: str
    components: list[str]
    source: str


@dataclass
class AntiPCStack:
    """Arquitectura completa AntiPC — punto de entrada único."""

    profile: ExecutionProfile = field(default_factory=ExecutionProfile.cluster)
    runtime: IndustrialRuntime = field(default_factory=IndustrialRuntime)
    auth: NetworkAuthGate = field(default_factory=NetworkAuthGate)
    user_id: str = "antipc-node-1"
    network_granted: bool = False
    fabric: UdpFabric | None = None
    _hub_count: int = 0

    def __post_init__(self) -> None:
        self.runtime.profile = self.profile
        self.bootstrap()
        self.fabric = UdpFabric(auth=self.auth, user_id=self.user_id)

    def bootstrap(self) -> None:
        plugins = [
            K3HashPlugin(),
            K3FilePlugin(),
            K3DedupPlugin(),
            MdcPhasePlugin(),
            MdcReglaPlugin(),
            MdcFactorPlugin(),
        ]
        self.runtime.bootstrap(*plugins)

    def architecture_map(self) -> list[LayerInfo]:
        return [
            LayerInfo(Layer.L0_FOUNDATION, "Fundación HASHTOOLCODE",
                      ["K3_HASH", "hash_engine", "HMAC NetworkAuthGate"],
                      str(source_origin())),
            LayerInfo(Layer.L1_MDC, "Regla Mecánica / MDC",
                      ["MDC_PHASE", "MDC_REGLA", "MDC_FACTOR", "ReglaCalculo"],
                      "HASHTOOLCODE/cpp y dll/biblio"),
            LayerInfo(Layer.L2_RUNTIME, "Flow Kernel",
                      ["References", "KnowledgeResolver", "EventBus", "PluginManager"],
                      "gptcomputing.txt v0.1.0"),
            LayerInfo(Layer.L3_NETWORK, "Red distribuida",
                      ["UDP hubs", "WORK/RESULT protocol", "permisos HMAC"],
                      "udp_benchmark.py / hub_node.py"),
            LayerInfo(Layer.L4_INDUSTRIAL, "Industrial",
                      ["Telemetry CSV", "ExecutionProfile", "K3_FILE dedup"],
                      "industrial_demo.py"),
        ]

    def request_network_permission(self) -> bool:
        """Handshake HMAC local (precondición antes de abrir red)."""
        uid, nonce = self.auth.issue_challenge(self.user_id)
        response = self.auth.client_response(uid, nonce)
        ok = self.auth.verify_response(uid, response)
        if ok:
            self.network_granted = True
            self.runtime.kernel.network_permission = True
            self.runtime.profile.network_permission = True
        return ok

    def start_network(self, hubs: int = 4) -> dict:
        """
        L3: arranca hubs UDP, autentica cada uno con HMAC, enlaza al stack.
        Requiere request_network_permission() previo.
        """
        if not self.network_granted:
            raise PermissionError("red no autorizada: ejecuta request_network_permission()")
        if self.fabric is None:
            self.fabric = UdpFabric(auth=self.auth, user_id=self.user_id)

        self.fabric.start_hubs(hubs)
        self._hub_count = hubs
        authed = self.fabric.authenticate_all()
        self.runtime.telemetry.log(
            "network_start", "UdpFabric", f"hubs={hubs}", f"auth={authed}/{hubs}"
        )
        return {"hubs": hubs, "authenticated": authed, "total": hubs}

    def dispatch_remote_hash(self, payloads: list[bytes]) -> dict[int, int]:
        """Despacha HASH a hubs autenticados; resultado vía UDP."""
        if not self.fabric or not self.network_granted:
            raise PermissionError("red no activa")
        results = self.fabric.dispatch_hashes(payloads)
        for seq, digest in results.items():
            self.runtime.telemetry.log(
                "remote_hash", "K3_HASH", f"seq={seq}", "FOUND_REMOTE",
                alu_saved=True,
            )
        return results

    def stop_network(self) -> None:
        if self.fabric:
            self.fabric.shutdown()
        self._hub_count = 0

    def pipeline_remote_files(self, paths: list[str]) -> list[dict]:
        """Lee ficheros localmente, envía bytes a hubs para HASH remoto."""
        payloads: list[bytes] = []
        meta: list[dict] = []
        for path in paths:
            try:
                with open(path, "rb") as f:
                    data = f.read()
                payloads.append(data)
                meta.append({"path": path, "size": len(data)})
            except OSError:
                continue

        if not payloads:
            return []

        digests = self.dispatch_remote_hash(payloads)
        out = []
        for i, m in enumerate(meta):
            d = digests.get(i, 0)
            m["digest_hex"] = f"{d:08X}"
            m["source"] = "L3_UDP_HUB"
            out.append(m)
        return out

    def pipeline_audit_file(self, path: str) -> dict:
        """Pipeline industrial: fichero → K3_FILE → K3_HASH → MDC_PHASE."""
        fp = self.runtime.fingerprint_file(path)
        if fp is None:
            return {"error": "K3_FILE failed", "path": path}

        size = fp.metadata.get("size", 0)
        nm = (size.to_bytes(8, "little") + (size % 997 + 1).to_bytes(8, "little"))
        phase_ref = self.runtime._run_plugin("MDC_PHASE", [fp], nm)
        return {
            "path": path,
            "k3_digest": fp.metadata.get("digest_hex"),
            "size": size,
            "mdc_fase": phase_ref.metadata.get("fase") if phase_ref else None,
        }

    def pipeline_analyze_number(self, n: int) -> dict:
        """Pipeline MDC: número → MDC_FACTOR + MDC_REGLA (toy)."""
        payload = n.to_bytes(8, "little") + (n % 1000 + 1).to_bytes(8, "little")
        empty = self.runtime.kernel.ingest_raw(b"", label="SEED")
        factor_ref = self.runtime._run_plugin("MDC_FACTOR", [empty], payload)
        regla_ref = self.runtime._run_plugin("MDC_REGLA", [empty], payload)
        return {
            "n": n,
            "factor": factor_ref.metadata.get("factor") if factor_ref else None,
            "regla_producto": regla_ref.metadata.get("producto_regla") if regla_ref else None,
        }

    def export_architecture_report(self, path: Path) -> None:
        lines = [
            "=" * 72,
            "  AntiPC — ARQUITECTURA UNIFICADA",
            "  HASHTOOLCODE + MDC + Flow Kernel + Red + Industrial",
            "=" * 72,
            "",
            f"  Hash backend: {get_backend()}",
            f"  Red autorizada: {self.network_granted}",
            "",
        ]
        for layer in self.architecture_map():
            lines.append(f"  [{layer.layer.name}] {layer.name}")
            lines.append(f"    Fuente: {layer.source}")
            for c in layer.components:
                lines.append(f"      - {c}")
            lines.append("")
        lines.append("  Plugins registrados:")
        for name in self.runtime.manager.plugins:
            info = self.runtime.manager.plugins[name].info()
            lines.append(f"      {name} v{info.version} ({info.plugin_type.name})")
        lines.append("=" * 72)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines), encoding="utf-8")