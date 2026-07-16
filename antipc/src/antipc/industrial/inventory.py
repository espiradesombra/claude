"""
Inventario industrial AntiPC — qué hace, qué tiene, estado operativo.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_REPO = _ROOT.parent.parent


@dataclass
class InventorySection:
    title: str
    items: list[str] = field(default_factory=list)


@dataclass
class AntiPCInventory:
    generated_at: str
    cli_version: str
    native_version: str | None
    native_backend: str
    dll_present: bool
    dll_path: str | None
    sections: list[InventorySection] = field(default_factory=list)
    health: dict = field(default_factory=dict)
    industrial_readiness: list[str] = field(default_factory=list)
    gaps: list[str] = field(default_factory=list)


def _native_status() -> tuple[str, str | None, bool, str | None]:
    try:
        from native_engine import get_backend, load_native, status_report

        load_native()
        backend = get_backend()
        dll_path = None
        native_ver = None
        lib = load_native()
        if lib and hasattr(lib, "antipc_native_version"):
            native_ver = lib.antipc_native_version().decode("ascii", errors="replace")
        for line in status_report().splitlines():
            if line.startswith("DLL path"):
                dll_path = line.split(":", 1)[1].strip()
                if dll_path == "—":
                    dll_path = None
        present = dll_path is not None and Path(dll_path).is_file()
        return backend, native_ver, present, dll_path
    except Exception as exc:
        return f"error:{exc}", None, False, None


def _health_checks() -> dict:
    checks: dict[str, str] = {}
    t0 = time.perf_counter()
    try:
        from hash_engine import k3_hash_buffer
        from k3_hash import k3_hash_buffer as py_hash

        sample = b"AntiPC-industrial-audit"
        h1 = k3_hash_buffer(sample)
        h2 = py_hash(sample)
        checks["k3_parity"] = "OK" if h1 == h2 else f"MISMATCH {h1:08X}!={h2:08X}"
    except Exception as exc:
        checks["k3_parity"] = f"FAIL:{exc}"

    try:
        from native_engine import mdc_factor, sieve_count

        f = mdc_factor(10473029)
        s = sieve_count(1000)
        checks["mdc_factor_10473029"] = str(f)
        checks["sieve_1000"] = str(s)
        checks["native_alu_ms"] = str(int((time.perf_counter() - t0) * 1000))
    except Exception as exc:
        checks["native_alu"] = f"FAIL:{exc}"

    try:
        from architecture.network_auth import NetworkAuthGate

        gate = NetworkAuthGate()
        uid, nonce = gate.issue_challenge("audit")
        resp = gate.client_response(uid, nonce)
        checks["hmac_auth"] = "OK" if gate.verify_response(uid, resp) else "FAIL"
    except Exception as exc:
        checks["hmac_auth"] = f"FAIL:{exc}"

    plugins_ok = 0
    for name in ("K3_HASH", "K3_FILE", "K3_DEDUP", "MDC_FACTOR", "MDC_PHASE", "MDC_REGLA"):
        try:
            from runtime.plugin_manager import PluginManager
            from plugins.k3_plugin import K3HashPlugin
            from plugins.k3_file_plugin import K3FilePlugin
            from plugins.k3_dedup_plugin import K3DedupPlugin
            from plugins.mdc_factor_plugin import MdcFactorPlugin
            from plugins.mdc_phase_plugin import MdcPhasePlugin
            from plugins.mdc_regla_plugin import MdcReglaPlugin

            mapping = {
                "K3_HASH": K3HashPlugin,
                "K3_FILE": K3FilePlugin,
                "K3_DEDUP": K3DedupPlugin,
                "MDC_FACTOR": MdcFactorPlugin,
                "MDC_PHASE": MdcPhasePlugin,
                "MDC_REGLA": MdcReglaPlugin,
            }
            p = mapping[name]()
            if p.info().name == name:
                plugins_ok += 1
        except Exception:
            pass
    checks["plugins_loaded"] = f"{plugins_ok}/6"

    return checks


def build_inventory(cli_version: str = "0.10.0-cmd") -> AntiPCInventory:
    backend, native_ver, dll_ok, dll_path = _native_status()
    health = _health_checks()

    sections = [
        InventorySection(
            "Núcleo y fórmula",
            [
                "P_util(N) = N·E(N) + K(N) — comunicación → cómputo reutilizable",
                "FlowKernel + KnowledgeResolver + Reference inmutable",
                "Microkernel KOPs: Boot, Hash, Identity",
                "ExecutionProfile: portable | server | cluster | laboratory | industrial_audit",
            ],
        ),
        InventorySection(
            "Arquitecturas A–E (benchmark)",
            [
                "A — Convencional: memcpy + cola bloqueante",
                "B — Lock-free / SlotRing recvinto (ZypyZape)",
                "C — Distribuido: partición hubs + cache K(N)",
                "D — Grafcet: lotes + matriz existencia + redundancia K3",
                "E — UDP real: hub_node.py + WORK/RESULT + HMAC",
            ],
        ),
        InventorySection(
            "Plugins industriales (IPlugin)",
            [
                "K3_HASH — digest buffer determinista",
                "K3_FILE — huella fichero (tamaño + hash streaming C)",
                "K3_DEDUP — grupos duplicados (tamaño, hash)",
                "MDC_FACTOR — factorización toy ≤10 dígitos",
                "MDC_PHASE — fase/curvatura modular",
                "MDC_REGLA — regla mecánica escala log",
            ],
        ),
        InventorySection(
            "antipc_native.dll (C v0.7.0-c)",
            [
                "k3hash + k3_suite: file, fingerprint, redundant, heavy, hamming",
                "MDC: factor, scan_trains, ksweep classic/predict",
                "Cribas: Eratóstenes, híbrida VMA, modular 6k±1, Criva π(x)/x",
                "Newton rápido + oráculo MEcuation",
                "Geo: converge binaria, Aleatorovix, masivo crypt π/e",
                "K3 stream XOR acordeón (33×1)",
            ],
        ),
        InventorySection(
            "CLI antipc (comandos)",
            [
                "version | boot | mk hash",
                "hash | wave | mechanical",
                "mdc: factor | analyze | ksweep | visual",
                "native: status | bench",
                "criba | newton | geo | geo-masivo",
                "k3: hash | file | verify | dedup | simil | search | redundant | heavy | hamming",
                "network demo | network bench — B slot-ring → D Grafcet",
                "game demo | game cluster | game bots — MMO masivo",
                "gemelo list | info | run — ZypyZape, Quijote, Kilòmetre",
            ],
        ),
        InventorySection(
            "Red y MMO",
            [
                "UdpFabric L3: spawn hubs, AUTH HMAC, WORK/RESULT",
                "hub_node.py — heavy K3 remoto + cache",
                "bd_pipeline — ingest UDP puerto 3333",
                "game_protocol 64B — estado jugador MM",
                "game cluster — heavy por shard; --hub-hosts para LAN remota",
                "AntiPCStack — capas L0..L4 unificadas",
            ],
        ),
        InventorySection(
            "Pipelines industriales",
            [
                "Auditoría ficheros: K3_FILE → resolver → K3_DEDUP + telemetría CSV/JSON",
                "industrial_demo.py — monolítico vs industrial vs 2ª pasada K(N)",
                "demo_stack_udp.py — hash remoto en hubs vs local",
                "Telemetría: IndustrialTelemetry → CSV + informe JSON",
            ],
        ),
        InventorySection(
            "Scripts Windows (scripts/)",
            [
                "01_instalar | 02_benchmark | 03_benchmark_udp | 04_udp_solo",
                "05_verificar_k3 | 08_build_k3hash | 21_build_antipc_native",
                "09_demo_industrial | 16_antipc_cmd | 17_build_antipc_exe",
                "20_network_demo | 22_game_mmo_demo | 23_game_cluster_demo",
            ],
        ),
        InventorySection(
            "Libros 1–6 — métodos VMA (H8)",
            [
                "L1 Números i numeritos — criba 6k±1, desmemoriada, salto √n",
                "L2 Números oTra VeZ — dos primos, MRAUV, híbrida",
                "L3 Sigo en mis Trece — Sofí, Criva, pitagórico, K=9/24",
                "L4 Encriptación/E. Verde — geo, gemelos, K3 stream",
                "L5 (2v+3) MDC — factor, analyze, ksweep, jerk",
                "L6 Implicaciones — newton, método-V*, teoremas 21–26",
                "CLI: antipc libro list | metodos | info | run",
                "criba --desmemoriada · mdc jerk",
            ],
        ),
        InventorySection(
            "DeepSeek 6 2026 + teoremas (H9 — v0.14)",
            [
                "deepseekjun26/ — tablas Libros 1–4, chats jun 2026",
                "filestot l5/ — mdc_v15–v23, ksweep_predictiu, benchmark",
                "PY L5/ — deepseek_python_20260328_*.py",
                "deepseek/bridge.py — registro scripts + docs",
                "teoremas/bridge.py — INDICE_MAESTRO 01–30",
                "mdc_lib/mdc_jerk.py — pinça 4+4 invariante jerk",
                "CLI: antipc deepseek list | info | run",
                "CLI: antipc teorema list | info",
            ],
        ),
        InventorySection(
            "Gemelos digitales (H7 — repo VMA)",
            [
                "ZypyZape — J·dω/dt, swing, slot_ring → arquitectura B AntiPC",
                "zypyzape/viability_udp.py — bench UDP A vs B",
                "Quijote — J_i(r), ω·J̇, Cp(λ,β) NREL 5MW — repo/quijote/",
                "gemelo v4.8 + v5 control activo inercia",
                "Kilòmetre — flotación F_n, E_k cinética — repo/kilometre;(soles_bateria)/",
                "CLI: antipc gemelo list | info | run",
            ],
        ),
        InventorySection(
            "Integraciones Desktop VMA",
            [
                "vma-methods — cribas, Newton, Criva",
                "vma-k3 / HASHTOOLCODE — k3hash.c",
                "encriptacionGeometrica — Aleatorovix, convergencia",
                "mdc_lib — analyze, visual, factoritzacio",
                "antipc-port-c — logs integración v04–v11",
            ],
        ),
    ]

    readiness = [
        "Hash K3 auditable (C + Python paridad)",
        "Reuse knowledge buffer (2ª pasada sin ALU)",
        "Telemetría CSV/JSON exportable",
        "Red con permiso HMAC obligatorio",
        "Plugins deterministas y firmables",
        "DLL nativa unificada para ALU crítica",
        "Demos reproducibles por .bat",
        "MMO cluster con hubs por shard",
    ]

    gaps = [
        "README/LEEME revisar tras cada release mayor (actual: v0.14)",
        "Criba desmemoriada solo Python (vma-methods) — sin port C aún",
        "MDC jerk/mdc_v23 factorización toy — no RSA industrial",
        "Sin TLS/HTTPS — UDP laboratorio; no producción WAN",
        "K3 no criptográfico — no contraseñas ni firma legal",
        "MDC factor toy ≤10 dígitos — no RSA real",
        "Sin orquestador Kubernetes/Docker oficial",
        "Tests automatizados parciales (microkernel, wave)",
        "antipc.exe PyInstaller — verificar DLL empaquetada en deploy",
    ]

    if health.get("k3_parity") != "OK":
        gaps.insert(0, "ALERTA: paridad K3 engine vs Python rota")
    if not dll_ok:
        gaps.insert(0, "ALERTA: antipc_native.dll no cargada")

    return AntiPCInventory(
        generated_at=time.strftime("%Y-%m-%dT%H:%M:%S"),
        cli_version=cli_version,
        native_version=native_ver,
        native_backend=backend,
        dll_present=dll_ok,
        dll_path=dll_path,
        sections=sections,
        health=health,
        industrial_readiness=readiness,
        gaps=gaps,
    )


def format_inventory_text(inv: AntiPCInventory) -> str:
    lines = [
        "=" * 72,
        "  INVENTARIO INDUSTRIAL AntiPC",
        f"  Generado: {inv.generated_at}",
        "=" * 72,
        "",
        f"  CLI          : {inv.cli_version}",
        f"  Native       : {inv.native_version or '—'}",
        f"  Backend      : {inv.native_backend}",
        f"  DLL          : {'OK' if inv.dll_present else 'NO'} {inv.dll_path or ''}",
        "",
        "  SALUD OPERATIVA",
        "-" * 72,
    ]
    for k, v in inv.health.items():
        lines.append(f"    {k:24} {v}")
    lines.extend(["", "  CAPACIDADES", "-" * 72])
    for sec in inv.sections:
        lines.append(f"\n  [{sec.title}]")
        for item in sec.items:
            lines.append(f"    · {item}")
    lines.extend(["", "  LISTO PARA INDUSTRIA", "-" * 72])
    for r in inv.industrial_readiness:
        lines.append(f"    ✓ {r}")
    lines.extend(["", "  BRECHAS / ROADMAP", "-" * 72])
    for g in inv.gaps:
        lines.append(f"    △ {g}")
    lines.extend(["", "=" * 72, "  VMA — AntiPC / 33×1", "=" * 72])
    return "\n".join(lines)


def export_inventory(
    *,
    cli_version: str = "0.9.0-cmd",
    txt_path: Path | None = None,
    json_path: Path | None = None,
) -> AntiPCInventory:
    inv = build_inventory(cli_version=cli_version)
    text = format_inventory_text(inv)

    if txt_path:
        txt_path.parent.mkdir(parents=True, exist_ok=True)
        txt_path.write_text(text, encoding="utf-8")

    if json_path:
        json_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "generated_at": inv.generated_at,
            "cli_version": inv.cli_version,
            "native_version": inv.native_version,
            "native_backend": inv.native_backend,
            "dll_present": inv.dll_present,
            "dll_path": inv.dll_path,
            "health": inv.health,
            "sections": {s.title: s.items for s in inv.sections},
            "industrial_readiness": inv.industrial_readiness,
            "gaps": inv.gaps,
        }
        json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    return inv


def main(argv: list[str] | None = None) -> int:
    out_txt = _ROOT / "output" / "inventario_industrial.txt"
    out_json = _ROOT / "output" / "inventario_industrial.json"
    if argv and len(argv) > 1:
        out_txt = Path(argv[1])
    inv = export_inventory(txt_path=out_txt, json_path=out_json)
    print(format_inventory_text(inv))
    print(f"\n  Exportado: {out_txt}")
    print(f"  Exportado: {out_json}")
    return 0 if inv.dll_present and inv.health.get("k3_parity") == "OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())