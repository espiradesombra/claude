"""
Puente AntiPC → gemelos digitales del repo VMA.

  · ZypyZape  — inercia sintética eólica + viabilidad UDP (arquitectura B)
  · Quijote   — modulación de inercia por masa desplazable en canal
  · Kilòmetre — almacenamiento cinético por flotación (fluidos densos)
"""

from __future__ import annotations

import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

_PKG = Path(__file__).resolve().parents[1]
_REPO_ANTIPC = _PKG.parents[1]
_DESKTOP_REPO = _REPO_ANTIPC.parent


@dataclass
class GemeloScript:
    variant: str
    path: Path
    description: str
    requires: tuple[str, ...] = ("numpy", "matplotlib")


@dataclass
class GemeloInfo:
    key: str
    title: str
    domain: str
    equations: list[str] = field(default_factory=list)
    antipc_link: str = ""
    scripts: list[GemeloScript] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)


def _first_existing(*candidates: Path) -> Path | None:
    for p in candidates:
        if p.is_file():
            return p
    return None


def _build_registry() -> dict[str, GemeloInfo]:
    zyp_viability = _REPO_ANTIPC / "zypyzape" / "viability_udp.py"
    zyp_emitter = _REPO_ANTIPC / "zypyzape" / "hub_emitter.py"
    zyp_explain = _REPO_ANTIPC / "zypyzape" / "EXPLICACION_CIENTIFICA.txt"

    quijote_v48 = _first_existing(
        _DESKTOP_REPO / "quijote" / "zypyzape_twin_v4_8_quijote.py",
        _DESKTOP_REPO / "just run" / "gemelos" / "zypyzape_twin_v4_8_quijote.py",
        _DESKTOP_REPO / "Quijotee" / "zypyzape_twin_v4_8_quijote.py",
    )
    quijote_v5 = _first_existing(
        _DESKTOP_REPO / "quijote" / "gemelo_v5_quijote_control.py",
        _DESKTOP_REPO / "just run" / "gemelos" / "gemelo_v5_quijote_control.py",
        _DESKTOP_REPO / "Quijotee" / "gemelo_v5_quijote_control.py",
    )
    km_base = _DESKTOP_REPO / "kilometre;(soles_bateria)"
    km_v15 = _first_existing(
        km_base / "kilometre_v15_fluids_densos.py",
        km_base / "kilometre_v14_flota.py",
        km_base / "kilometre_dimensions.py",
    )

    scripts_zyp: list[GemeloScript] = []
    if zyp_viability.is_file():
        scripts_zyp.append(
            GemeloScript(
                "viability",
                zyp_viability,
                "Bench UDP A vs B — slot ring + hubs (viabilidad AntiPC)",
            )
        )
    if zyp_emitter.is_file():
        scripts_zyp.append(
            GemeloScript("emitter", zyp_emitter, "Emisor UDP hub simulado (puerto 3333)")
        )

    scripts_q: list[GemeloScript] = []
    if quijote_v48:
        scripts_q.append(
            GemeloScript(
                "v48",
                quijote_v48,
                "Gemelo v4.8 — Quijote masa desplazable + pitch NREL 5MW",
            )
        )
    if quijote_v5:
        scripts_q.append(
            GemeloScript(
                "v5",
                quijote_v5,
                "Gemelo v5 — control activo J(t) con término −ω·J̇",
            )
        )

    scripts_km: list[GemeloScript] = []
    if km_v15:
        scripts_km.append(
            GemeloScript(
                km_v15.stem.replace("kilometre_", ""),
                km_v15,
                "Comparativa fluidos densos — flotación y almacenamiento cinético",
            )
        )

    return {
        "zypyzape": GemeloInfo(
            key="zypyzape",
            title="ZypyZape — gemelo eólico + red AntiPC",
            domain="Energía eólica · inercia sintética · UDP lock-free",
            equations=[
                "J·dω/dt = T_aero − T_gen − T_roz + T_ZZ + T_IS",
                "Swing: 2H·df/dt = P_mec − P_elec − D·Δf",
                "Ring SPSC: coste B ≈ P·(C_kernel + H_ring)  (sin C_user extra)",
            ],
            antipc_link="Arquitectura B — bd_pipeline, game_server, network demo",
            scripts=scripts_zyp,
            outputs=[
                str(_REPO_ANTIPC / "output" / "zypyzape_viability.json"),
                str(zyp_explain) if zyp_explain.is_file() else "",
            ],
        ),
        "quijote": GemeloInfo(
            key="quijote",
            title="Quijote — modulación de inercia por masa en canal",
            domain="Actuador mecánico · batería cinética en palas / canal radial",
            equations=[
                "J_i(r) = J_G + N_b·m_q·r_q[i]²",
                "J·ω̇ = T_net − ω·J̇   (término cinemático Quijote)",
                "F_c = m_q·ω²·r_q  — control P: v_slide = f(Δω, Δf, Δv_viento)",
                "Cp(λ,β) — superficie 2D calibrada NREL 5MW (Jonkman 2009)",
            ],
            antipc_link="Validación dinámica del concepto ZypyZape; no enlazado UDP aún",
            scripts=scripts_q,
            outputs=[
                "zypyzape_v48_quijote.png",
                "gemelo_v5_quijote_control.png",
                "zypyzape_v4_validacion.png",
            ],
        ),
        "kilometre": GemeloInfo(
            key="kilometre",
            title="Kilòmetre — almacenamiento por flotación (soles/batería)",
            domain="Hidráulica densa · tubos km · flotador He/vacío",
            equations=[
                "F_n = (ρ_fluid − ρ_obj)·V_obj·g",
                "E_k = ½·I_tot·ω_max²,  I_tot = I_mec + I_added_mass",
                "P_nom = α·F_n·R·ω_max/2",
                "Comparativa ρ: salmuera, Fe+aceite, Ga, Hg (solo referencia)",
            ],
            antipc_link="Dimensionamiento energético complementario a Quijote/ZypyZape",
            scripts=scripts_km,
            outputs=[
                "kilometre_v15_fluids_densos.png",
                "kilometre_dimensions.png",
            ],
        ),
    }


GEMELO_REGISTRY: dict[str, GemeloInfo] = _build_registry()


def list_gemelos() -> list[GemeloInfo]:
    return list(GEMELO_REGISTRY.values())


def resolve_script(key: str, variant: str | None = None) -> tuple[GemeloInfo, GemeloScript]:
    info = GEMELO_REGISTRY.get(key.lower())
    if info is None:
        raise KeyError(f"gemelo desconocido: {key}")
    if not info.scripts:
        raise FileNotFoundError(f"ningún script disponible para {key}")
    if variant is None:
        return info, info.scripts[0]
    variant_l = variant.lower()
    for s in info.scripts:
        if s.variant.lower() == variant_l:
            return info, s
    raise KeyError(f"variante '{variant}' no encontrada para {key}")


def format_gemelo_info(info: GemeloInfo) -> str:
    lines = [
        f"[{info.key}] {info.title}",
        f"  Dominio     : {info.domain}",
        f"  AntiPC      : {info.antipc_link}",
        "",
        "  Ecuaciones / modelo:",
    ]
    for eq in info.equations:
        lines.append(f"    · {eq}")
    lines.append("")
    lines.append("  Scripts:")
    if info.scripts:
        for s in info.scripts:
            ok = "OK" if s.path.is_file() else "NO"
            lines.append(f"    · {s.variant:10} [{ok}] {s.path}")
            lines.append(f"      {s.description}")
    else:
        lines.append("    (ninguno encontrado en repo)")
    if info.outputs:
        lines.append("")
        lines.append("  Salidas típicas:")
        for o in info.outputs:
            if o:
                lines.append(f"    · {o}")
    return "\n".join(lines)


@dataclass
class GemeloRunResult:
    key: str
    variant: str
    script: Path
    returncode: int
    elapsed_s: float
    stdout: str
    stderr: str


def run_gemelo(
    key: str,
    *,
    variant: str | None = None,
    extra_args: list[str] | None = None,
    timeout_s: float = 300.0,
) -> GemeloRunResult:
    info, script = resolve_script(key, variant)
    if not script.path.is_file():
        raise FileNotFoundError(f"script no encontrado: {script.path}")

    cmd = [sys.executable, str(script.path)]
    if extra_args:
        cmd.extend(extra_args)

    t0 = time.perf_counter()
    proc = subprocess.run(
        cmd,
        cwd=str(script.path.parent),
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    elapsed = time.perf_counter() - t0

    return GemeloRunResult(
        key=info.key,
        variant=script.variant,
        script=script.path,
        returncode=proc.returncode,
        elapsed_s=elapsed,
        stdout=proc.stdout or "",
        stderr=proc.stderr or "",
    )