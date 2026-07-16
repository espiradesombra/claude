"""
Corpus DeepSeek 6 2026 → comandos AntiPC.

Fuentes:
  repo/deepseekjun26/          — chats y tablas Libros 1–4 (jun 2026)
  repo/filestot l5/            — MDC v15–v23, ksweep, benchmark
  repo/PY L5/                  — deepseek_python_20260328_*.py
  repo/cosas/deepseek 6-5-26.txt
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

_DEEPSEEK_DIR = _DESKTOP_REPO / "deepseekjun26"
_L5_SCRIPTS = _DESKTOP_REPO / "filestot l5"
_PY_L5 = _DESKTOP_REPO / "PY L5"


@dataclass
class DeepSeekScript:
    id: str
    path: Path
    description: str
    libro: int | None = None
    antipc_cmd: str = ""


@dataclass
class DeepSeekEntry:
    key: str
    titulo: str
    tema: str
    fuente: str
    ecuaciones: list[str] = field(default_factory=list)
    scripts: list[DeepSeekScript] = field(default_factory=list)
    docs: list[Path] = field(default_factory=list)


def _existing(path: Path) -> Path | None:
    return path if path.is_file() else None


def _collect_txt_docs() -> list[Path]:
    if not _DEEPSEEK_DIR.is_dir():
        return []
    return sorted(_DEEPSEEK_DIR.glob("*.txt"))


def _script(id_: str, rel: Path, desc: str, *, libro: int | None = None, cmd: str = "") -> DeepSeekScript | None:
    for base in (_L5_SCRIPTS, _PY_L5, _DESKTOP_REPO):
        p = base / rel
        if p.is_file():
            return DeepSeekScript(id_, p, desc, libro=libro, antipc_cmd=cmd)
    return None


def _build_registry() -> dict[str, DeepSeekEntry]:
    scripts_mdc: list[DeepSeekScript] = []
    for spec in [
        ("mdc-v23", Path("mdc_v23.py"), "MDC v23 bisección asimétrica 25/75", 5, "antipc deepseek run mdc-v23 --n N"),
        ("mdc-v17", Path("mdc_v17.py"), "MDC v17 pinça doble 4+4 + jerk", 5, "antipc mdc jerk --n N"),
        ("ksweep-py", Path("ksweep_predictiu.py"), "K-sweep predictivo entero (Python puro)", 5, "antipc mdc ksweep N"),
        ("benchmark-frontera", Path("benchmark_frontera.py"), "Benchmark frontera MDC vs clásico", 5, ""),
        ("record-mundial", Path("record_mundial.py"), "Récord factorización demo", 5, ""),
        ("mdc-parabola", Path("mdc_parabola_hybrid.py"), "Híbrido parabólico + K-sweep", 5, ""),
    ]:
        s = _script(spec[0], spec[1], spec[2], libro=spec[3], cmd=spec[4])
        if s:
            scripts_mdc.append(s)

    for spec in [
        ("ds-115c63", Path("deepseek_python_20260328_115c63.py"), "DeepSeek — medir_tramo V/A/jerk", 5),
        ("ds-c8a6e9", Path("deepseek_python_20260328_c8a6e9.py"), "DeepSeek — pinza doble 16 puntos", 5),
        ("ds-ff181e", Path("deepseek_python_20260328_ff181e.py"), "DeepSeek — derivadas snap", 5),
        ("mdc-libro5", Path("mdc_libro5.py"), "MDC Libro 5 consolidado", 5),
    ]:
        s = _script(spec[0], spec[1], spec[2], libro=spec[3])
        if s:
            scripts_mdc.append(s)

    docs = _collect_txt_docs()
    tabla = _existing(_DEEPSEEK_DIR / "Nuevo Documento de texto (15).txt")
    cosas = _existing(_DESKTOP_REPO / "cosas" / "deepseek 6-5-26.txt")

    return {
        "mdc-u": DeepSeekEntry(
            key="mdc-u",
            titulo="MDC-U — Método Diofántico Cinemático Unificado",
            tema="Jerk 4+4, salto balístico m−V/A, topología primorial, K-sweep",
            fuente=str(_L5_SCRIPTS),
            ecuaciones=[
                "d(m) = N mod (2(2m+3)) / (2(2m+3)) − ½",
                "J = d₃ − 3d₂ + 3d₁ − d₀  (invariante parabólico)",
                "Δm ≈ −V/A  (salto balístico si |J|<ε)",
                "r(m) = N % (2D), D=2m+3 — K-sweep predictivo",
            ],
            scripts=scripts_mdc,
            docs=[p for p in [tabla, cosas] if p],
        ),
        "libros-tabla": DeepSeekEntry(
            key="libros-tabla",
            titulo="DeepSeek jun 2026 — tabla hitos Libros 1–4",
            tema="Líneas matemáticas extraídas, criba desmemoriada, MRAUV, ZypyZape",
            fuente=str(_DEEPSEEK_DIR),
            ecuaciones=[
                "Criba desmemoriada: AND booleano entre listas replicadas",
                "D_pred ≈ D₀ + V₀·Δ + ½a₀·Δ²  (MRAUV densidad)",
                "Δm = −v/a  (salto balístico Libro 2)",
                "D_{n+1} = (1−s)D_n + s·T  (Criva iterativa)",
            ],
            scripts=[],
            docs=docs[:8] + ([tabla] if tabla and tabla not in docs[:8] else []),
        ),
        "criba-aleatorovix": DeepSeekEntry(
            key="criba-aleatorovix",
            titulo="Aleatorovix + criba desmemoriada (DeepSeek)",
            tema="Organismo digital, renacimiento memset, sincronía 121",
            fuente=str(_DEEPSEEK_DIR),
            ecuaciones=[
                "Criba desmemoriada → antipc criba --desmemoriada",
                "Aleatorovix → antipc geo-masivo",
                "K3 stream → antipc k3 stream-xor",
            ],
            scripts=[],
            docs=[p for p in docs if "15" in p.name or "14" in p.name][:3],
        ),
    }


DEEPSEEK_REGISTRY: dict[str, DeepSeekEntry] = _build_registry()


def list_deepseek() -> list[DeepSeekEntry]:
    return list(DEEPSEEK_REGISTRY.values())


def resolve_script(key: str, script_id: str | None = None) -> tuple[DeepSeekEntry, DeepSeekScript]:
    entry = DEEPSEEK_REGISTRY.get(key.lower())
    if entry is None:
        raise KeyError(f"entrada deepseek desconocida: {key}")
    if not entry.scripts:
        raise FileNotFoundError(f"{key} no tiene scripts ejecutables")
    if script_id is None:
        return entry, entry.scripts[0]
    sid = script_id.lower()
    for s in entry.scripts:
        if s.id.lower() == sid:
            return entry, s
    raise KeyError(f"script '{script_id}' no encontrado en {key}")


def format_deepseek_info(entry: DeepSeekEntry) -> str:
    lines = [
        f"[{entry.key}] {entry.titulo}",
        f"  Tema   : {entry.tema}",
        f"  Fuente : {entry.fuente}",
        "",
        "  Ecuaciones:",
    ]
    for eq in entry.ecuaciones:
        lines.append(f"    · {eq}")
    lines.append("")
    lines.append("  Scripts:")
    if entry.scripts:
        for s in entry.scripts:
            ok = "OK" if s.path.is_file() else "NO"
            lines.append(f"    · {s.id:18} [{ok}] {s.path.name}")
            lines.append(f"      {s.description}")
            if s.antipc_cmd:
                lines.append(f"      → {s.antipc_cmd}")
    else:
        lines.append("    (ninguno)")
    if entry.docs:
        lines.append("")
        lines.append("  Documentos:")
        for d in entry.docs[:12]:
            lines.append(f"    · {d}")
        if len(entry.docs) > 12:
            lines.append(f"    … +{len(entry.docs) - 12} más en {_DEEPSEEK_DIR}")
    return "\n".join(lines)


def format_deepseek_table() -> str:
    lines = [
        "=" * 88,
        "  DEEPSEEK 6 2026 → AntiPC",
        "=" * 88,
        f"  {'Clave':<16}  {'Scripts':>7}  {'Docs':>5}  Título",
        "-" * 88,
    ]
    for e in list_deepseek():
        lines.append(
            f"  {e.key:<16}  {len(e.scripts):>7}  {len(e.docs):>5}  {e.titulo}"
        )
    lines.append("=" * 88)
    return "\n".join(lines)


@dataclass
class DeepSeekRunResult:
    key: str
    script_id: str
    script: Path
    returncode: int
    elapsed_s: float
    stdout: str
    stderr: str


def run_script(
    key: str,
    script_id: str | None = None,
    *,
    extra_args: list[str] | None = None,
    timeout_s: float = 180.0,
) -> DeepSeekRunResult:
    entry, script = resolve_script(key, script_id)
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

    return DeepSeekRunResult(
        key=entry.key,
        script_id=script.id,
        script=script.path,
        returncode=proc.returncode,
        elapsed_s=elapsed,
        stdout=proc.stdout or "",
        stderr=proc.stderr or "",
    )