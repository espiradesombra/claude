"""
Fichas teoremas VMA (repo/teoremas) → CLI AntiPC.

Lee INDICE_MAESTRO.txt y fichas 01–30.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

_PKG = Path(__file__).resolve().parents[1]
_REPO_ANTIPC = _PKG.parents[1]
_TEOREMAS = _REPO_ANTIPC.parent / "teoremas"


@dataclass
class TeoremaInfo:
    numero: int
    id: str
    estado: str
    archivo: str

    @property
    def path(self) -> Path:
        return _TEOREMAS / self.archivo


_INDICE_LINE = re.compile(
    r"^\s*(\d{1,2})\s*\|\s*(\S+)\s*\|\s*\[(\w+)\]\s*\|\s*(.+)$"
)


def _parse_indice() -> list[TeoremaInfo]:
    indice = _TEOREMAS / "INDICE_MAESTRO.txt"
    if not indice.is_file():
        return _fallback_registry()

    out: list[TeoremaInfo] = []
    for line in indice.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _INDICE_LINE.match(line)
        if not m:
            continue
        num = int(m.group(1))
        if num > 30:
            continue
        archivo = m.group(4).strip().split("|")[0].strip()
        out.append(
            TeoremaInfo(
                numero=num,
                id=m.group(2),
                estado=m.group(3),
                archivo=archivo,
            )
        )
    return out if out else _fallback_registry()


def _fallback_registry() -> list[TeoremaInfo]:
    items: list[TeoremaInfo] = []
    for p in sorted(_TEOREMAS.glob("[0-9][0-9]_*.txt")):
        num = int(p.name[:2])
        tid = p.stem.split("_", 1)[-1][:24]
        items.append(TeoremaInfo(numero=num, id=tid, estado="?", archivo=p.name))
    return items


_REGISTRY: list[TeoremaInfo] | None = None


def list_teoremas() -> list[TeoremaInfo]:
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = _parse_indice()
    return list(_REGISTRY)


def resolve_teorema(numero: int | None = None, teorema_id: str | None = None) -> TeoremaInfo:
    items = list_teoremas()
    if numero is not None:
        for t in items:
            if t.numero == numero:
                return t
        raise KeyError(f"teorema {numero} no encontrado")
    if teorema_id:
        tid = teorema_id.lower().replace("-", "_")
        for t in items:
            if t.id.lower() == tid or tid in t.archivo.lower():
                return t
        raise KeyError(f"teorema id '{teorema_id}' no encontrado")
    raise ValueError("indica número o id")


def read_teorema_body(info: TeoremaInfo, *, max_lines: int = 80) -> str:
    p = info.path
    if not p.is_file():
        return f"(ficha no encontrada: {p})"
    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) <= max_lines:
        return "\n".join(lines)
    head = "\n".join(lines[:max_lines])
    return f"{head}\n… ({len(lines) - max_lines} líneas más en {p.name})"


def format_teorema_info(info: TeoremaInfo, *, body: bool = True) -> str:
    lines = [
        f"Teorema {info.numero:02d} — {info.id}",
        f"  Estado  : [{info.estado}]",
        f"  Archivo : {info.path}",
    ]
    if body:
        lines.append("")
        lines.append(read_teorema_body(info))
    return "\n".join(lines)


def format_teoremas_table() -> str:
    items = list_teoremas()
    lines = [
        "=" * 72,
        "  TEOREMAS VMA — índice maestro",
        "=" * 72,
        f"  {'Nº':>3}  {'ID':<22}  {'Estado':<8}  Archivo",
        "-" * 72,
    ]
    for t in items:
        ok = "OK" if t.path.is_file() else "NO"
        lines.append(
            f"  {t.numero:>3}  {t.id:<22}  [{t.estado:<6}]  {t.archivo}  [{ok}]"
        )
    lines.append("=" * 72)
    return "\n".join(lines)