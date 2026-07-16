"""
Puente AntiPC → goldbach/goldbach_vma.py (paper sympy, Libro 1/6).

No duplica el artículo completo; carga el módulo canónico del repo.
"""

from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

_REPO = Path(__file__).resolve().parents[4]
_GB_PATH = _REPO / "goldbach" / "goldbach_vma.py"
_mod = None


def _load() -> object:
    global _mod
    if _mod is not None:
        return _mod
    if not _GB_PATH.is_file():
        raise FileNotFoundError(f"goldbach_vma.py no encontrado: {_GB_PATH}")
    spec = importlib.util.spec_from_file_location("goldbach_vma", _GB_PATH)
    if spec is None or spec.loader is None:
        raise ImportError("no se pudo cargar goldbach_vma")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["goldbach_vma"] = mod
    spec.loader.exec_module(mod)
    _mod = mod
    return mod


@dataclass
class GoldbachAnalyze:
    n2: int
    r: int
    pp: int
    pc: int
    cp: int
    cc: int
    goldbach_cierta: bool
    brecha: int
    descomps: list[tuple[int, int]]


def analyze(n2: int) -> GoldbachAnalyze:
    if n2 < 4 or n2 % 2 != 0:
        raise ValueError("2n debe ser par y >= 4")
    gb = _load()
    card = gb.conjuntos_cardinalidad(n2)
    ie = gb.propiedades_ie_exactas(n2)
    return GoldbachAnalyze(
        n2=n2,
        r=card["R"],
        pp=ie["PP"],
        pc=ie["PC"],
        cp=ie["CP"],
        cc=ie["CC"],
        goldbach_cierta=ie["goldbach_cierta"],
        brecha=card["brecha"],
        descomps=card["descomps"][:8],
    )


def verify_range(n2_max: int = 200, *, step: int = 2) -> tuple[bool, list[int]]:
    gb = _load()
    fallos: list[int] = []
    for n2 in range(6, n2_max + 1, step):
        ie = gb.propiedades_ie_exactas(n2)
        if not ie["goldbach_cierta"]:
            fallos.append(n2)
    return len(fallos) == 0, fallos


def format_analyze(a: GoldbachAnalyze) -> str:
    des = ", ".join(f"({p},{q})" for p, q in a.descomps[:6])
    if len(a.descomps) > 6:
        des += "..."
    return "\n".join(
        [
            f"Goldbach VMA — 2n={a.n2}",
            f"  R descomps: {a.r}",
            f"  PP/PC/CP/CC: {a.pp}/{a.pc}/{a.cp}/{a.cc}",
            f"  Goldbach    : {'SÍ' if a.goldbach_cierta else 'NO'}",
            f"  Brecha §4   : {a.brecha}",
            f"  Muestra     : {des or '—'}",
        ]
    )


def format_verify(n2_max: int, ok_all: bool, fallos: list[int]) -> str:
    if ok_all:
        return f"Goldbach VMA verify — 2n≤{n2_max}: OK (PP>0 en todos los pares)"
    muestra = fallos[:10]
    return "\n".join(
        [
            f"Goldbach VMA verify — 2n≤{n2_max}: FALLO",
            f"  Contraejemplos ({len(fallos)}): {muestra}",
        ]
    )