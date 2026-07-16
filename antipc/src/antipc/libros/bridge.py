"""
Métodos matemáticos Libros 1–6 (corpus VMA) → comandos AntiPC.

Fuentes canónicas:
  Libro 1  txt l5/recopilacion de libros.txt  (Números i numeritos)
  Libro 2  Libro2 Números otra vez/
  Libro 3  Sigo en mis trece (+ lee arbusto)
  Libro 4  lee arbusto/libro4_encriptacion_convergencias.py, encriptacionGeometrica/
  Libro 5  Libro5 Factorizacion con 2v+3/, PY L5/mdc_libro5.py
  Libro 6  Libro6 NewtonRapido..., teoremas/21–26
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
_TEOREMAS = _DESKTOP_REPO / "teoremas"


@dataclass
class LibroMetodo:
    id: str
    nombre: str
    ecuacion: str
    antipc: str
    estado: str  # integrado | parcial | pendiente
    notas: str = ""


@dataclass
class LibroInfo:
    numero: int
    titulo: str
    tema: str
    fuente: str
    teoremas: list[str] = field(default_factory=list)
    metodos: list[LibroMetodo] = field(default_factory=list)


def _build_registry() -> dict[int, LibroInfo]:
    return {
        1: LibroInfo(
            numero=1,
            titulo="Números i numeritos",
            tema="Bases, primos, Andrica, Goldbach, criba desmemoriada",
            fuente=str(_DESKTOP_REPO / "txt l5" / "recopilacion de libros.txt"),
            teoremas=["01 dos primos", "03 salto máximo", "06 Goldbach VMA", "16 Andrica"],
            metodos=[
                LibroMetodo(
                    "L1-desmemoriada",
                    "Criba desmemoriada",
                    "AND booleano entre listas replicadas (posición = n)",
                    "antipc criba --desmemoriada",
                    "integrado",
                    "vma-methods CribaDesmemoriada (Python bridge v0.14)",
                ),
                LibroMetodo(
                    "L1-6k-rejilla",
                    "Rejilla 6k±1",
                    "{2·3·n±1} − {30m±5, 6·7·ℓ±7, …}",
                    "antipc criba --modular6k",
                    "integrado",
                    "DLL 0.10.0-c criba_modular6k",
                ),
                LibroMetodo(
                    "L1-salto",
                    "Salto máximo",
                    "SaltoMaximo(n) ≈ √n (peor caso 6k±1)",
                    "antipc salto ventana | verify",
                    "integrado",
                    "mdc_lib/salto_maximo.py — filtrado filesclaude 6-5 2026-07-16",
                ),
                LibroMetodo(
                    "L1-goldbach-asim",
                    "Goldbach por asimetría",
                    "f(x)·2 ≤ f(2x) − acumularFALLO",
                    "antipc goldbach analyze | verify",
                    "integrado",
                    "goldbach/goldbach_vma.py + bridge 2026-07-16",
                ),
                LibroMetodo(
                    "L1-eratostenes",
                    "Eratóstenes optimizado",
                    "Segmentación + límites",
                    "antipc criba --limit N",
                    "integrado",
                ),
            ],
        ),
        2: LibroInfo(
            numero=2,
            titulo="Números oTra VeZ",
            tema="Dos primos, MRAUV, criba híbrida, desfase espejo",
            fuente=str(_DESKTOP_REPO / "Libro2 Números otra vez"),
            teoremas=["01 dos primos", "02 intervalo", "14 MRAUV Goldbach"],
            metodos=[
                LibroMetodo(
                    "L2-dos-primos",
                    "Criterio dos primos",
                    "I(n): L(n)−m(n) ≥ 2 ⇒ ≥2 primos",
                    "antipc dos-primos verify | tabla | audit",
                    "integrado",
                    "mdc_lib/dos_primos.py + mrauv/datos/ anexoE 2026-07-16",
                ),
                LibroMetodo(
                    "L2-mrauv",
                    "MRAUV densidad",
                    "D ≈ D₀ + V₀·Δ + ½a₀·Δ² + (⅙)j₀·Δ³",
                    "antipc mrauv calibrar | densidad | goldbach",
                    "integrado",
                    "mdc_lib/mrauv.py — filtrado archivos/ 2026-07-16",
                ),
                LibroMetodo(
                    "L2-criba-hibrida",
                    "Criba híbrida",
                    "Ascendente + descendente con semilla",
                    "antipc criba --hibrida",
                    "integrado",
                    "DLL 0.4.0-c",
                ),
                LibroMetodo(
                    "L2-espejo",
                    "Desfase espejo",
                    "(n % x) simetría con restos espejo",
                    "antipc mdc analyze",
                    "parcial",
                    "Fase dual A/B en MDC analyze",
                ),
            ],
        ),
        3: LibroInfo(
            numero=3,
            titulo="Sigo en mis Trece",
            tema="Sofí, Criva, densidad fractal, primos decrecen K=9/24",
            fuente="lee arbusto + teoremas 04–06, 28–30",
            teoremas=["04 Sofí L4−L2", "05 U2 infinito", "30 primos decrecen K"],
            metodos=[
                LibroMetodo(
                    "L3-sofi",
                    "Estructuras Sofí",
                    "L1={6k+5}, U2=L1\\(L3∪L4) ⊆ LSG",
                    "antipc sofi classify",
                    "integrado",
                    "mdc_lib/sofi.py — filtrado filesclaude 6-5 2026-07-16",
                ),
                LibroMetodo(
                    "L3-criva",
                    "Criva iterativa",
                    "D_{n+1}=(1−s)D_n + s·T, T=π(x)/x",
                    "antipc criba + vma-methods criva",
                    "integrado",
                    "criva.c en DLL + bridge Python",
                ),
                LibroMetodo(
                    "L3-densidad-fractal",
                    "Densidad por capas",
                    "π(x) ≈ x·Σ D₀/2ⁿ·sigmoid(x−10ⁿ)",
                    "vma-methods criva --compare",
                    "parcial",
                ),
                LibroMetodo(
                    "L3-K-9-24",
                    "Primos decrecen K",
                    "K=⌊⌊√n⌋·9/24⌋ en ventana local",
                    "teoremas/30_primos_decrecen_K_9_24.txt",
                    "parcial",
                    "Enlaza ficha 01 dos primos",
                ),
                LibroMetodo(
                    "L3-pitagorico",
                    "Método pitagórico visual",
                    "(2x+3)(2y+3)=n por rectas en plano",
                    "antipc mdc visual + analyze",
                    "integrado",
                    "Ficha 29 teoremas",
                ),
            ],
        ),
        4: LibroInfo(
            numero=4,
            titulo="Encriptación y Energía Verde",
            tema="Convergencia geométrica, K3, ZypyZape, Quijote, Kilòmetre",
            fuente=str(_DESKTOP_REPO / "lee arbusto" / "libro4_encriptacion_convergencias.py"),
            teoremas=["13 zona densa e", "17 convergencia", "19 PhaseAmplifier K3"],
            metodos=[
                LibroMetodo(
                    "L4-convergencia",
                    "Convergencia binaria",
                    "Pares 10/11/00/01 → perímetro Thales",
                    "antipc geo",
                    "integrado",
                    "DLL geo_converge; Decimal prec≥100",
                ),
                LibroMetodo(
                    "L4-aleatorovix",
                    "Aleatorovix π/e",
                    "Perímetro × π ofuscado × e convergente",
                    "antipc geo-masivo",
                    "integrado",
                    "DLL k3_geo_aleatorovix",
                ),
                LibroMetodo(
                    "L4-k3-stream",
                    "Motor acordeón K3",
                    "L+=v+1, v+=2; XOR (L^v)×golden",
                    "antipc k3 stream-xor",
                    "integrado",
                    "base=33, rel=1",
                ),
                LibroMetodo(
                    "L4-zypyzape",
                    "Gemelo ZypyZape",
                    "J·dω/dt = T_aero − T_gen + T_ZZ",
                    "antipc gemelo run zypyzape",
                    "integrado",
                ),
                LibroMetodo(
                    "L4-quijote",
                    "Actuador Quijote",
                    "J_i(r)=J_G+N_b·m_q·r²; J·ω̇=T−ω·J̇",
                    "antipc gemelo run quijote --variant v5|v10",
                    "integrado",
                ),
                LibroMetodo(
                    "L4-gemelo-v10",
                    "Gemelo v10 Cp dinámico",
                    "Cp(λ,β) NREL 5MW; pitch asíncrono por pala",
                    "antipc gemelo run quijote --variant v10",
                    "integrado",
                    "2026/gemelo_v10_cp_dinamic.py — filtrado filesclaude 6-5",
                ),
                LibroMetodo(
                    "L4-kilometre",
                    "Kilòmetre flotación",
                    "F_n=(ρ_f−ρ_o)·V·g; E_k=½Iω²",
                    "antipc gemelo run kilometre",
                    "integrado",
                ),
            ],
        ),
        5: LibroInfo(
            numero=5,
            titulo="Factorización con (2v+3)(2y+3)",
            tema="MDC, pinza, numeritos, jerk, K-sweep, dientes de sierra",
            fuente=str(_DESKTOP_REPO / "Libro5 Factorizacion con 2v+3"),
            teoremas=["09 MDC teorema1", "10 finestra local", "11 atrappament", "12 cambio signo", "20 cofactor"],
            metodos=[
                LibroMetodo(
                    "L5-regla-mecanica",
                    "Regla mecánica",
                    "(2v+3)(2b+3) = N",
                    "antipc mdc factor",
                    "integrado",
                    "Toy ≤10 dígitos",
                ),
                LibroMetodo(
                    "L5-dos-trenes",
                    "Dos trenes MDC",
                    "Colisiones (x,y) + scan entero C",
                    "antipc mdc analyze",
                    "integrado",
                    "mdc_trains.c DLL 0.4.0-c",
                ),
                LibroMetodo(
                    "L5-ksweep",
                    "K-sweep predictivo",
                    "Barrido D con predicción entera",
                    "antipc mdc ksweep",
                    "integrado",
                    "mdc_ksweep.c DLL 0.5.0-c",
                ),
                LibroMetodo(
                    "L5-jerk",
                    "Invariante Jerk",
                    "|J|<ε ⇒ salto balístico m−v/a",
                    "antipc mdc jerk --n N",
                    "integrado",
                    "mdc_jerk.py + mdc_v17/v23 DeepSeek L5",
                ),
                LibroMetodo(
                    "L5-numeritos",
                    "Numeritos / restos",
                    "Onda de restos + dientes de sierra",
                    "antipc mdc analyze --n N",
                    "parcial",
                ),
                LibroMetodo(
                    "L5-pinza",
                    "Pinza / atrappament",
                    "Δ(S)=S²+6S−(N−9)=k²",
                    "antipc discriminant factor | trajectory",
                    "integrado",
                    "mdc_lib/discriminant.py — filtrado filesclaude 6-5 2026-07-16",
                ),
            ],
        ),
        6: LibroInfo(
            numero=6,
            titulo="Pseudocódigo de Implicaciones",
            tema="Sofí 15 impl, Método-V, MRAUV, Goldbach hueco, Newton, plantillas 1–7",
            fuente=str(_DESKTOP_REPO / "Libro6 NewtonRapido en busqueda con oraculo"),
            teoremas=["21–26 libro6", "25 método-V", "22 Sofí 15 líneas"],
            metodos=[
                LibroMetodo(
                    "L6-newton",
                    "Newton rápido + oráculo",
                    "MEcuation: familias cuadrados,cubos,potencia…",
                    "antipc newton E -f cuadrados",
                    "integrado",
                    "newton_rapido.c DLL 0.5.0-c",
                ),
                LibroMetodo(
                    "L6-metodo-V",
                    "Método-V factorial",
                    "n>2v−2>m>v−1; V=v·n!/m!",
                    "teoremas/25_metodo_V_implicaciones.txt",
                    "parcial",
                    "Heurística gemelos/Mersenne/Haus",
                ),
                LibroMetodo(
                    "L6-sofi-impl",
                    "Sofí 15 implicaciones",
                    "SI n∈L1 ∧ n∉L3 ∧ n∉L4 ⇒ n∈U2",
                    "teoremas/22 + formal/F22",
                    "parcial",
                ),
                LibroMetodo(
                    "L6-goldbach-hueco",
                    "Goldbach hueco L−m",
                    "Oscilaciones locales dominadas por hueco",
                    "teoremas/24_goldbach_hueco_implicaciones.txt",
                    "parcial",
                ),
                LibroMetodo(
                    "L6-plantillas",
                    "Plantillas 1–7",
                    "Esquemas SI A ⇒ B universales",
                    "teoremas/26_plantillas_1_a_7.txt",
                    "parcial",
                ),
                LibroMetodo(
                    "L6-fermat",
                    "Fermat alineación",
                    "∏F_k = F_{n+1}−2; CRT residuo privilegiado",
                    "antipc fermat align",
                    "integrado",
                    "mdc_lib/fermat_modular.py — filtrado diamante/ 2026-07-16",
                ),
            ],
        ),
    }


LIBRO_REGISTRY: dict[int, LibroInfo] = _build_registry()


def list_libros() -> list[LibroInfo]:
    return [LIBRO_REGISTRY[i] for i in sorted(LIBRO_REGISTRY)]


def list_metodos(*, libro: int | None = None) -> list[tuple[LibroInfo, LibroMetodo]]:
    out: list[tuple[LibroInfo, LibroMetodo]] = []
    for info in list_libros():
        if libro is not None and info.numero != libro:
            continue
        for m in info.metodos:
            out.append((info, m))
    return out


def resolve_metodo(libro: int, metodo_id: str | None = None) -> tuple[LibroInfo, LibroMetodo]:
    info = LIBRO_REGISTRY.get(libro)
    if info is None:
        raise KeyError(f"libro {libro} no existe (rango 1–6)")
    if not info.metodos:
        raise FileNotFoundError(f"libro {libro} sin métodos registrados")
    if metodo_id is None:
        return info, info.metodos[0]
    mid = metodo_id.lower()
    for m in info.metodos:
        if m.id.lower() == mid or m.nombre.lower() == mid:
            return info, m
    raise KeyError(f"método '{metodo_id}' no encontrado en libro {libro}")


def format_libro_info(info: LibroInfo) -> str:
    lines = [
        f"Libro {info.numero}: {info.titulo}",
        f"  Tema    : {info.tema}",
        f"  Fuente  : {info.fuente}",
        f"  Teoremas: {', '.join(info.teoremas) if info.teoremas else '—'}",
        "",
        "  Métodos:",
    ]
    for m in info.metodos:
        lines.append(f"    [{m.estado:10}] {m.id}")
        lines.append(f"      {m.nombre}")
        lines.append(f"      {m.ecuacion}")
        lines.append(f"      → {m.antipc}")
        if m.notas:
            lines.append(f"      ({m.notas})")
    return "\n".join(lines)


def format_metodos_table(*, libro: int | None = None) -> str:
    rows = list_metodos(libro=libro)
    lines = [
        "=" * 88,
        "  MÉTODOS LIBROS 1–6 VMA → AntiPC",
        "=" * 88,
        f"  {'Lib':>3}  {'ID':<18}  {'Estado':<10}  Método / comando AntiPC",
        "-" * 88,
    ]
    for info, m in rows:
        lines.append(
            f"  {info.numero:>3}  {m.id:<18}  {m.estado:<10}  {m.nombre}"
        )
        lines.append(f"       {m.ecuacion}")
        lines.append(f"       → {m.antipc}")
    lines.append("=" * 88)
    return "\n".join(lines)


# Comandos AntiPC ejecutables por id de método (solo estado integrado/parcial con CLI)
_RUN_BY_ID: dict[str, list[str]] = {
    "L1-desmemoriada": ["criba", "--limit", "10000", "--desmemoriada"],
    "L1-salto": ["salto", "ventana", "1000"],
    "L1-goldbach-asim": ["goldbach", "analyze", "100"],
    "L1-6k-rejilla": ["criba", "--limit", "50000", "--modular6k"],
    "L1-eratostenes": ["criba", "--limit", "10000"],
    "L2-criba-hibrida": ["criba", "--limit", "50000", "--hibrida"],
    "L2-mrauv": ["mrauv", "calibrar", "100000"],
    "L2-dos-primos": ["dos-primos", "verify", "10000"],
    "L2-espejo": ["mdc", "analyze", "--n", "1147"],
    "L3-criva": ["criba", "--limit", "10000"],
    "L3-sofi": ["sofi", "classify", "5000"],
    "L3-pitagorico": ["mdc", "analyze", "--n", "1147"],
    "L4-convergencia": ["geo", "--demo"],
    "L4-aleatorovix": ["geo-masivo", "--text", "vma-aleatorovix", "--semilla", "43210"],
    "L4-k3-stream": ["k3", "stream-xor", "--text", "libro4-33x1"],
    "L4-zypyzape": ["gemelo", "run", "zypyzape", "--hubs", "2", "--packets", "5000"],
    "L4-quijote": ["gemelo", "run", "quijote", "--variant", "v5"],
    "L4-gemelo-v10": ["gemelo", "run", "quijote", "--variant", "v10", "--no-plot"],
    "L4-kilometre": ["gemelo", "run", "kilometre"],
    "L5-regla-mecanica": ["mdc", "factor", "1147"],
    "L5-dos-trenes": ["mdc", "analyze", "--n", "1147"],
    "L5-ksweep": ["mdc", "ksweep", "1147"],
    "L5-jerk": ["mdc", "jerk", "1147"],
    "L5-numeritos": ["mdc", "analyze", "--n", "10473029"],
    "L5-pinza": ["discriminant", "factor", "11449"],
    "L6-newton": ["newton", "121", "-f", "cuadrados"],
    "L6-fermat": ["fermat", "align", "3"],
}


@dataclass
class MetodoRunResult:
    libro: int
    metodo_id: str
    comando: list[str]
    returncode: int
    elapsed_s: float
    stdout: str
    stderr: str


def run_metodo(
    libro: int,
    metodo_id: str | None = None,
    *,
    extra_n: int | None = None,
) -> MetodoRunResult:
    info, metodo = resolve_metodo(libro, metodo_id)
    argv = _RUN_BY_ID.get(metodo.id)
    if argv is None:
        raise RuntimeError(
            f"método {metodo.id} no tiene comando ejecutable "
            f"(estado={metodo.estado}); ver: {metodo.antipc}"
        )

    argv = list(argv)
    if extra_n is not None:
        for i, a in enumerate(argv):
            if a == "--n" and i + 1 < len(argv):
                argv[i + 1] = str(extra_n)
            if a.isdigit() and i > 0 and argv[i - 1] in ("factor", "ksweep"):
                argv[i] = str(extra_n)

    cli = _PKG / "cli.py"
    comando = [sys.executable, str(cli), *argv]
    t0 = time.perf_counter()
    proc = subprocess.run(
        comando,
        cwd=str(_PKG),
        capture_output=True,
        text=True,
        timeout=120.0,
    )
    elapsed = time.perf_counter() - t0
    return MetodoRunResult(
        libro=libro,
        metodo_id=metodo.id,
        comando=argv,
        returncode=proc.returncode,
        elapsed_s=elapsed,
        stdout=proc.stdout or "",
        stderr=proc.stderr or "",
    )