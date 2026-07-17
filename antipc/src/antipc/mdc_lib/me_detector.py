"""
MEcuation Factor Detector — oráculo local por SVD (Libro 6 / Newton Rápido).

Fuente filtrada:
  repo/py/me_detector.py
  filesclaude 6-5/me_detector.py

Dado una familia de enteros E, construye features logarítmicas, centra y
aplica SVD. Si κ = σ_min/σ_max colapsa y el bootstrap es robusto, hay
MEcuation local (oráculo de guess para Newton).

Estado: [HEUR] — detector empírico; no teorema.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

FEATURE_NAMES = [
    "log(E)",
    "log(E)/2",
    "log(E)/3",
    "log(E)²",
    "√log(E)",
    "log₂(E)",
    "log₃(E)",
    "log₅(E)",
    "log₆(E)",
    "log₁₀(E)",
    "log(log(E))",
    "E mod 6",
    "E mod 30",
]


def es_primo(n: int) -> bool:
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def primos_hasta(n: int) -> list[int]:
    return [p for p in range(2, n + 1) if es_primo(p)]


def feature_map(e: float) -> np.ndarray | None:
    if e <= 0:
        return None
    le = math.log(e)
    return np.array(
        [
            le,
            le / 2,
            le / 3,
            le**2,
            math.sqrt(le),
            le / math.log(2),
            le / math.log(3),
            le / math.log(5),
            le / math.log(6),
            le / math.log(10),
            math.log(le) if le > 1 else 0.0,
            e % 6,
            e % 30,
        ],
        dtype=float,
    )


def generar_familia(
    regla: str,
    n_elementos: int = 15,
    *,
    k: int = 6,
    lista: list[int] | None = None,
    primos_lim: int = 500,
) -> list[tuple[int, dict[str, Any]]]:
    primos = primos_hasta(primos_lim)
    datos: list[tuple[int, dict[str, Any]]] = []

    if regla == "cuadrados":
        for p in primos[:n_elementos]:
            datos.append((p**2, {"tipo": "p^2", "p": p}))
    elif regla == "cubos":
        for p in primos[:n_elementos]:
            datos.append((p**3, {"tipo": "p^3", "p": p}))
    elif regla == "semiprimos":
        count = 0
        for i, p in enumerate(primos):
            for q in primos[i + 1 :]:
                if count >= n_elementos:
                    break
                datos.append((p * q, {"tipo": "p*q", "p": p, "q": q}))
                count += 1
            if count >= n_elementos:
                break
    elif regla == "semiprimos_cercanos":
        for i in range(min(n_elementos, len(primos) - 1)):
            p, q = primos[i], primos[i + 1]
            datos.append((p * q, {"tipo": "p*q_cercanos", "p": p, "q": q}))
    elif regla == "k_primo":
        for p in primos[:n_elementos]:
            datos.append((k * p, {"tipo": f"{k}*p", "k": k, "p": p}))
    elif regla == "mersenne":
        count = 0
        for n in range(2, 60):
            if count >= n_elementos:
                break
            datos.append((2**n - 1, {"tipo": "2^n-1", "n": n}))
            count += 1
    elif regla == "potencia_mixta":
        count = 0
        for p in primos:
            for q in primos:
                if q == p:
                    continue
                if count >= n_elementos:
                    break
                datos.append((p**2 * q, {"tipo": "p^2*q", "p": p, "q": q}))
                count += 1
            if count >= n_elementos:
                break
    elif regla == "sophie_germain":
        count = 0
        for p in primos:
            if count >= n_elementos:
                break
            if es_primo(2 * p + 1):
                datos.append((p, {"tipo": "sophie", "p": p, "2p+1": 2 * p + 1}))
                count += 1
    elif regla == "custom":
        for e in lista or []:
            datos.append((int(e), {"tipo": "custom"}))
    else:
        raise ValueError(f"familia desconocida: {regla}")
    return datos


FAMILIAS_SCAN = [
    ("cuadrados", {}),
    ("cubos", {}),
    ("semiprimos_cercanos", {}),
    ("semiprimos", {}),
    ("k_primo", {"k": 2}),
    ("k_primo", {"k": 6}),
    ("mersenne", {}),
    ("sophie_germain", {}),
    ("potencia_mixta", {}),
]


@dataclass
class MeDetectResult:
    existe: bool
    kappa: float
    var_r1: float
    deg_frac: float
    n_datos: int
    v1: list[float]
    sigma_head: list[float]
    interpretaciones: list[tuple[str, str]] = field(default_factory=list)
    top_coef: list[tuple[str, float]] = field(default_factory=list)


def detectar_mecuation(
    datos: list[tuple[int, dict[str, Any]]],
    *,
    tau: float = 1e-4,
    bootstrap_n: int = 300,
    consensus: float = 0.90,
    seed: int = 42,
) -> MeDetectResult:
    if len(datos) < 3:
        raise ValueError("se necesitan al menos 3 observaciones")
    valores = [d[0] for d in datos]
    y = np.array([feature_map(e) for e in valores], dtype=float)
    y_c = y - y.mean(axis=0)
    _u, sigma, vt = np.linalg.svd(y_c, full_matrices=False)
    kappa = float(sigma[-1] / sigma[0]) if sigma[0] > 0 else 1.0
    var_r1 = float(sigma[0] ** 2 / (sigma**2).sum())
    v1 = vt[0]

    rng = np.random.default_rng(seed)
    deg_count = 0
    epsilons = [0.001, 0.005, 0.01]
    for _ in range(bootstrap_n):
        eps = float(rng.choice(epsilons))
        ruido = 1 + rng.uniform(-eps, eps, size=len(valores))
        vals_p = [v * r for v, r in zip(valores, ruido)]
        yp = np.array([feature_map(e) for e in vals_p], dtype=float)
        yp_c = yp - yp.mean(axis=0)
        _up, sp, _vp = np.linalg.svd(yp_c, full_matrices=False)
        kp = float(sp[-1] / sp[0]) if sp[0] > 0 else 1.0
        if kp < tau:
            deg_count += 1
    deg_frac = deg_count / bootstrap_n
    existe = kappa < tau and deg_frac >= consensus

    idx_sorted = np.argsort(np.abs(v1))[::-1]
    top_coef = [
        (FEATURE_NAMES[int(i)], float(v1[int(i)]))
        for i in idx_sorted[:5]
        if abs(v1[int(i)]) > 0.05
    ]
    interp = interpretar_vector(v1)

    return MeDetectResult(
        existe=existe,
        kappa=kappa,
        var_r1=var_r1,
        deg_frac=deg_frac,
        n_datos=len(datos),
        v1=[float(x) for x in v1],
        sigma_head=[float(x) for x in sigma[:5]],
        interpretaciones=interp,
        top_coef=top_coef,
    )


def interpretar_vector(v1: np.ndarray | list[float]) -> list[tuple[str, str]]:
    c = np.asarray(v1, dtype=float)
    out: list[tuple[str, str]] = []
    if abs(c[1]) > 1e-6:
        ratio = float(c[0] / c[1])
        if abs(ratio - 2.0) < 0.15:
            out.append(("cuadrado", "E = p² → factor = √E = exp(log(E)/2)"))
        if abs(ratio - 3.0) < 0.2:
            out.append(("cubo", "E = p³ → factor = ∛E = exp(log(E)/3)"))
    if abs(c[1]) > 0.5 and abs(c[0]) > 0.5:
        out.append(
            (
                "potencia_par",
                f"posible potencia par, ratio c0/c1={c[0] / c[1]:.3f}",
            )
        )
    return out


def oracle_j_inicial(e: float, familia: str = "cuadrados") -> float:
    """Guess j para Newton Rápido según familia (oráculo MEcuation)."""
    if e <= 0:
        raise ValueError("E debe ser > 0")
    le = math.log(e)
    if familia == "cuadrados":
        return le / 2
    if familia == "cubos":
        return le / 3
    if familia == "general":
        return 1.0
    return le / 2


def scan_familias(
    n: int = 12,
    *,
    bootstrap_n: int = 200,
    seed: int = 42,
) -> list[tuple[str, MeDetectResult]]:
    resumen: list[tuple[str, MeDetectResult]] = []
    for regla, kwargs in FAMILIAS_SCAN:
        nombre = regla + (f"(k={kwargs['k']})" if "k" in kwargs else "")
        datos = generar_familia(regla, n_elementos=n, **kwargs)
        if len(datos) < 3:
            continue
        res = detectar_mecuation(
            datos, bootstrap_n=bootstrap_n, seed=seed
        )
        resumen.append((nombre, res))
    return resumen


def format_detect(nombre: str, r: MeDetectResult) -> str:
    estado = "ME local" if r.existe else "sin ME global"
    lines = [
        f"Detector MEcuation — {nombre} (n={r.n_datos})",
        f"  κ (colapso)  : {r.kappa:.2e}  {'COLAPSO' if r.kappa < 1e-4 else 'no colapsa'}",
        f"  Var r=1      : {r.var_r1 * 100:.1f}%",
        f"  Bootstrap    : {r.deg_frac:.2f}",
        f"  Estado       : {estado}",
        f"  σ head       : {[round(x, 4) for x in r.sigma_head]}",
    ]
    if r.top_coef:
        lines.append("  Coeficientes :")
        for name, c in r.top_coef:
            lines.append(f"    {name:16s}  c={c:+.6f}")
    for _tag, desc in r.interpretaciones:
        lines.append(f"  → {desc}")
    return "\n".join(lines)


def format_scan(rows: list[tuple[str, MeDetectResult]]) -> str:
    lines = [
        "Scan de familias — MEcuation (SVD + bootstrap)",
        f"  {'familia':30s}  {'κ':>10}  {'boot':>6}  estado",
        "  " + "-" * 62,
    ]
    for nombre, r in rows:
        flag = "ME" if r.existe else "--"
        lines.append(
            f"  {nombre:30s}  {r.kappa:10.1e}  {r.deg_frac:6.2f}  {flag}"
        )
        for _tag, desc in r.interpretaciones:
            lines.append(f"    ↳ {desc}")
    return "\n".join(lines)


def format_oracle(e: float, familia: str) -> str:
    j = oracle_j_inicial(e, familia)
    lines = [
        f"Oráculo Newton — familia={familia}  E={e}",
        f"  j_inicial  : {j:.12f}",
    ]
    if familia == "cuadrados":
        p = math.exp(j)
        lines.append(f"  factor √E  : {p:.12f}")
        lines.append(f"  check p²   : {p * p:.6f}")
    elif familia == "cubos":
        p = math.exp(j)
        lines.append(f"  factor ∛E  : {p:.12f}")
        lines.append(f"  check p³   : {p ** 3:.6f}")
    return "\n".join(lines)
