#!/usr/bin/env python3
"""
MEcuation / conexión nula — bootstrap + SVD con ruido ±ε
Extraído de: Pseudocódigo de implicaciones (.eml / chat dot .html)
Autor corpus: Víctor Manzanares Alberola
"""
from __future__ import annotations

import math
import random
from typing import List, Sequence, Tuple

import numpy as np
from numpy.linalg import svd


def build_feature_vector(E: float) -> np.ndarray:
    """Features logarítmicas: f3 ≡ 0 tautológicamente; f1,f2 capturan escala."""
    f1 = math.log(E)
    f2 = math.log(math.sqrt(E))
    f3 = f1 - 2.0 * f2
    return np.array([f1, f2, f3], dtype=float)


def center_columns(Y: np.ndarray) -> np.ndarray:
    return Y - np.mean(Y, axis=0)


def test_connection_null(
    observations: Sequence[float],
    epsilons: Sequence[Sequence[float]],
    tau: float = 1e-6,
    bootstrap_N: int = 300,
    seed: int | None = None,
) -> float:
    """
    Fracción de réplicas bootstrap donde σ_min/σ_max < tau (colapso).
    Muestrea cada feature: base[j] + U(-ε_ij, +ε_ij).
    """
    rng = random.Random(seed)
    x = len(observations)
    d = len(build_feature_vector(observations[0]))
    deg_count = 0

    for _ in range(bootstrap_N):
        Ys = np.zeros((x, d))
        for i in range(x):
            base = build_feature_vector(observations[i])
            samp = np.array(
                [
                    base[j] + rng.uniform(-epsilons[i][j], epsilons[i][j])
                    for j in range(d)
                ]
            )
            Ys[i, :] = samp
        Yc = center_columns(Ys)
        try:
            S = svd(Yc, compute_uv=False)
        except Exception:
            S = np.linalg.svd(Yc, compute_uv=False)
        if len(S) == 0:
            continue
        kappa = S[-1] / (S[0] if S[0] != 0 else 1e-30)
        if kappa < tau:
            deg_count += 1
    return deg_count / bootstrap_N


def propose_new_observations(
    observations: Sequence[float],
    epsilons: Sequence[Sequence[float]],
    k: int = 4,
) -> List[Tuple[float, List[float]]]:
    """Guía búsqueda: ortogonal a PC1, extremos log, reducción de ε."""
    Ys = np.array([build_feature_vector(E) for E in observations])
    Yc = center_columns(Ys)
    _, _, Vt = svd(Yc, full_matrices=False)
    pc1 = Vt[0] if Vt.shape[0] > 0 else np.zeros(Yc.shape[1])
    mean_vec = np.mean(Ys, axis=0)
    default_eps = [0.1, 0.05, 0.02]

    candidates: List[Tuple[float, List[float]]] = []

    ort_dir = np.random.normal(size=mean_vec.shape)
    ort_dir = ort_dir - (ort_dir.dot(pc1)) * pc1
    if np.linalg.norm(ort_dir) == 0:
        ort_dir = np.random.normal(size=mean_vec.shape)
    ort_dir = ort_dir / np.linalg.norm(ort_dir)
    new_log = mean_vec[0] + 0.25 * ort_dir[0]
    candidates.append((math.exp(new_log), default_eps.copy()))

    logs = [math.log(E) for E in observations]
    minlog, maxlog = min(logs), max(logs)
    candidates.append((math.exp(minlog - 0.3), default_eps.copy()))
    candidates.append((math.exp(maxlog + 0.3), default_eps.copy()))
    candidates.append(
        (observations[0], [e * 0.3 for e in epsilons[0]])
    )

    return candidates[:k]


def run_iterative_search(
    observations: List[float],
    epsilons: List[List[float]],
    tau: float = 1e-6,
    bootstrap_N: int = 300,
    consensus: float = 0.95,
    max_x: int = 20,
    add_k: int = 4,
    seed: int = 42,
) -> dict:
    round_idx = 0
    obs = list(observations)
    eps = [list(e) for e in epsilons]

    while True:
        round_idx += 1
        deg_frac = test_connection_null(obs, eps, tau=tau, bootstrap_N=bootstrap_N, seed=seed)
        if deg_frac >= consensus:
            return {
                "collapse": True,
                "rounds": round_idx,
                "n_obs": len(obs),
                "deg_frac": deg_frac,
                "observations": obs,
            }
        if len(obs) >= max_x:
            return {
                "collapse": False,
                "rounds": round_idx,
                "n_obs": len(obs),
                "deg_frac": deg_frac,
                "observations": obs,
            }
        for Enew, eps_new in propose_new_observations(obs, eps, k=add_k):
            obs.append(Enew)
            eps.append(eps_new)


def main() -> None:
    # Caso canónico: cuadrados perfectos 5², 7², 11²
    obs = [25.0, 49.0, 121.0]
    eps_examples = [
        [0.4, 0.2, 0.01],
        [0.44, 0.22, 0.01],
        [0.444, 0.222, 0.01],
    ]

    deg = test_connection_null(obs, eps_examples, bootstrap_N=300, seed=42)
    print(f"Cuadrados {{25,49,121}}: deg_frac={deg:.4f}")

    result = run_iterative_search(obs, eps_examples, bootstrap_N=200, seed=42)
    print(f"Iterativo: collapse={result['collapse']}, rounds={result['rounds']}, "
          f"n={result['n_obs']}, deg_frac={result['deg_frac']:.4f}")


if __name__ == "__main__":
    main()