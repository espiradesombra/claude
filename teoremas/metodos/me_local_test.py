#!/usr/bin/env python3
"""
ME locales (cuadrática/cúbica) y test de NO-unificación — bootstrap + SVD
Extraído de: Pseudocódigo de implicaciones / Nuevo Documento de texto2.txt
Autor corpus: Víctor Manzanares Alberola
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from math import exp, log
from typing import Dict, List, Optional, Tuple

import numpy as np

PRIMES_SQUARE = [5, 7, 11]
PRIMES_CUBE = [5, 7, 11]

EPS_SQ = 0.15
EPS_CB = 0.18

TAU = 5e-3
CONSENSUS = 0.90
BOOTSTRAP_N = 900
SEED = 2026


@dataclass
class LogObs:
    y: float
    eps: float
    tag: str = ""


@dataclass
class MEQTesterV3:
    tau: float = TAU
    consensus: float = CONSENSUS
    bootstrap_N: int = BOOTSTRAP_N
    seed: int = SEED
    obs: List[LogObs] = field(default_factory=list)

    def add_log_obs(self, y: float, eps: float, tag: str = "") -> None:
        self.obs.append(LogObs(y=y, eps=eps, tag=tag))

    @staticmethod
    def _near_k_power_residual(E: float, k: int) -> float:
        x = E ** (1.0 / k)
        r = round(x)
        return abs(x - r)

    def _features_from_E(self, E: float) -> np.ndarray:
        y = log(max(E, 1e-300))
        d2 = self._near_k_power_residual(E, 2)
        d3 = self._near_k_power_residual(E, 3)
        d5 = self._near_k_power_residual(E, 5)
        f1 = y
        f2 = np.log(1.0 + d2)
        f3 = np.log(1.0 + d3)
        f4 = np.log(1.0 + d5)
        mn = min(d2, d3) + 1e-12
        mx = max(d2, d3) + 1e-12
        f5 = np.log(mx / mn)
        return np.array([f1, f2, f3, f4, f5], dtype=float)

    def _sample_matrix(
        self, rng: np.random.Generator, indices: Optional[List[int]] = None
    ) -> np.ndarray:
        rows = []
        src = self.obs if indices is None else [self.obs[i] for i in indices]
        for ob in src:
            y_s = rng.uniform(ob.y - ob.eps, ob.y + ob.eps)
            E_s = exp(y_s)
            rows.append(self._features_from_E(E_s))
        Y = np.vstack(rows)
        return Y - Y.mean(axis=0, keepdims=True)

    def _degenerate_fraction(self, indices: Optional[List[int]] = None) -> float:
        rng = np.random.default_rng(self.seed)
        deg = 0
        for _ in range(self.bootstrap_N):
            Yc = self._sample_matrix(rng, indices)
            if Yc.shape[0] < 3 or Yc.shape[1] < 2:
                continue
            S = np.linalg.svd(Yc, full_matrices=False, compute_uv=False)
            if S[0] <= 0:
                continue
            if (S[-1] / S[0]) < self.tau:
                deg += 1
        return deg / max(1, self.bootstrap_N)

    def local_tests(self, window_indices: List[List[int]]) -> List[Dict]:
        out = []
        for idx in window_indices:
            df = self._degenerate_fraction(indices=idx)
            out.append(
                {
                    "indices": idx,
                    "degenerate_fraction": df,
                    "local_ME_exists": df >= self.consensus,
                }
            )
        return out

    def global_unification_test(self) -> Dict:
        df = self._degenerate_fraction(indices=None)
        return {
            "degenerate_fraction_all": df,
            "global_ME_exists": df >= self.consensus,
        }


def build_power_logs(
    primes: List[int], k: int, eps: float, label: str
) -> List[Tuple[float, float, str]]:
    out = []
    for p in primes:
        E = float(p**k)
        out.append((log(E), eps, f"{int(E)}={p}^{k} [{label}]"))
    return out


def run_demo(
    primes_square: List[int] = PRIMES_SQUARE,
    primes_cube: List[int] = PRIMES_CUBE,
) -> Dict:
    tester = MEQTesterV3()

    sq_obs = build_power_logs(primes_square, k=2, eps=EPS_SQ, label="SQ")
    for y, eps, tag in sq_obs:
        tester.add_log_obs(y, eps, tag)

    cb_obs = build_power_logs(primes_cube, k=3, eps=EPS_CB, label="CB")
    offset = len(sq_obs)
    for y, eps, tag in cb_obs:
        tester.add_log_obs(y, eps, tag)

    sq_idx = list(range(len(sq_obs)))
    cb_idx = list(range(offset, offset + len(cb_obs)))

    local_results = tester.local_tests([sq_idx, cb_idx])
    global_result = tester.global_unification_test()
    non_unif = (
        sum(1 for r in local_results if r["local_ME_exists"]) > 0
        and not global_result["global_ME_exists"]
    )

    return {
        "params": {
            "tau": tester.tau,
            "consensus": tester.consensus,
            "bootstrap_N": tester.bootstrap_N,
            "seed": tester.seed,
            "eps_sq": EPS_SQ,
            "eps_cb": EPS_CB,
            "primes_square": primes_square,
            "primes_cube": primes_cube,
        },
        "local": local_results,
        "global": global_result,
        "non_unification": non_unif,
    }


def main() -> None:
    out = run_demo()
    print("ME locals quadrática/cúbica — no-unificación")
    for i, r in enumerate(out["local"], 1):
        print(
            f"  Finestra {i}: deg_frac={r['degenerate_fraction']:.3f} "
            f"ME_local={r['local_ME_exists']}"
        )
    g = out["global"]
    print(
        f"  Global: deg_frac_all={g['degenerate_fraction_all']:.3f} "
        f"ME_global={g['global_ME_exists']}"
    )
    print(f"  NO-UNIFICACIÓN: {out['non_unification']}")
    print("\nJSON_RESULT =")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()