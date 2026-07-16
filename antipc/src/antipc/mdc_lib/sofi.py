"""
Estructura Sofí VMA — L1/L3/L4/L2/U2 y candidatos Sophie Germain.

Filtrado desde filesclaude 6-5/sofi_structure.py y teoremas 04–05 (Libro 3).
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def is_type_a_composite(a: int) -> bool:
    if a < 25:
        return False
    limit = int(math.sqrt(a)) + 1
    f = 5
    while f <= limit:
        if a % f == 0 and (a // f) % 6 == 1:
            return True
        f += 6
    return False


def is_type_b_composite(a: int) -> bool:
    return is_type_a_composite(2 * a + 1)


def is_prime_trial(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True


@dataclass
class SofiResult:
    limit: int
    l1: int
    l3: int
    l4: int
    l2: int
    u2: int
    sg_in_u2: int
    prime_not_sg: int
    not_prime_in_u2: list[int]
    u2_subset_holds: bool
    u2_sample: list[int]
    sg_sample: list[int]


def classify_sofi(limit: int, *, verify: bool = True) -> SofiResult:
    l1_list = list(range(5, limit + 1, 6))
    l3: set[int] = set()
    l4: set[int] = set()
    for a in l1_list:
        if is_type_a_composite(a):
            l3.add(a)
        if is_type_b_composite(a):
            l4.add(a)
    l2 = l3 & l4
    u2 = set(l1_list) - l3 - l4

    not_prime: list[int] = []
    sg: list[int] = []
    prime_not_sg = 0
    if verify:
        for a in u2:
            if not is_prime_trial(a):
                not_prime.append(a)
            elif is_prime_trial(2 * a + 1):
                sg.append(a)
            else:
                prime_not_sg += 1

    return SofiResult(
        limit=limit,
        l1=len(l1_list),
        l3=len(l3),
        l4=len(l4),
        l2=len(l2),
        u2=len(u2),
        sg_in_u2=len(sg),
        prime_not_sg=prime_not_sg,
        not_prime_in_u2=not_prime,
        u2_subset_holds=len(not_prime) == 0,
        u2_sample=sorted(u2)[:12],
        sg_sample=sorted(sg)[:12],
    )


def format_sofi(r: SofiResult) -> str:
    lines = [
        f"Sofí VMA — límite {r.limit}",
        f"  |L1| = {r.l1:6d}  (6k−1 candidatos)",
        f"  |L3| = {r.l3:6d}  (compuestos tipo A)",
        f"  |L4| = {r.l4:6d}  (compuestos tipo B)",
        f"  |L2| = {r.l2:6d}  (doble compuesto)",
        f"  |U2| = {r.u2:6d}  (residual SGP)",
        f"  SGP en U2 : {r.sg_in_u2}",
        f"  U2 ⊆ LSG  : {'SÍ' if r.u2_subset_holds else 'NO'}",
    ]
    if r.sg_sample:
        lines.append(f"  Muestra SGP: {r.sg_sample}")
    if r.not_prime_in_u2:
        lines.append(f"  ⚠ no primos en U2: {r.not_prime_in_u2[:8]}")
    return "\n".join(lines)