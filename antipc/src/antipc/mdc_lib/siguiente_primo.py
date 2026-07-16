"""
Siguiente primo VMA — acumuladores + Karnaugh 3-paso (teorema 07).

Filtrado desde filesclaude 6-5/siguiente_primo.py.
"""

from __future__ import annotations

from dataclasses import dataclass


def _is_prime_ref(n: int) -> bool:
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


def siguiente_primo(inicio: int) -> int | None:
    """Próximo primo estrictamente mayor que `inicio` (inicio debe ser primo)."""
    n = inicio
    ny = n - 1
    m = n + 1
    t = ny * ny
    tt = n * n
    nt = m * m
    n += 1
    m += 1
    ny += 1
    antp1 = False
    ant2p1 = False
    antp2 = False
    safety = inicio * 20 + 1000

    for _ in range(safety):
        t3 = t % m
        nt2 = nt % n
        t2 = t % n
        tt2 = tt % n
        t1 = t % ny
        nt1 = nt % ny

        core1 = (t3 > 0) and (nt2 == 0)
        paso1 = core1 or ant2p1
        new_ant2p1 = core1 or antp1
        new_antp1 = core1
        paso3 = antp2 and (t1 > 0) and (nt1 + nt2 == 0)
        new_antp2 = paso1 and (t2 > 0) and (tt2 > 0) and (t3 == 0)

        if paso3:
            detected = n - 1
            if detected > inicio:
                return detected

        t *= ny
        tt *= n
        nt *= m
        ny += 1
        n += 1
        m += 1
        antp1 = new_antp1
        ant2p1 = new_ant2p1
        antp2 = new_antp2

    return None


def primeros_n_primos(count: int, start: int = 2) -> list[int]:
    if count < 1:
        return []
    out: list[int] = []
    current = start if start <= 2 else start
    if current <= 2:
        out.append(2)
        current = 2
    elif not _is_prime_ref(current):
        while not _is_prime_ref(current):
            current += 1
        out.append(current)
    else:
        out.append(current)

    while len(out) < count:
        nxt = siguiente_primo(current)
        if nxt is None:
            break
        out.append(nxt)
        current = nxt
    return out


@dataclass
class ValidateRow:
    p: int
    expected: int
    computed: int | None
    ok: bool


def validate_chain(n_primes: int = 50) -> tuple[bool, list[ValidateRow]]:
    ref: list[int] = []
    candidate = 2
    while len(ref) < n_primes + 1:
        if _is_prime_ref(candidate):
            ref.append(candidate)
        candidate += 1

    rows: list[ValidateRow] = []
    all_ok = True
    for i in range(min(n_primes, len(ref) - 1)):
        p = ref[i]
        expected = ref[i + 1]
        computed = siguiente_primo(p)
        ok = computed == expected
        if not ok:
            all_ok = False
        rows.append(ValidateRow(p=p, expected=expected, computed=computed, ok=ok))
    return all_ok, rows


def format_next(inicio: int, nxt: int | None) -> str:
    if nxt is None:
        return f"Siguiente primo — desde p={inicio}: no encontrado (límite seguridad)"
    return f"Siguiente primo — desde p={inicio}: {nxt}"


def format_validate(all_ok: bool, rows: list[ValidateRow]) -> str:
    lines = [
        f"Validación siguiente_primo — {len(rows)} saltos",
        "  [HEUR] Conjetura 07 — algoritmo Karnaugh; no garantía formal",
    ]
    for row in rows[:12]:
        mark = "OK" if row.ok else "FALLO"
        lines.append(f"  {row.p:>6} → {row.expected:>6}  alg={row.computed}  [{mark}]")
    if len(rows) > 12:
        lines.append(f"  ... ({len(rows) - 12} más)")
    lines.append(f"  Resultado: {'OK' if all_ok else 'FALLO'}")
    return "\n".join(lines)