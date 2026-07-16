"""Teorema Amplificador de Fase K=3 XOR — verificación (THEOREM_PhaseAmplifier_K3_XOR.md)."""


def shift(state: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    v1, f1, v2, f2 = state
    m = {(0, 0): (1, 1), (0, 1): (1, 0), (1, 1): (0, 0), (1, 0): (0, 1)}
    v1n, v2n = m[(v1, v2)]
    return (v1n, v2n, f1, f2)


def toffoli(state: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    v1, f1, v2, f2 = state
    return (v1, f1, v1 & v2, 0)


def t3(state: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    def xor4(a, b):
        return tuple(a[i] ^ b[i] for i in range(4))

    return xor4(xor4(toffoli(state), toffoli(shift(state))), toffoli(shift(shift(state))))


def next_state(state: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    t = t3(state)
    return tuple(state[i] ^ t[i] for i in range(4))


def run_demo(steps: int = 6) -> list[tuple[int, tuple, tuple, int]]:
    s_a = (0, 0, 1, 0)
    s_b = (0, 1, 1, 0)
    rows = []
    for t in range(steps):
        dist = sum(a != b for a, b in zip(s_a, s_b))
        rows.append((t, s_a, s_b, dist))
        s_a, s_b = next_state(s_a), next_state(s_b)
    return rows