
import math

def is_candidate(v):
    return v % 6 == 1 or v % 6 == 5

def criba_hibrida(N):
    P = [0] * (N + 1)
    for v in range(2, N + 1):
        if v in [2, 3] or is_candidate(v):
            P[v] = 1
    M = N // 2
    for p in range(5, M + 1):
        if P[p] != 1:
            continue
        j = p * p
        salto1 = 2 * p
        salto2 = 4 * p
        toggle = True
        while j <= M:
            P[j] = 0
            j += salto1 if toggle else salto2
            toggle = not toggle
    for p in range(5, int(math.sqrt(N)) + 1):
        if P[p] != 1:
            continue
        j = (N // p) * p
        while j > M:
            if is_candidate(j):
                P[j] = 0
            j -= p
    return P
