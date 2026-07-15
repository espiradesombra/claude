
def siguiente_primo_heuristica(p0):
    n = p0
    m = n + 1
    ny = n - 1
    t = ny**2
    tt = n**2
    nt = m**2
    restos = [
        t % ny, t % n, t % m,
        tt % ny, tt % n, tt % m,
        nt % ny, nt % n, nt % m
    ]
    bools = [r == 0 for r in restos]
    paso1 = bools[0] or bools[3]
    paso2 = bools[1] and not bools[4]
    paso3 = bools[2] ^ bools[5]
    return paso1 and paso2 and paso3
