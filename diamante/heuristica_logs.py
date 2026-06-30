
import math

def log_iterativo(a, b, max_iter=100, tol=1e-15):
    if a <= 0 or b <= 0 or b == 1:
        return float('nan'), 0, 'dominio'

    j = 1.0
    d1 = d2 = d3 = 0.0
    f1, f2, f3, f4, f5, f6 = 1.2, 3.4, -0.8, 1.9, -0.8, 1.9

    for k in range(max_iter):
        j_old = j
        j = (a + (a / b**j) - 1) / a
        d3, d2, d1 = d2, d1, j - j_old
        j *= (1 + f2 * abs(d1))
        denom = max(abs(d3 - d1), 1e-12)
        factor = f3 + f4 * (d3 - d2) / denom
        factor = max(min(factor, 3.0), -2.0)
        j = math.copysign(abs(j)**abs(factor), j)
        if abs(j - j_old) < tol:
            return j, k + 1, 'ok'
    return j, max_iter, 'no_converge'
