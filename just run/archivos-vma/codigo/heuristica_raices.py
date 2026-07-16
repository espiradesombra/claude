
import math

def raiz_iterativa(a, n, max_iter=100, tol=1e-15):
    if a < 0 or n == 0:
        return float('nan'), 0, 'dominio'

    x = max(a**(1/n), 1e-6)
    d1 = d2 = d3 = 0.0
    f2, f3, f4 = 3.4, -0.8, 1.9

    for k in range(max_iter):
        x_old = x
        x = (a + a / x**n - 1) / a
        d3, d2, d1 = d2, d1, x - x_old
        x *= (1 + f2 * abs(d1))
        denom = max(abs(d3 - d1), 1e-12)
        factor = f3 + f4 * (d3 - d2) / denom
        factor = max(min(factor, 3.0), -2.0)
        x = math.copysign(abs(x)**abs(factor), x)
        if abs(x - x_old) < tol:
            return x, k + 1, 'ok'
    return x, max_iter, 'no_converge'
