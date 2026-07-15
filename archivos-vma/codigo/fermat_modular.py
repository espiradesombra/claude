
def fermat(n):
    return 2**(2**n) + 1

def base_B(n):
    prod = 1
    for k in range(n + 1):
        prod *= fermat(k)
    return 2 * prod

def residuo_privilegiado(n):
    mods = [fermat(k) for k in range(n + 1)]
    from sympy.ntheory.modular import crt
    residues = [2] * len(mods)
    r, _ = crt(mods, residues)
    return r

def alineacion_modular(n):
    Bn = base_B(n)
    rn = residuo_privilegiado(n)
    return [(m, fermat(m) % Bn == rn % Bn) for m in range(n + 1, n + 6)]
