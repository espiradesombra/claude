
from sympy import isprime, primerange

def generar_sofi_conjuntos(limite):
    L1 = [a for a in range(5, limite + 1, 6)]
    L3 = [a for a in L1 if any(a == (6*k - 1)*(6*h + 1) for k in range(1, a) for h in range(1, a))]
    L4 = [a for a in L1 if any(2*a + 1 == (6*j - 1)*(6*g + 1) for j in range(1, a) for g in range(1, a))]
    L2 = list(set(L3) & set(L4))
    LSG = [p for p in L1 if isprime(p) and isprime(2*p + 1)]
    U2 = list(set(L1) - set(L3) - set(L4))
    return {
        'L1': L1,
        'L3': L3,
        'L4': L4,
        'L2': L2,
        'LSG': LSG,
        'U2': U2
    }
