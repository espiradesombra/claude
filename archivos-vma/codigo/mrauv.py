
import math

E_MINUS_2 = math.e - 2.0

inv = [1.0]
for i in range(1, 2001):
    inv.append(inv[-1] / i)

def L(n):
    return int(math.sqrt(n + 3)) + 7

def K(n):
    return int(int(math.sqrt(n)) * 9 / 24)

def m_sum(n):
    k = max(2, K(n))
    k = min(k, 2000)
    s = sum(inv[2:k + 1])
    return math.sqrt(n + 3) * s

def calibrar_mrauv(n0):
    dn = int(2 * math.sqrt(n0))
    pts = [n0, n0 + dn, n0 + 2 * dn]
    med = [L(x) - m_sum(x) for x in pts]
    med0, med1, med2 = med
    V0 = 1.0 / med0
    a0 = V0 / med1
    j = a0 / med2 if med2 != 0 else 0.0
    D0 = 1.0 / max(1, int(math.log(max(3, n0))))
    D_pred = D0 + V0 * dn + 0.5 * a0 * (dn ** 2)
    primes_est = int(round(D_pred * dn))
    return {
        "n0": n0,
        "dn": dn,
        "med": med,
        "V0": V0,
        "a0": a0,
        "j": j,
        "D_pred": D_pred,
        "primes_est": primes_est
    }
