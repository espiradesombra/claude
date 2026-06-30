# anexoF_mrauv_calibrador.py
# Uso: python anexoF_mrauv_calibrador.py  n0=100000  dn=632
# Estima D_pred, primos en el tramo y L(n)-m(n) en tres puntos para calibrar V0,a0

import sys, math

E_MINUS_2 = math.e - 2.0

def L(n):
    return int(math.sqrt(n+3)) + 7

def K(n):
    return int(int(math.sqrt(n)) * 9/24)

# inversos factoriales hasta 2000
inv = [1.0]
for i in range(1,2001):
    inv.append(inv[-1]/i)

def m_sum(n):
    k = max(2, K(n))
    k = min(k, 2000)
    s = sum(inv[2:k+1])
    return math.sqrt(n+3) * s

if __name__=='__main__':
    # parse args
    n0 = 100000
    dn = int(2*math.sqrt(n0))
    for a in sys.argv[1:]:
        if a.startswith('n0='): n0=int(a.split('=',1)[1])
        if a.startswith('dn='): dn=int(a.split('=',1)[1])

    # tres puntos: n0, n0+dn, n0+2dn
    pts = [n0, n0+dn, n0+2*dn]
    med = []
    for x in pts:
        med.append(L(x) - m_sum(x))

    med0, med1, med2 = med
    V0 = 1.0/med0
    a0 = V0/med1
    j  = a0/med2 if med2!=0 else 0.0

    D0 = 1.0/max(1, int(math.log(max(3,n0)))) # marcador; sustituye por densidad real si se conoce
    D_pred = D0 + V0*dn + 0.5*a0*(dn**2)
    primes_est = int(round(D_pred * dn))

    print(f"n0={n0}, dn={dn}")
    print(f"med0={med0:.6f}, med1={med1:.6f}, med2={med2:.6f}")
    print(f"V0={V0:.6e}, a0={a0:.6e}, j={j:.6e}")
    print(f"D_pred={D_pred:.6f}, primes_range_est={primes_est}")
