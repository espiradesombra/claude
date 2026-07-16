import math

E_MINUS_2 = math.e - 2.0

def L(n):
    return int(math.sqrt(n+3)) + 7

def K(n):
    return int(int(math.sqrt(n)) * 9/24)

# precompute factorial inverses up to 200
fac = [1.0]
for i in range(1,201):
    fac.append(fac[-1]/i)


def m(n):
    Kmax = max(2, K(n))
    Kmax = min(Kmax, 200)
    root = math.sqrt(n+3)
    s = 0.0
    for i in range(2, Kmax+1):
        s += root * fac[i]
    return s

if __name__ == '__main__':
    for n in [1000, 5000, 20000, 100000]:
        Ln = L(n)
        mn = m(n)
        print(n, Ln, f"{mn:.10f}", f"{(Ln-mn):.10f}")
