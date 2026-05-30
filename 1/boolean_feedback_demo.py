def xor(a, b):
    return tuple([ai ^ bi for ai, bi in zip(a, b)])

def shift(v1, v2):
    # Gray code shift as in the 1477 method
    gray = {(0,0):(1,1), (0,1):(1,0), (1,1):(0,0), (1,0):(0,1)}
    return gray[(v1, v2)]

def toffoli(v1, f1, v2, f2):
    return (v1, f1, v1 & v2, 0)

def sum_toffoli(s, K=3):
    v1, f1, v2, f2 = s
    states = [s]
    for _ in range(K-1):
        v1, f1, v2, f2 = shift(v1, v2)[0], f1, shift(v1, v2)[1], f2
        states.append((v1, f1, v2, f2))
    total = [0]*4
    for st in states:
        t = toffoli(*st)
        total = [x^y for x,y in zip(total, t)]
    return tuple(total)

def feedback_step(state):
    return xor(state, sum_toffoli(state))

# Ejemplo: dos estados que difieren solo en f1
SA = (0,0,1,0)
SB = (0,1,1,0)
cycles = 6

print("Ciclo\tEstadoA\tEstadoB\tDiferencias")
for i in range(cycles+1):
    dif = sum([a!=b for a,b in zip(SA,SB)])
    print(f"{i}\t{SA}\t{SB}\t{dif}")
    SA = feedback_step(SA)
    SB = feedback_step(SB)