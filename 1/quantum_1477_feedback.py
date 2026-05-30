import numpy as np

# Tablas Gray-shift para dos bits (qubits)
def gray_shift(state):
    # input: tuple (v1, v2)
    mapping = {(0,0):(1,1), (0,1):(1,0), (1,1):(0,0), (1,0):(0,1)}
    return mapping[state]

def kron(*args):
    """Kronecker product for multiple matrices (for multi-qubit gates)."""
    result = args[0]
    for m in args[1:]:
        result = np.kron(result, m)
    return result

# Toffoli "cuántico" gate (3-qubit standard); aquí 2-qubits, le adaptamos para demo:
def quantum_toffoli():
    # Here, for demo: swap |11> <-> |10>, keep others
    U = np.eye(4, dtype=complex)
    # swap [2] (|10>) and [3] (|11>)
    U[2,2] = 0
    U[2,3] = 1
    U[3,2] = 1
    U[3,3] = 0
    # This is NOT the true Toffoli, but serves as a 2-Q generalization
    return U

def apply_shift(psi):
    """Permute the amplitudes: shift basis positions as in Gray code on (v1,v2)."""
    perm = [3,2,0,1] # (00)->(11), (01)->(10), (10)->(00), (11)->(01)
    return psi[perm]

def feedback_1477_quantum(psi, K=3, gate=None):
    # Aplica Toffoli cuantica a K trayectorias desfasadas y suma (coherente)
    Kstates = [psi]
    temp = psi
    for _ in range(K-1):
        temp = apply_shift(temp)
        Kstates.append(temp)
    # Suma coherente de acciones de Toffoli sobre cada trayectoria desfasada
    toff = quantum_toffoli() if gate is None else gate
    snaps = [toff @ vec for vec in Kstates]
    # Suma lineal (coherente)
    return sum(snaps)

def normalize(psi):
    return psi / np.linalg.norm(psi)

def measure_prob(psi):
    return np.abs(psi)**2

# Estados iniciales: |01> y (|01> + |11>)/sqrt(2)  # Diferencia de fase relativa en el bit más significativo
psiA = np.zeros(4, dtype=complex); psiA[1] = 1.0  # |01>
psiB = np.zeros(4, dtype=complex); psiB[1] = 1/np.sqrt(2); psiB[3] = 1/np.sqrt(2)  # Superposición |01> y |11>

K = 3
cycles = 6

print("Ciclo | Estado       | Probabilidades             | Estado B         | Probabilidades B")
for i in range(cycles+1):
    probsA = measure_prob(psiA)
    probsB = measure_prob(psiB)
    print(f"{i:2d}    | {np.round(psiA,2)} | {np.round(probsA,2)} | {np.round(psiB,2)} | {np.round(probsB,2)}")
    psiA = feedback_1477_quantum(psiA, K)
    psiB = feedback_1477_quantum(psiB, K)
    psiA = normalize(psiA)
    psiB = normalize(psiB)