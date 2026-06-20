import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 4 qubits: |v1, f1, v2, f2> -> índice = v1*8 + f1*4 + v2*2 + f2
# ------------------------------------------------------------

def index(v1, f1, v2, f2):
    return (v1 << 3) | (f1 << 2) | (v2 << 1) | f2

def basis_state(v1, f1, v2, f2):
    vec = np.zeros(16, dtype=complex)
    vec[index(v1, f1, v2, f2)] = 1.0
    return vec

def gray_shift(v1, v2):
    """Gray code shift as in 1477 method"""
    mapping = {(0,0):(1,1), (0,1):(1,0), (1,1):(0,0), (1,0):(0,1)}
    return mapping[(v1, v2)]

def apply_gray_shift_to_state(psi):
    """Permute amplitudes according to Gray shift on (v1, v2)"""
    new_psi = np.zeros(16, dtype=complex)
    for v1 in (0,1):
        for f1 in (0,1):
            for v2 in (0,1):
                for f2 in (0,1):
                    v1p, v2p = gray_shift(v1, v2)
                    old_idx = index(v1, f1, v2, f2)
                    new_idx = index(v1p, f1, v2p, f2)
                    new_psi[new_idx] += psi[old_idx]
    return new_psi

def apply_toffoli(psi):
    """Classical Toffoli: (v1, f1, v2, f2) -> (v1, f1, v1&v2, 0)"""
    new_psi = np.zeros(16, dtype=complex)
    for v1 in (0,1):
        for f1 in (0,1):
            for v2 in (0,1):
                for f2 in (0,1):
                    v1_out = v1
                    f1_out = f1
                    v2_out = v1 & v2
                    f2_out = 0
                    out_idx = index(v1_out, f1_out, v2_out, f2_out)
                    in_idx = index(v1, f1, v2, f2)
                    new_psi[out_idx] += psi[in_idx]
    return new_psi

def feedback_1477_quantum(psi, K=3):
    """Quantum coherent feedback: sum over K shifted paths, each Toffoli'd."""
    paths = [psi]
    current = psi
    for _ in range(K-1):
        current = apply_gray_shift_to_state(current)
        paths.append(current)
    total = np.zeros(16, dtype=complex)
    for p in paths:
        total += apply_toffoli(p)
    return total

def normalize(psi):
    norm = np.linalg.norm(psi)
    if norm == 0:
        return psi
    return psi / norm

def measure_prob(psi):
    return np.abs(psi)**2

# ------------------------------------------------------------
# Estados iniciales: solo difieren en f1 (segundo bit)
# Estado A: (0,0,1,0) -> |0,0,1,0>  índice 2
# Estado B: (0,1,1,0) -> |0,1,1,0>  índice 6
# ------------------------------------------------------------
psiA = basis_state(0,0,1,0)   # |0,0,1,0>
psiB = basis_state(0,1,1,0)   # |0,1,1,0>   # diferencia en f1

K = 3
cycles = 8
history_a = [psiA.copy()]
history_b = [psiB.copy()]
pa, pb = psiA.copy(), psiB.copy()

for _ in range(cycles):
    pa = feedback_1477_quantum(pa, K)
    pb = feedback_1477_quantum(pb, K)
    pa = normalize(pa)
    pb = normalize(pb)
    history_a.append(pa.copy())
    history_b.append(pb.copy())

# ------------------------------------------------------------
# Visualización
# ------------------------------------------------------------
def plot_4qubit_evolution(history_a, history_b):
    steps = range(len(history_a))
    labels = [f'|{v1}{f1}{v2}{f2}⟩' for v1 in (0,1) for f1 in (0,1) for v2 in (0,1) for f2 in (0,1)]
    fig, axes = plt.subplots(1, 2, figsize=(16,5))
    for i, (hist, title) in enumerate([(history_a, 'A: |0,0,1,0⟩'), (history_b, 'B: |0,1,1,0⟩')]):
        probs = np.array([measure_prob(state) for state in hist])
        for j in range(16):
            axes[i].plot(steps, probs[:,j], label=labels[j])
        axes[i].set_title(title)
        axes[i].set_xlabel('Ciclo')
        axes[i].set_ylabel('Probabilidad')
        axes[i].set_ylim(0,1)
        axes[i].legend(ncol=2, fontsize=6)
        axes[i].grid(True)
    plt.tight_layout()
    plt.show()

plot_4qubit_evolution(history_a, history_b)

# ------------------------------------------------------------
# Tabla de probabilidades (solo bases con prob > 0.01)
# ------------------------------------------------------------
print("Ciclo | Estado A (probabilidades)                 | Estado B")
for i in range(len(history_a)):
    probs_a = measure_prob(history_a[i])
    probs_b = measure_prob(history_b[i])
    nonzero_a = [(f'|{v1}{f1}{v2}{f2}⟩', probs_a[idx]) 
                 for v1 in (0,1) for f1 in (0,1) for v2 in (0,1) for f2 in (0,1)
                 if probs_a[idx] > 0.01]
    nonzero_b = [(f'|{v1}{f1}{v2}{f2}⟩', probs_b[idx]) 
                 for idx, (v1,f1,v2,f2) in enumerate([(a,b,c,d) for a in (0,1) for b in (0,1) for c in (0,1) for d in (0,1)])
                 if probs_b[idx] > 0.01]
    print(f"{i:2d} | {nonzero_a} | {nonzero_b}")

# ------------------------------------------------------------
# Comparación con el método booleano original (XOR feedback)
# ------------------------------------------------------------
def boolean_toffoli(v1, f1, v2, f2):
    return (v1, f1, v1 & v2, 0)

def boolean_gray_shift(v1, v2):
    mapping = {(0,0):(1,1), (0,1):(1,0), (1,1):(0,0), (1,0):(0,1)}
    return mapping[(v1, v2)]

def boolean_sum_toffoli(state, K=3):
    v1, f1, v2, f2 = state
    states = [state]
    for _ in range(K-1):
        v1, v2 = boolean_gray_shift(v1, v2)
        states.append((v1, f1, v2, f2))
    total = [0,0,0,0]
    for s in states:
        t = boolean_toffoli(*s)
        total = [total[i] ^ t[i] for i in range(4)]
    return tuple(total)

def boolean_feedback_step(state, K=3):
    return tuple([state[i] ^ boolean_sum_toffoli(state, K)[i] for i in range(4)])

# Ejecutar booleano
SA = (0,0,1,0)
SB = (0,1,1,0)
statesA = [SA]
statesB = [SB]
for _ in range(cycles):
    SA = boolean_feedback_step(SA, K)
    SB = boolean_feedback_step(SB, K)
    statesA.append(SA)
    statesB.append(SB)

print("\nMétodo booleano (XOR feedback):")
print("Ciclo | Estado A      | Estado B      | Diferencias")
for i in range(cycles+1):
    diff = sum(a != b for a,b in zip(statesA[i], statesB[i]))
    print(f"{i:2d}    | {statesA[i]} | {statesB[i]} | {diff}")