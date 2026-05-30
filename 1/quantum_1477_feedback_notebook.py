import numpy as np
import matplotlib.pyplot as plt

# ----------- Definiciones de puertas y funciones base -----------

def gray_shift(state):
    mapping = {(0,0):(1,1), (0,1):(1,0), (1,1):(0,0), (1,0):(0,1)}
    return mapping[state]

def apply_shift(psi):
    # Shift amplitudes siguiendo Gray code en base 2 qubits (|00>,|01>,|10>,|11>)
    perm = [3,2,0,1]  # (00)->(11), (01)->(10), (10)->(00), (11)->(01)
    return psi[perm]

def quantum_toffoli():
    # "Toffoli" para 2 qubits (aquí simplemente un swap entre |10> y |11>)
    U = np.eye(4, dtype=complex)
    U[2,2] = 0
    U[2,3] = 1
    U[3,2] = 1
    U[3,3] = 0
    return U

def quantum_cnot():
    # CNOT clásico: flips target (v2) si v1==1
    U = np.eye(4, dtype=complex)
    U[2,2] = 0
    U[2,3] = 1
    U[3,2] = 1
    U[3,3] = 0
    return U

def feedback_1477_quantum(psi, K=3, gate=None):
    # Suma de K trayectorias desfasadas, cada una procesada con la puerta elegida
    Kstates = [psi]
    temp = psi
    for _ in range(K-1):
        temp = apply_shift(temp)
        Kstates.append(temp)
    toff = quantum_toffoli() if gate is None else gate
    snaps = [toff @ vec for vec in Kstates]
    summed = sum(snaps)
    return summed

def normalize(psi):
    return psi / np.linalg.norm(psi)

def measure_prob(psi):
    return np.abs(psi)**2

# ----------- Estados iniciales -----------

# A: |01>
psiA = np.zeros(4, dtype=complex); psiA[1] = 1.0
# B: superposición |01> + |11>
psiB = np.zeros(4, dtype=complex); psiB[1] = 1/np.sqrt(2); psiB[3] = 1/np.sqrt(2)

K = 3                  # Número de caminos desfasados
cycles = 8             # Número de ciclos de feedback
used_gate = None       # Usa None o quantum_toffoli() o quantum_cnot() para experimentar

history_a = [psiA.copy()]
history_b = [psiB.copy()]
pa, pb = psiA.copy(), psiB.copy()

# ----------- Evolución y registro de historial -----------

for _ in range(cycles):
    pa = feedback_1477_quantum(pa, K=K, gate=used_gate)
    pb = feedback_1477_quantum(pb, K=K, gate=used_gate)
    pa = normalize(pa)
    pb = normalize(pb)
    history_a.append(pa.copy())
    history_b.append(pb.copy())

# ----------- Visualización gráfica -----------

def plot_evolution(history_a, history_b):
    steps = range(len(history_a))
    fig, axes = plt.subplots(1, 2, figsize=(14,4))
    for i, (hist, title) in enumerate([
        (history_a, 'Estado A (|01⟩)'),
        (history_b, 'Estado B (|01⟩ + |11⟩)/√2')
        ]):
        probs = np.array([[np.abs(amp)**2 for amp in state] for state in hist])
        for j in range(4):
            axes[i].plot(steps, probs[:,j], label=f'|{j:02b}⟩')
        axes[i].set_title(title)
        axes[i].set_xlabel('Ciclo')
        axes[i].set_ylabel('Probabilidad')
        axes[i].set_ylim(0,1)
        axes[i].legend()
        axes[i].grid(True)
    plt.tight_layout()
    plt.show()

plot_evolution(history_a, history_b)

# ----------- Tabla de probabilidades ciclo a ciclo -----------

print("Ciclo | Estado A [|00⟩ |01⟩ |10⟩ |11⟩]   || Estado B [|00⟩ |01⟩ |10⟩ |11⟩]")
for i in range(len(history_a)):
    pa_probs = np.round(np.abs(history_a[i])**2,2)
    pb_probs = np.round(np.abs(history_b[i])**2,2)
    print(f"{i:2d}    | {pa_probs} || {pb_probs}")

# ----------- Ejercicios sugeridos -----------

print("\nEjercicios sugeridos:")
print("1. Cambia el número de caminos 'K' (por ej., K=2 o K=4) y observa las gráficas.")
print("2. Usa quantum_cnot() como puerta y compara la evolución con Toffoli.")
print("3. Prueba estados iniciales diferentes, por ejemplo superposición en ambos qubits.")
print("4. ¿Cuándo desaparece una probabilidad completamente (efecto de la interferencia)? Experimenta y explica.")

# ----------- Modo booleano comparativo (opcional) -------------

def boolean_feedback_cycle(state, K=3):
    # Forma compacta del XOR-feedback para el sistema (como antes)
    def shift_bool(v1, v2):  # Gray
        mapping = {(0,0):(1,1), (0,1):(1,0), (1,1):(0,0), (1,0):(0,1)}
        return mapping[(v1,v2)]
    def toffoli_bool(v1,f1,v2,f2):
        return (v1, f1, v1 & v2, 0)
    Kstates = [state]
    v1, f1, v2, f2 = state
    for _ in range(K-1):
        v1, v2 = shift_bool(v1, v2)
        Kstates.append((v1, f1, v2, f2))
    # XOR suma
    res = [0,0,0,0]
    for s in Kstates:
        t = toffoli_bool(*s)
        for i in range(4):
            res[i] ^= t[i]
    # Feedback XOR
    return tuple([si^ti for si,ti in zip(state, res)])

# Ejemplo comparativo booleano
print("\nMétodo booleano comparativo:")
statesA = [(0,0,1,0)]
statesB = [(0,1,1,0)]
for _ in range(cycles):
    statesA.append(boolean_feedback_cycle(statesA[-1], K=K))
    statesB.append(boolean_feedback_cycle(statesB[-1], K=K))

print("Ciclo | Estado A | Estado B | Diferencias")
for i in range(len(statesA)):
    diff = sum([a!=b for a,b in zip(statesA[i],statesB[i])])
    print(f"{i:2d}    | {statesA[i]} | {statesB[i]} | {diff}")

# -------------- FIN --------------