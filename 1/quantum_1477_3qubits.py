import numpy as np
import matplotlib.pyplot as plt

# --------- Toffoli 3-qubit (CCNOT) ---------
def toffoli_3qubit():
    U = np.eye(8, dtype=complex)
    # Intercambia |110> <-> |111>
    U[6,6] = 0
    U[6,7] = 1
    U[7,6] = 1
    U[7,7] = 0
    return U

# --------- Shifter Gray para 3 bits (valor) ---------
def gray_shift_3bits(statevec):
    # Shift Gray cyclical en los 3 bits menos significativos (estados |xyz>)
    # Puedes ajustar el shift para experimentar diferentes "desfasajes"
    # Mapea índices (|000>,|001>,|010>,|011>,|100>,|101>,|110>,|111>)
    perm = [7,6,4,5,1,0,2,3]  # ejemplo: 000->111, 001->110, ..., 110->010, 111->011
    return statevec[perm]

# --------- Feedback quantum-1477 para 3 qubits ---------
def feedback_1477_quantum(psi, K=3, gate=None):
    Kstates = [psi]
    temp = psi
    for _ in range(K-1):
        temp = gray_shift_3bits(temp)
        Kstates.append(temp)
    gate = toffoli_3qubit() if gate is None else gate
    snaps = [gate @ vec for vec in Kstates]
    res = sum(snaps)
    return res

def normalize(psi):
    norm = np.linalg.norm(psi)
    return psi/norm if norm != 0 else psi

def measure_prob(psi):
    return np.abs(psi)**2

# --------- Estados iniciales ---------
# Inspirados en la codificación 1477: |v1,f1,v2> (ejemplo)
# Estado base (sin fase extra)
psiA = np.zeros(8, dtype=complex); psiA[2] = 1.0  # |010>
# Estado igual pero con fase relativa en f1 (qubit 1)
psiB = np.zeros(8, dtype=complex); psiB[2] = 1/np.sqrt(2); psiB[6] = 1/np.sqrt(2)  # superposición |010> y |110>

cycles = 8
K = 3  # caminos explorados
history_a = [psiA.copy()]
history_b = [psiB.copy()]
pa, pb = psiA.copy(), psiB.copy()

# --------- Evolución y registro ---------
for _ in range(cycles):
    pa = feedback_1477_quantum(pa, K=K)
    pb = feedback_1477_quantum(pb, K=K)
    pa = normalize(pa)
    pb = normalize(pb)
    history_a.append(pa.copy())
    history_b.append(pb.copy())

# --------- Graficar ---------
def plot_3qubit_evolution(history_a, history_b):
    steps = range(len(history_a))
    fig, axes = plt.subplots(1, 2, figsize=(16,5))
    states_labels = [f'|{i:03b}⟩' for i in range(8)]
    for i, (hist, title) in enumerate([
        (history_a, 'A: |010⟩'),
        (history_b, 'B: (|010⟩+|110⟩)/√2')
    ]):
        probs = np.array([[np.abs(amp)**2 for amp in state] for state in hist])
        for j in range(8):
            axes[i].plot(steps, probs[:,j], label=states_labels[j])
        axes[i].set_title(title)
        axes[i].set_xlabel('Ciclo')
        axes[i].set_ylabel('Probabilidad')
        axes[i].set_ylim(0,1)
        axes[i].legend(ncol=2)
        axes[i].grid(True)
    plt.tight_layout()
    plt.show()

plot_3qubit_evolution(history_a, history_b)

# --------- Tabla resumen ----------
print("Ciclo | Estado A [|000⟩ ... |111⟩]        | Estado B [|000⟩ ... |111⟩]")
for i in range(len(history_a)):
    pa_probs = np.round(measure_prob(history_a[i]),2)
    pb_probs = np.round(measure_prob(history_b[i]),2)
    print(f"{i:2d}    | {pa_probs} | {pb_probs}")

# --------- Ejercicios sugeridos para 3 qubits ----------
print("\nEjercicios para experimentar:")
print("1. Cambia el estado inicial a otras posiciones/fases (ej. más de dos componentes en superposición).")
print("2. Cambia K en feedback_1477_quantum() y ve cómo varía la propagación de la 'fase'.")
print("3. Ajusta 'gray_shift_3bits' para probar otras permutaciones/caminos.")
print("4. Analiza con qué parámetro inicial la amplitud/probabilidad de un estado base se anula por interferencia.")