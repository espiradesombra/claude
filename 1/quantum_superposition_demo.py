import numpy as np

# Estados base
ket0 = np.array([1,0])
ket1 = np.array([0,1])

# Estado clásico: solo |0>
psi_a = ket0.copy()
psi_b = (ket0 + ket1)/np.sqrt(2)  # Superposición

def hadamard():
    return (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]])

def show_state(psi, label):
    print(f"{label}: {psi}  Probabilidad 0: {np.abs(psi[0])**2:.2f}  1: {np.abs(psi[1])**2:.2f}")

# Evolución
H = hadamard()
print("\nEvolución cuántica con Hadamard:")
for step in range(3):
    psi_a = H @ psi_a
    psi_b = H @ psi_b
    show_state(psi_a, f"Paso {step+1} - Estado A (partía de |0>)")
    show_state(psi_b, f"Paso {step+1} - Estado B (partía de superposición)")