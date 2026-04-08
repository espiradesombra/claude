import numpy as np
import matplotlib.pyplot as plt

# ==============================
# ZYPYZAPE DIGITAL TWIN (SIMPLE)
# ==============================

# Number of turbines (5 module)
N = 5

# Physical parameters
J = 5e6              # inertia
Kd = 2e5             # damping gain
Kc = 5e5             # coupling strength

# Time params
dt = 0.01
T = 20
steps = int(T/dt)

# State variables
omega = np.ones(N) * 1.5   # initial speeds
omega_hist = []

# Connectivity (ring + central pair)
K = np.zeros((N, N))

# central pair (0,1)
K[0,1] = K[1,0] = Kc

# ring (2,3,4)
K[2,3] = K[3,2] = Kc
K[3,4] = K[4,3] = Kc
K[4,2] = K[2,4] = Kc

# external excitation (only ring)
def external_torque(t):
    T_ext = np.zeros(N)
    T_ext[2] = 1e6 * np.sin(2*np.pi*0.5*t)
    T_ext[3] = -1e6 * np.sin(2*np.pi*0.5*t)
    return T_ext

# simulation loop
for step in range(steps):
    t = step * dt
    omega_cm = np.mean(omega)
    
    domega = np.zeros(N)
    
    for i in range(N):
        # coupling term
        coupling = 0
        for j in range(N):
            coupling += K[i,j] * (omega[j] - omega[i])
        
        # damping (your key idea)
        damping = -Kd * (omega[i] - omega_cm)
        
        # external torque
        T_ext = external_torque(t)[i]
        
        # dynamics
        domega[i] = (coupling + damping + T_ext) / J
    
    omega += domega * dt
    omega_hist.append(omega.copy())

omega_hist = np.array(omega_hist)

# ==============================
# PLOT RESULTS
# ==============================

plt.figure()
for i in range(N):
    plt.plot(omega_hist[:, i], label=f'Turbine {i}')

plt.title('ZypyZape - Rotor Speeds')
plt.xlabel('Time step')
plt.ylabel('Omega')
plt.legend()
plt.grid()

plt.show()
