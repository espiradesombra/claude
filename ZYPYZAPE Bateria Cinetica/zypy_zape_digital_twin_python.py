import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# ZYPYZAPE DIGITAL TWIN (PHYSICS UPGRADED)
# ==========================================

# --- PARAMETERS ---
N = 5
R = 60.0                     # rotor radius (m)
A = np.pi * R**2
rho = 1.225                  # air density
J = 5e6

Kd = 2e5                     # damping (control)
Kc = 5e5                     # coupling

omega_max = 2.5

# --- TIME ---
dt = 0.01
T = 60
steps = int(T/dt)

# --- STATE ---
omega = np.ones(N) * 1.2
theta = np.zeros(N)

omega_hist = []
power_aero_hist = []
power_elec_hist = []

# --- CONNECTIVITY ---
K = np.zeros((N, N))

# central pair
K[0,1] = K[1,0] = Kc

# ring
K[2,3] = K[3,2] = Kc
K[3,4] = K[4,3] = Kc
K[4,2] = K[2,4] = Kc

# --- WIND MODEL ---
def wind_speed(t):
    return 12 + 2*np.sin(0.1*t)   # simple variable wind

# --- Cp(lambda) MODEL ---
def Cp(lmbda):
    # simplified realistic curve
    return 0.22*(116/lmbda - 5)*np.exp(-12.5/lmbda)

# --- AERODYNAMIC TORQUE ---
def aero_torque(omega_i, v):
    if omega_i < 0.1:
        return 0
    lmbda = omega_i * R / v
    cp = max(0, Cp(lmbda))
    P = 0.5 * rho * A * v**3 * cp
    return P / omega_i

# --- ELECTRICAL CONTROL (ZYPYZAPE) ---
def zypyzape_control(i, t, omega_i, omega_cm):
    # base damping
    T = -Kd * (omega_i - omega_cm)

    # excitation on ring (ZypyZape cycle)
    if i >= 2:
        T += 5e5 * np.sin(2*np.pi*1*t)
    
    return T

# --- LOSSES ---
def friction(omega_i):
    return 1e4 * omega_i + 2e3 * omega_i**2

# --- SIMULATION ---
E_aero = 0
E_elec = 0

for step in range(steps):
    t = step * dt
    v = wind_speed(t)
    omega_cm = np.mean(omega)
    
    domega = np.zeros(N)
    
    for i in range(N):
        # --- aerodynamic torque ---
        T_aero = aero_torque(omega[i], v)
        P_aero = T_aero * omega[i]
        
        # --- coupling (more physical: phase-based) ---
        coupling = 0
        for j in range(N):
            if K[i,j] != 0:
                coupling += K[i,j] * np.sin(theta[j] - theta[i])
        
        # --- control torque ---
        T_ctrl = zypyzape_control(i, t, omega[i], omega_cm)
        P_elec = T_ctrl * omega[i]
        
        # --- friction ---
        T_fric = friction(omega[i])
        
        # --- dynamics ---
        domega[i] = (T_aero + coupling + T_ctrl - T_fric) / J
        
        # --- energy accounting ---
        E_aero += P_aero * dt
        E_elec += P_elec * dt
    
    omega += domega * dt
    theta += omega * dt
    
    omega = np.clip(omega, 0, omega_max)
    
    omega_hist.append(omega.copy())
    power_aero_hist.append(P_aero)
    power_elec_hist.append(P_elec)

omega_hist = np.array(omega_hist)

# ==========================================
# RESULTS
# ==========================================

print("\n--- ENERGY BALANCE ---")
print(f"Energy from wind: {E_aero:.2e} J")
print(f"Electrical energy: {E_elec:.2e} J")
print(f"Net: {E_elec - E_aero:.2e} J")

# --- PLOT ---
plt.figure()
for i in range(N):
    plt.plot(omega_hist[:, i], label=f'Turbine {i}')

plt.title('Rotor Speeds (Improved Twin)')
plt.xlabel('Time')
plt.ylabel('Omega')
plt.legend()
plt.grid()

plt.show()
