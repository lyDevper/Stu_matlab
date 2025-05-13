import numpy as np
import scipy.linalg

def motor_kalman_matrices(Jm, Bm, Km, Ke, R, L, Ts):
    """
    Compute continuous and discretized state-space matrices for a revolute DC motor Kalman filter model.
    Uses state: [theta, omega, m_L, i]
    
    Parameters:
    - Jm: Inertia (kg·m²)
    - Bm: Viscous friction (N·m·s/rad)
    - Km: Torque constant (N·m/A)
    - Ke: Back EMF constant (V·s/rad)
    - R: Armature resistance (Ohm)
    - L: Armature inductance (H)
    - Ts: Sampling time (s)
    
    Returns:
    - A (continuous), B (continuous), G (disturbance), C (output matrix)
    - Ad (discrete), Bd (discrete)
    """
    Bm_Jm = Bm / Jm
    inv_Jm = 1 / Jm
    Km_Jm = Km / Jm
    Ke_L = Ke / L
    R_L = R / L
    inv_L = 1 / L

    A = np.array([
        [0, 1, 0, 0],
        [0, -Bm_Jm, inv_Jm, Km_Jm],
        [0, 0, 0, 0],
        [0, -Ke_L, 0, -R_L]
    ])

    B = np.array([
        [0],
        [0],
        [0],
        [inv_L]
    ])

    G = np.array([
        [0],
        [0],
        [1],
        [0]
    ])

    C = np.array([
        [1, 0, 0, 0],
        [0, 0, 0, 1]
    ])

    # Discretization using matrix exponential
    n_states = A.shape[0]
    M = np.block([
        [A, B],
        [np.zeros((1, n_states + 1))]
    ])

    Md = scipy.linalg.expm(M * Ts)
    Ad = Md[:n_states, :n_states]
    Bd = Md[:n_states, n_states:n_states+1]

    return A, B, G, C, Ad, Bd

Jm = 0.0502
Bm = 0.1831
Km = 0.4795
Ke = 0.492
R = 0.983
L = 0.0254
Ts = 0.001  # 1 ms

A, B, G, C, Ad, Bd = motor_kalman_matrices(Jm, Bm, Km, Ke, R, L, Ts)
print("Continuous A matrix:\n", A)
print("Continuous B matrix:\n", B)
print("Ad matrix:\n", Ad)
print("Bd matrix:\n", Bd)
