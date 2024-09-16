import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define constants
g = 9.81  # gravity (m/s^2)

def equations_of_motion(t, y, N, lengths, masses):
    # y contains [theta_1, theta_2, ..., theta_N, theta_dot_1, theta_dot_2, ..., theta_dot_N]
    thetas = y[:N]
    theta_dots = y[N:]
    
    theta_ddots = np.zeros(N)  # To store accelerations
    
    # Mass matrix and forcing terms
    M = np.zeros((N, N))  # Mass matrix
    F = np.zeros(N)       # Forcing vector
    
    # Fill M and F based on the equations of motion (Yesilyurt's method)
    for j in range(N):
        sigma_jk_sum = 0
        for k in range(j, N):
            sigma_jk_sum += masses[k]
        
        M[j, j] = (sum(masses[j:]) * lengths[j] ** 2)
        
        for k in range(j + 1, N):
            cos_term = np.cos(thetas[j] - thetas[k])
            M[j, k] = masses[k] * lengths[j] * lengths[k] * cos_term
            M[k, j] = M[j, k]
        
        F[j] = -sigma_jk_sum * g * lengths[j] * np.sin(thetas[j])
        
        # Add nonlinear interaction terms
        for k in range(N):
            if j != k:
                F[j] += masses[k] * lengths[j] * lengths[k] * np.sin(thetas[j] - thetas[k]) * theta_dots[j] * theta_dots[k]
    
    # Solve for angular accelerations
    theta_ddots = np.linalg.solve(M, F)
    
    return np.concatenate([theta_dots, theta_ddots])

# Function to simulate the N-point pendulum
def simulate_pendulum(N, lengths, masses, initial_conditions, t_span, t_eval):
    y0 = np.concatenate([initial_conditions["thetas"], initial_conditions["theta_dots"]])
    sol = solve_ivp(equations_of_motion, t_span, y0, args=(N, lengths, masses), t_eval=t_eval)
    return sol

# Example setup for a 3-point pendulum
N = 3
lengths = [1.0, 1.0, 1.0]  # lengths of the rods
masses = [1.0, 1.0, 1.0]   # masses of the pendulum bobs
initial_conditions = {
    "thetas": [np.pi / 4, np.pi / 6, np.pi / 3],  # initial angles (in radians)
    "theta_dots": [0.0, 0.0, 0.0]                # initial angular velocities
}
t_span = (0, 10)  # time span for the simulation
t_eval = np.linspace(0, 10, 1000)  # time points to evaluate

# Run the simulation
sol = simulate_pendulum(N, lengths, masses, initial_conditions, t_span, t_eval)

# Plot the results
plt.figure(figsize=(10, 6))
for i in range(N):
    plt.plot(sol.t, sol.y[i, :], label=f'Theta {i+1}')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.legend()
plt.title(f'N-point Pendulum Simulation for N={N}')
plt.show()
