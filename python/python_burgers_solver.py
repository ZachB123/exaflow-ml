import numpy as np
import matplotlib.pyplot as plt

# Parameters (should match C++ values)
nx = 101
nt = 2000
dx = 2.0 / (nx - 1)
nu = 0.01
dt = 0.001
x = np.linspace(0, 2, nx)

# Initial condition (step function)
u = np.ones(nx)
u[int(50):int(100)] = 2.0  # velocity step

# Time-stepping loop (finite difference scheme)
for n in range(nt):
    u_new = u.copy()
    for i in range(1, nx - 1):
        u_new[i] = u[i] - u[i] * dt / dx * (u[i] - u[i - 1]) + nu * dt / (dx * dx) * (u[i + 1] - 2 * u[i] + u[i - 1])
    u = u_new

# Save the Python output for comparison (optional)
np.savetxt('burgers_output_python.dat', np.column_stack((x, u)))

# Plotting the results
plt.plot(x, u, label="Python")