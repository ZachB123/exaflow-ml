import numpy as np

from constants import *
from burgers_solution import BurgersSolution


def verify_initial_condition_reconstruction(sample):
    solution = BurgersSolution(sample)
    print(solution)

    #testing u0() vs. csv data for timestep_00000
    x_csv, u_csv = solution.get_time_step(0)
    u_reconstructed = np.array([solution.initial_condition(x) for x in x_csv])

    max_abs_error = np.max(np.abs(u_reconstructed - u_csv))
    mean_abs_error = np.mean(np.abs(u_reconstructed - u_csv))

    for i in range(10):
        print(f"x={x_csv[i]:.2f}  u0={u_reconstructed[i]:.5f}  timestep_00000={u_csv[i]:.5f}")

    print("Max abs error:", max_abs_error)
    print("Mean abs error:", mean_abs_error)
