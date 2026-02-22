from burgers_solution import BurgersSolution
import numpy as np

if __name__ == "__main__":
    solution = BurgersSolution("sample_000000")
    print(solution)
    # print(solution.get_u(10, 0.02 + (1e-5 / 2)))

    # print(solution.get_u(10, 0))
    # print(solution.initial_condition(10))

    #testing u0() vs. csv data for timestep_00000
    x_csv, u_csv = solution.get_time_step_data(0)
    u_reconstructed = np.array([solution.initial_condition(x) for x in x_csv])

    max_abs_error = np.max(np.abs(u_reconstructed - u_csv))
    mean_abs_error = np.mean(np.abs(u_reconstructed - u_csv))

    for i in range(10):
        print(f"x={x_csv[i]:.2f}  u0={u_reconstructed[i]:.5f}  timestep_00000={u_csv[i]:.5f}")
    
    print("Max abs error:", max_abs_error)
    print("Mean abs error:", mean_abs_error)