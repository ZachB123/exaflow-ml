from burgers_solution import BurgersSolution
import numpy as np

if __name__ == "__main__":
    solution = BurgersSolution("sample_000000")
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

    valid_count = 0
    none_count = 0

    for item in solution.requires_artificial_viscosity_generator():
        if item is None:
            none_count += 1
        else:
            valid_count += 1

    total = valid_count + none_count

    print(f"Valid yields: {valid_count}")
    print(f"None yields: {none_count}")
    print(f"Total: {total}")

    if total > 0:
        print(f"Valid ratio: {valid_count / total:.6f}")
        print(f"None ratio: {none_count / total:.6f}")
        
    # for abs_ux, dx, u_i, u_ip1, u_im1 in solution.requires_artificial_viscosity_generator():
    #     print(
    #         f"(abs_ux={abs_ux:.6f}, "
    #         f"dx={dx:.6f}, "
    #         f"u_i={u_i:.6f}, "
    #         f"u_ip1={u_ip1:.6f}, "
    #         f"u_im1={u_im1:.6f})"
    #     )