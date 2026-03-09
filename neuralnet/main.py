import numpy as np
from tqdm import tqdm
import joblib
from multiprocessing import Pool, cpu_count

from nn import train_model
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


def get_training_data_folder_names():
    return sorted([item.name for item in DEFAULT_TRAINING_DATA_DIR.iterdir() if item.is_dir() and (item / METADATA_FILENAME).exists()])


def requires_artificial_viscosity(dx, u_i_minus_1, u_i_plus_1):
    ux = (u_i_plus_1 - u_i_minus_1) / (2.0 * dx)
    return ux < 0


# def reverse_engineer_cq(dt, dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1, nu=0, eps=1e-3):
#     uxx = u_i_plus_1 - 2*u_i + u_i_minus_1
#     ux_diff = abs(u_i_plus_1 - u_i_minus_1)
#
#     if abs(uxx) < eps or ux_diff < eps:
#         return None
#
#     return (2 / (dx * ux_diff)) * ((dx**2 / (dt * uxx)) * (u_next_i - u_i + u_i * (dt/dx) * (u_i_plus_1 - u_i_minus_1)) - nu)


# def reverse_engineer_cq(dt, dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1, nu=0):
#     # todo store nu in the metadata
#     # sometimes got divide by 0
#     uxx = (u_i_plus_1 - 2 * u_i + u_i_minus_1)
#     if uxx == 0:
#         return None
#
#     ux_diff = abs(u_i_plus_1 - u_i_minus_1)
#     if ux_diff == 0:
#         return None
#
#     return (2 / (dx * ux_diff)) * ((dx**2 / (dt * uxx)) * (u_next_i - u_i + u_i * (dt/dx) * (u_i_plus_1 - u_i_minus_1)) - nu)

def reverse_engineer_cq(dt, dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1, nu=0):
    # todo store nu in the metadata
    # sometimes got divide by 0
    if (u_i_plus_1 - 2*u_i + u_i_minus_1) == 0:
        return None

    return (2 / (dx * abs(u_i_plus_1 - u_i_minus_1))) * ((dx**2 / (dt * (u_i_plus_1 - 2 * u_i + u_i_minus_1))) * (u_next_i - u_i + u_i * (dt/dx) * (u_i_plus_1 - u_i_minus_1)) - nu)

def get_feature_matrices_for_sample(sample_name):
    X_rows = []
    y_rows = []

    burgers_solution = BurgersSolution(sample_name)
    fine_dt = burgers_solution.time_step_size
    fine_dx = burgers_solution.spatial_step_size
    coarse_dt = fine_dt * COARSENESS_MULTIPLIER
    coarse_dx = fine_dx * COARSENESS_MULTIPLIER
    coarse_num_timesteps = int(burgers_solution.max_time // coarse_dt)
    coarse_num_domain_points = int(burgers_solution.domain_length // coarse_dx)

    # Collect data from multiple timesteps
    for time_step in range(min(coarse_num_timesteps, 50)):
        # just ignoring the boundary condition for now cuz I need to ask yuvi about it
        # and it probably won't make a difference for training
        for spatial_step in range(1, coarse_num_domain_points - 1):
            curr_time = time_step * coarse_dt
            next_time = (time_step + 1) * coarse_dt
            curr_x = spatial_step * coarse_dx
            prev_x = (spatial_step - 1) * coarse_dx
            next_x = (spatial_step + 1) * coarse_dx

            u_i = burgers_solution.get_u(curr_x, curr_time)
            u_next_i = burgers_solution.get_u(curr_x, next_time)
            u_i_minus_1 = burgers_solution.get_u(prev_x, curr_time)
            u_i_plus_1 = burgers_solution.get_u(next_x, curr_time)

            if requires_artificial_viscosity(coarse_dx, u_i_minus_1, u_i_plus_1):
                cq = reverse_engineer_cq(coarse_dt, coarse_dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1)
                if cq is not None and not np.isnan(cq) and not np.isinf(cq):
                    # cq = np.clip(cq, -100, 100)
                    # if abs(cq) > 1000000000:
                    #     print(cq)
                    # we can't use u_next_i as a feature because we will never have that when running the sim normally
                    X_rows.append([coarse_dt, coarse_dx, u_i, u_i_minus_1, u_i_plus_1])
                    y_rows.append(cq)

    X = np.array(X_rows)
    y = np.array(y_rows)
    return X, y


if __name__ == "__main__":
    samples = get_training_data_folder_names()

    print(f"Found {len(samples)} samples. Extracting features...")
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(get_feature_matrices_for_sample, samples)

    X_matrices, y_matrices = zip(*results)

    X = np.vstack([m for m in X_matrices if m.size > 0])
    y = np.concatenate([m for m in y_matrices if m.size > 0])

    model, X_scaler = train_model(X, y)

