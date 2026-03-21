import numpy as np
from multiprocessing import Pool, cpu_count

from constants import *
from burgers_solution import BurgersSolution


def get_training_data_folder_names():
    return sorted([item.name for item in DEFAULT_TRAINING_DATA_DIR.iterdir() if item.is_dir() and (item / METADATA_FILENAME).exists()])


def requires_artificial_viscosity(dx, u_i_minus_1, u_i_plus_1):
    ux = (u_i_plus_1 - u_i_minus_1) / (2.0 * dx)
    return ux < 0


def reverse_engineer_cq(dt, dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1, nu=0):

    if not requires_artificial_viscosity(dx, u_i_minus_1, u_i_plus_1):
        return None

    numerator = (
        u_next_i - u_i
        + u_i * (dt / (2 * dx)) * (u_i_plus_1 - u_i_minus_1)
        - nu * (dt / dx**2) * (u_i_plus_1 - 2 * u_i + u_i_minus_1)
    )

    denominator = (
        abs((u_i_plus_1 - u_i_minus_1) / (2 * dx))
        * dt
        * (u_i_plus_1 - 2 * u_i + u_i_minus_1)
    )

    if abs(denominator) < CQ_DENOMINATOR_EPSILON:
        return None  # would explode

    cq = numerator / denominator

    if abs(cq) >= CQ_MAX_MAGNITUDE:
        return None

    return cq


def get_feature_matrices_for_sample(sample_name, seed):
    X_rows = []
    y_rows = []

    burgers_solution = BurgersSolution(sample_name)
    nu = burgers_solution.nu
    fine_dt = burgers_solution.time_step_size
    fine_dx = burgers_solution.spatial_step_size
    coarse_dx = fine_dx * COARSENESS_MULTIPLIER
    coarse_dt_advection = (ALPHA * coarse_dx) / burgers_solution.max_u
    coarse_dt_diffusion = (BETA * (coarse_dx**2)) / burgers_solution.nu
    coarse_dt = min(coarse_dt_advection, coarse_dt_diffusion)

    print(f"Fine: dt={fine_dt:.6e}, dx={fine_dx:.6e} | Coarse: dt={coarse_dt:.6e}, dx={coarse_dx:.6e}")

    coarse_num_timesteps = int(burgers_solution.max_time // coarse_dt)
    coarse_num_domain_points = int(burgers_solution.domain_length // coarse_dx)

    # Randomly sample timesteps, then from each timestep randomly sample spatial points
    rng = np.random.default_rng(seed)
    all_timesteps = np.arange(coarse_num_timesteps - 1)
    sampled_timesteps = rng.choice(all_timesteps, size=min(MAX_TIMESTEPS_PER_SOLUTION, len(all_timesteps)), replace=False)
    all_spatial_steps = np.arange(U_RADIUS, coarse_num_domain_points - U_RADIUS)

    for time_step in sampled_timesteps:
        sampled_spatial_steps = rng.choice(all_spatial_steps, size=min(MAX_SPATIAL_POINTS_PER_TIMESTEP, len(all_spatial_steps)), replace=False)
        for spatial_step in sampled_spatial_steps:
            curr_time = time_step * coarse_dt
            next_time = (time_step + 1) * coarse_dt
            curr_x = spatial_step * coarse_dx

            additional_u_features = []
            for offset in range(-U_RADIUS, U_RADIUS + 1):
                x = (spatial_step + offset) * coarse_dx
                additional_u_features.append(burgers_solution.get_u(x, curr_time))

            u_i = additional_u_features[U_RADIUS]
            u_i_minus_1 = additional_u_features[U_RADIUS - 1]
            u_i_plus_1 = additional_u_features[U_RADIUS + 1]
            u_next_i = burgers_solution.get_u(curr_x, next_time)

            if requires_artificial_viscosity(coarse_dx, u_i_minus_1, u_i_plus_1):
                cq = reverse_engineer_cq(coarse_dt, coarse_dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1, nu)
                if cq is not None and not np.isnan(cq) and not np.isinf(cq):
                    del_u_del_x  = (u_i_plus_1 - u_i_minus_1) / (2 * coarse_dx)
                    del_u_del_x_squared = (u_i_plus_1 - 2 * u_i + u_i_minus_1) / (coarse_dx**2)
                    # we can't use u_next_i as a feature because we will never have that when running the sim normally
                    X_rows.append([coarse_dt, coarse_dx, del_u_del_x, del_u_del_x_squared] + additional_u_features + [nu])
                    y_rows.append(cq)

    return np.array(X_rows), np.array(y_rows)


if __name__ == "__main__":
    samples = get_training_data_folder_names()

    print(f"Found {len(samples)} samples. Extracting features...")
    seeds = [SEED + i for i in range(len(samples))]
    with Pool(processes=cpu_count()) as pool:
        results = pool.starmap(get_feature_matrices_for_sample, zip(samples, seeds))

    X_matrices, y_matrices = zip(*results)

    X = np.vstack([m for m in X_matrices if m.size > 0])
    y = np.concatenate([m for m in y_matrices if m.size > 0])

    print(f"X shape: {X.shape} | y shape: {y.shape}")
    print(f"y mean: {y.mean():.4f}, std: {y.std():.4f}, min: {y.min():.4f}, max: {y.max():.4f}")

    np.savez(FEATURE_MATRICES_PATH, X=X, y=y)
    print(f"Saved to {FEATURE_MATRICES_PATH}")
