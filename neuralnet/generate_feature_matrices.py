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
    """Vectorized version: accepts scalars or numpy arrays. Returns cq or array of cqs (None for invalid)."""
    
    # Compute ux for all points (vectorized)
    ux = (u_i_plus_1 - u_i_minus_1) / (2.0 * dx)
    
    # Only compute cq where ux < 0 (requires artificial viscosity)
    requires_av = ux < 0
    
    # Initialize cq array with None values (as nan for now)
    if isinstance(u_i, np.ndarray):
        cq_result = np.full_like(u_i, np.nan, dtype=float)
    else:
        # Scalar case
        if not requires_av:
            return None
        cq_result = None

    if isinstance(u_i, np.ndarray):
        # Vectorized path: u_i, u_next_i, etc. are arrays
        numerator = (
            u_next_i - u_i
            + u_i * (dt / (2 * dx)) * (u_i_plus_1 - u_i_minus_1)
            - nu * (dt / dx**2) * (u_i_plus_1 - 2 * u_i + u_i_minus_1)
        )

        denominator = (
            np.abs(ux) 
            * dt
            * (u_i_plus_1 - 2 * u_i + u_i_minus_1)
        )

        # Only compute where ux < 0
        valid_mask = (np.abs(denominator) >= CQ_DENOMINATOR_EPSILON) & requires_av
        
        cq_array = np.full_like(u_i, np.nan, dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            cq_array[valid_mask] = numerator[valid_mask] / denominator[valid_mask]
        
        # Apply magnitude threshold
        cq_result = cq_array.copy()
        cq_result[np.abs(cq_array) >= CQ_MAX_MAGNITUDE] = np.nan
        
        return cq_result
    else:
        # Scalar path (original logic)
        numerator = (
            u_next_i - u_i
            + u_i * (dt / (2 * dx)) * (u_i_plus_1 - u_i_minus_1)
            - nu * (dt / dx**2) * (u_i_plus_1 - 2 * u_i + u_i_minus_1)
        )

        denominator = (
            abs(ux)
            * dt
            * (u_i_plus_1 - 2 * u_i + u_i_minus_1)
        )

        if abs(denominator) < CQ_DENOMINATOR_EPSILON:
            return None

        cq = numerator / denominator

        if abs(cq) >= CQ_MAX_MAGNITUDE:
            return None

        return cq


def get_feature_matrices_for_sample(sample_name, seed):
    import time
    start_sample = time.time()
    
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

    # Vectorized loop: process all spatial samples for a timestep in one shot
    for time_step in sampled_timesteps:
        sampled_spatial_steps = rng.choice(all_spatial_steps, size=min(MAX_SPATIAL_POINTS_PER_TIMESTEP, len(all_spatial_steps)), replace=False)
        
        # Current and next time indices (in actual time, not grid indices)
        curr_time = time_step * coarse_dt
        next_time = (time_step + 1) * coarse_dt
        
        # Convert actual time to fine-grid indices (same logic as BurgersSolution.get_u)
        t_index_float = curr_time / fine_dt
        t_index_lower = int(np.floor(t_index_float))
        t_index_upper = int(np.ceil(t_index_float))
        if t_index_upper >= burgers_solution.time_steps:
            t_index_upper = burgers_solution.time_steps - 1
        
        t_next_index_float = next_time / fine_dt
        t_next_index_lower = int(np.floor(t_next_index_float))
        t_next_index_upper = int(np.ceil(t_next_index_float))
        if t_next_index_upper >= burgers_solution.time_steps:
            t_next_index_upper = burgers_solution.time_steps - 1
        
        # Build all x coordinates we need for this timestep (center + radius points)
        all_x_coords = []
        
        for spatial_step in sampled_spatial_steps:
            for offset in range(-U_RADIUS, U_RADIUS + 1):
                x = (spatial_step + offset) * coarse_dx
                all_x_coords.append(x)
        
        if len(all_x_coords) == 0:
            continue
        
        all_x_coords = np.array(all_x_coords)
        
        # Vectorized spatial interpolation for current and next times
        # Get u values at all x positions, with temporal interpolation
        u_curr_lower = np.interp(all_x_coords, burgers_solution._x_array, burgers_solution._u[t_index_lower, :])
        if t_index_lower == t_index_upper:
            u_curr_all = u_curr_lower
        else:
            u_curr_upper = np.interp(all_x_coords, burgers_solution._x_array, burgers_solution._u[t_index_upper, :])
            weight = t_index_float - t_index_lower
            u_curr_all = u_curr_lower * (1 - weight) + u_curr_upper * weight
        
        # Same for next time
        u_next_lower = np.interp(all_x_coords, burgers_solution._x_array, burgers_solution._u[t_next_index_lower, :])
        if t_next_index_lower == t_next_index_upper:
            u_next_all = u_next_lower
        else:
            u_next_upper = np.interp(all_x_coords, burgers_solution._x_array, burgers_solution._u[t_next_index_upper, :])
            weight_next = t_next_index_float - t_next_index_lower
            u_next_all = u_next_lower * (1 - weight_next) + u_next_upper * weight_next
        
        # Now reorganize u values by spatial step (group into center + radius offsets per step)
        num_points = len(sampled_spatial_steps)
        batch_X = []
        batch_y = []
        
        for pt_idx, spatial_step in enumerate(sampled_spatial_steps):
            # Extract the u values for this spatial step (center + radius points)
            start_idx = pt_idx * (2 * U_RADIUS + 1)
            end_idx = start_idx + (2 * U_RADIUS + 1)
            
            additional_u_features = u_curr_all[start_idx:end_idx]
            u_i = additional_u_features[U_RADIUS]
            u_i_minus_1 = additional_u_features[U_RADIUS - 1]
            u_i_plus_1 = additional_u_features[U_RADIUS + 1]
            
            # u_next_i for the center point
            u_next_i = u_next_all[start_idx + U_RADIUS]
            
            # Check if artificial viscosity is required
            ux = (u_i_plus_1 - u_i_minus_1) / (2.0 * coarse_dx)
            if ux < 0:
                # Vectorized cq computation (even for single point, use scalar -> array -> scalar for consistency)
                cq = reverse_engineer_cq(coarse_dt, coarse_dx, u_i, u_next_i, u_i_minus_1, u_i_plus_1, nu)
                
                if cq is not None and not np.isnan(cq) and not np.isinf(cq):
                    del_u_del_x = ux
                    del_u_del_x_squared = (u_i_plus_1 - 2 * u_i + u_i_minus_1) / (coarse_dx**2)
                    batch_X.append([coarse_dt, coarse_dx, del_u_del_x, del_u_del_x_squared] + list(additional_u_features) + [nu])
                    batch_y.append(cq)
        
        # Batch extend (one extend per timestep instead of one append per point)
        if batch_X:
            X_rows.extend(batch_X)
            y_rows.extend(batch_y)

    elapsed = time.time() - start_sample
    print(f"  {sample_name}: {len(y_rows)} valid cq points extracted in {elapsed:.2f}s")
    
    return np.array(X_rows), np.array(y_rows)


if __name__ == "__main__":
    import time
    start_total = time.time()
    
    samples = get_training_data_folder_names()

    print(f"Found {len(samples)} samples. Extracting features...")
    seeds = [SEED + i for i in range(len(samples))]
    
    # Use chunksize to reduce worker overhead (chunks tasks to workers instead of one task per worker)
    chunksize = max(1, len(samples) // (cpu_count() * 4))  # Aim for ~4x more chunks than workers
    print(f"Using chunksize={chunksize} for {cpu_count()} CPUs")
    
    with Pool(processes=cpu_count()) as pool:
        results = pool.starmap(get_feature_matrices_for_sample, zip(samples, seeds), chunksize=chunksize)

    X_matrices, y_matrices = zip(*results)

    X = np.vstack([m for m in X_matrices if m.size > 0])
    y = np.concatenate([m for m in y_matrices if m.size > 0])

    print(f"X shape: {X.shape} | y shape: {y.shape}")
    print(f"y mean: {y.mean():.4f}, std: {y.std():.4f}, min: {y.min():.4f}, max: {y.max():.4f}")

    np.savez(FEATURE_MATRICES_PATH, X=X, y=y)
    print(f"Saved to {FEATURE_MATRICES_PATH}")
    
    elapsed_total = time.time() - start_total
    print(f"=== Total time: {elapsed_total:.2f}s ===" )
