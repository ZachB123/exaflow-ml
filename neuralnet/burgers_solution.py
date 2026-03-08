import json
import os
import numpy as np
import pandas as pd
from constants import *


class BurgersSolution:

    def __init__(self, sample_name, training_data_dir=DEFAULT_TRAINING_DATA_DIR):

        self.sample_name = sample_name
        self.sample_dir = os.path.join(training_data_dir, sample_name)

        if not os.path.exists(self.sample_dir):
            raise ValueError(f"Sample directory does not exist: {self.sample_dir}")

        self._cache = {}

        # get solution files
        self.solution_bin_path = os.path.join(self.sample_dir, SOLUTION_DATA_FILENAME)
        self.metadata_path = os.path.join(self.sample_dir, METADATA_FILENAME)
        if not os.path.exists(self.solution_bin_path):
            raise ValueError(f"Binary solution file not found: {self.solution_bin_path}")
        if not os.path.exists(self.metadata_path):
            raise ValueError(f"Binary solution metadata file not found: {self.metadata_path}")

        # Parse solution metadata (JSON)
        with open(self.metadata_path, "r") as f:
            self.metadata = json.load(f)

        self.config = self.metadata[CONFIG_KEY]
        self.solver = self.metadata[SOLVER_KEY]

        self.domain_length = (self.solver[NUM_DOMAIN_POINTS_KEY] - 1) * self.solver[SPATIAL_STEP_SIZE_KEY]

        self.spatial_step_size = self.solver[SPATIAL_STEP_SIZE_KEY]
        self.num_domain_points = self.solver[NUM_DOMAIN_POINTS_KEY]
        self.time_steps = self.solver[TIME_STEPS_KEY]
        self.time_step_size = self.solver[TIME_STEP_SIZE_KEY]
        self.domain_length = self.config[DOMAIN_LENGTH_KEY]
        self.max_time = (self.time_steps - 1) * self.time_step_size

        self._cache = {}

        try:
            bias = float(self.metadata[BIAS_KEY])
            terms = self.metadata[TERMS_KEY]
        except KeyError as e:
            raise ValueError(
                f"Missing required metadata field: {e.args[0]}"
            )

        if not isinstance(terms, list) or len(terms) == 0:
            raise ValueError(
                f"'{TERMS_KEY}' must be a non-empty list to reconstruct u0(x)."
            )

        self.initial_conditions_params = {
            BIAS_KEY: bias,
            TERMS_KEY: terms,
            AMPLITUDE_KEY: np.array(
                [t[AMPLITUDE_KEY] for t in terms], dtype=float
            ),
            FREQUENCY_KEY: np.array(
                [t[FREQUENCY_KEY] for t in terms], dtype=float
            ),
            PHASE_SHIFT_KEY: np.array(
                [t[PHASE_SHIFT_KEY] for t in terms], dtype=float
            ),
        }

        self.time_steps = int(self.metadata[SOLVER_KEY][TIME_STEPS_KEY])
        self.time_step_size = float(self.metadata[SOLVER_KEY][TIME_STEP_SIZE_KEY])
        self.max_time = (self.time_steps - 1) * self.time_step_size

        self.spatial_step_size = float(self.metadata[SOLVER_KEY][SPATIAL_STEP_SIZE_KEY])
        self.num_domain_points = int(self.metadata[SOLVER_KEY][NUM_DOMAIN_POINTS_KEY])
        self.domain_length = float(self.num_domain_points - 1) * self.spatial_step_size

        # Memory-mapped 2D array; does not load everything into RAM, only time steps as requested
        self._u = np.memmap(
            self.solution_bin_path,
            dtype=np.float64,
            mode="r",
            shape=(self.time_steps, self.num_domain_points)
        )

        # x_array will be the same for all time steps, generate now to return with get_time_step_data()
        self._x_array = self.spatial_step_size * np.arange(self.num_domain_points)


    def initial_condition(self, x):
        if x < 0 or x > self.domain_length:
            raise ValueError(
                f"x={x} is out of domain bounds [0, {self.domain_length}]"
            )

        params = self.initial_conditions_params

        return float(
            params[BIAS_KEY]
            + np.sum(
                params[AMPLITUDE_KEY]
                * np.sin(
                    params[FREQUENCY_KEY]
                    * (x - params[PHASE_SHIFT_KEY])
                )
            )
        )


    def get_time_step(self, time_step_index):
        
        # Check cache first
        if time_step_index in self._cache:
            return self._cache[time_step_index]

        # Check bounds
        if time_step_index < 0 or time_step_index >= self.time_steps:
            raise ValueError(
                f"Time step index {time_step_index} out of bounds "
                f"[0, {self.time_steps - 1}]"
            )

        # Pull u from memmap
        u_array = self._u[time_step_index, :].copy()

        x_array = self._x_array
        
        self._cache[time_step_index] = (x_array, u_array)
        return x_array, u_array


    def requires_artificial_viscosity_generator(self):
        # lists the points that require artificial viscosity to make stable
        for t_index in range(self.time_steps):

            x_t, u_t = self.get_time_step(t_index)

            for i in range(self.num_domain_points - 1):
                ux = (
                    (u_t[i + 1] - u_t[self.num_domain_points - 2]) / (2.0 * self.spatial_step_size)
                    if i == 0
                    else (u_t[i + 1] - u_t[i - 1]) / (2.0 * self.spatial_step_size)
                )

                if ux < 0:
                    if i == 0:
                        yield (t_index, abs(ux), self.spatial_step_size, u_t[i], u_t[i + 1], u_t[self.num_domain_points - 2])
                    else:
                        yield (t_index, abs(ux), self.spatial_step_size, u_t[i], u_t[i + 1], u_t[i - 1])
                else:
                    yield None

    def _interpolate_spatial(self, x, x_array, u_array):
        return np.interp(x, x_array, u_array)


    def get_u(self, x, t):
        if x < 0 or x > self.domain_length:
            raise ValueError(
                f"x={x} is out of domain bounds [0, {self.domain_length}]"
            )

        if t < 0 or t > self.max_time:
            raise ValueError(
                f"t={t} is out of time bounds [0, {self.max_time}]"
            )

        t_index_float = t / self.time_step_size
        t_index_lower = int(np.floor(t_index_float))
        t_index_upper = int(np.ceil(t_index_float))

        if t_index_upper >= self.time_steps:
            t_index_upper = self.time_steps - 1
            t_index_lower = t_index_upper

        x_array_lower, u_array_lower = self.get_time_step(t_index_lower)

        if t_index_lower == t_index_upper:
            return self._interpolate_spatial(
                x, x_array_lower, u_array_lower
            )

        x_array_upper, u_array_upper = self.get_time_step(t_index_upper)

        u_lower = self._interpolate_spatial(
            x, x_array_lower, u_array_lower
        )

        u_upper = self._interpolate_spatial(
            x, x_array_upper, u_array_upper
        )

        t_lower = t_index_lower * self.time_step_size
        t_upper = t_index_upper * self.time_step_size

        if t_upper == t_lower:
            return u_lower

        weight = (t - t_lower) / (t_upper - t_lower)

        return u_lower * (1 - weight) + u_upper * weight


    def clear_cache(self):
        self._cache.clear()


    def __repr__(self):
        return (
            f"BurgersSolution(sample='{self.sample_name}', "
            f"domain=[0, {self.domain_length:.2f}], "
            f"time=[0, {self.max_time:.4f}], "
            f"num_points={self.num_domain_points})"
        )
