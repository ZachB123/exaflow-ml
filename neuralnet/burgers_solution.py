import json
import os
import numpy as np
import pandas as pd
from .constants import *


class BurgersSolution:

    def __init__(self, sample_name, training_data_dir=DEFAULT_TRAINING_DATA_DIR):

        # Store sample info and folder paths
        self.sample_name = sample_name
        self.sample_dir = os.path.join(training_data_dir, sample_name)

        if not os.path.exists(self.sample_dir):
            raise ValueError(f"Sample directory does not exist: {self.sample_dir}")

        # Load metadata.json
        metadata_path = os.path.join(self.sample_dir, METADATA_FILENAME)

        if not os.path.exists(metadata_path):
            raise ValueError(f"Metadata file not found: {metadata_path}")

        with open(metadata_path, 'r') as f:
            self.metadata = json.load(f)

        # Extract config and solver info using constants
        self.config = self.metadata[CONFIG_KEY]
        self.solver = self.metadata[SOLVER_KEY]

        self.domain_length = (
            (self.config[NUM_DOMAIN_POINTS_KEY] - 1)
            * self.config[SPATIAL_STEP_SIZE_KEY]
        )

        self.spatial_step_size = self.solver[SPATIAL_STEP_SIZE_KEY]
        self.num_domain_points = self.solver[NUM_DOMAIN_POINTS_KEY]
        self.time_steps = self.solver[TIME_STEPS_KEY]
        self.time_step_size = self.solver[TIME_STEP_SIZE_KEY]
        self.max_time = (self.time_steps - 1) * self.time_step_size

        # Cache for on-demand timestep loading
        self._cache = {}

        # Validate and store initial condition parameters
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

    def calculateArtificialViscosity(self, u):
        """Compute artificial viscosity for a given timestep array u."""
        av = np.zeros_like(u)

        for i in range(1, len(u) - 1):
            if abs(u[i + 1] - u[i - 1]) > 0.1:
                av[i] = 0.01 * abs(u[i + 1] - u[i - 1])

        return av

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

    def _get_time_step(self, time_step_index):

        if time_step_index in self._cache:
            return self._cache[time_step_index]

        if time_step_index < 0 or time_step_index >= self.time_steps:
            raise ValueError(
                f"Time step index {time_step_index} out of bounds "
                f"[0, {self.time_steps - 1}]"
            )

        csv_filename = CSV_FILENAME_FORMAT.format(time_step_index)
        csv_path = os.path.join(self.sample_dir, csv_filename)

        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV file not found: {csv_path}")

        df = pd.read_csv(csv_path)

        x_array = df[X_COLUMN].values
        u_array = df[U_COLUMN].values

        self._cache[time_step_index] = (x_array, u_array)

        return x_array, u_array

    def non_zero_viscosity_points(self):
        """
        Yield timestep, x index, and artificial viscosity
        only where viscosity is required (ux < 0)
        """

        for t_index in range(self.time_steps):

            _, u = self._get_time_step(t_index)

            av_array = self.calculateArtificialViscosity(u)

            for i in range(1, self.num_domain_points - 1):

                ux = (
                    (u[i + 1] - u[i - 1])
                    / (2.0 * self.spatial_step_size)
                )

                if ux < 0:
                    yield t_index, i, av_array[i]

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

        x_array_lower, u_array_lower = self._get_time_step(t_index_lower)

        if t_index_lower == t_index_upper:
            return self._interpolate_spatial(
                x, x_array_lower, u_array_lower
            )

        x_array_upper, u_array_upper = self._get_time_step(t_index_upper)

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

    def get_time_step_data(self, time_step_index):
        return self._get_time_step(time_step_index)

    def __repr__(self):
        return (
            f"BurgersSolution(sample='{self.sample_name}', "
            f"domain=[0, {self.domain_length:.2f}], "
            f"time=[0, {self.max_time:.4f}], "
            f"num_points={self.num_domain_points})"
        )
