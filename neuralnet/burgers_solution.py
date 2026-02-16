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
        
        metadata_path = os.path.join(self.sample_dir, METADATA_FILENAME)
        if not os.path.exists(metadata_path):
            raise ValueError(f"Metadata file not found: {metadata_path}")
        
        with open(metadata_path, 'r') as f:
            self.metadata = json.load(f)
        
        self.config = self.metadata[CONFIG_KEY]
        self.solver = self.metadata[SOLVER_KEY]
        
        self.domain_length = self.config[DOMAIN_LENGTH_KEY]
        self.spatial_step_size = self.solver[SPATIAL_STEP_SIZE_KEY]
        self.num_domain_points = self.solver[NUM_DOMAIN_POINTS_KEY]
        
        self.time_steps = self.solver[TIME_STEPS_KEY]
        self.time_step_size = self.solver[TIME_STEP_SIZE_KEY]
        self.max_time = self.time_steps * self.time_step_size
        
        self._cache = {}

        
    def _get_time_step(self, time_step_index):
        if time_step_index in self._cache:
            return self._cache[time_step_index]
        
        if time_step_index < 0 or time_step_index >= self.time_steps:
            raise ValueError(f"Time step index {time_step_index} out of bounds [0, {self.time_steps-1}]")
        
        # CSV files are named with 5 digits: _00000, _00001, etc.
        csv_filename = CSV_FILENAME_FORMAT.format(time_step_index)
        csv_path = os.path.join(self.sample_dir, csv_filename)
        
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        
        df = pd.read_csv(csv_path)
        x_array = df[X_COLUMN].values
        u_array = df[U_COLUMN].values
        
        self._cache[time_step_index] = (x_array, u_array)
        
        return x_array, u_array
    

    def _interpolate_spatial(self, x, x_array, u_array):
        return np.interp(x, x_array, u_array)
    

    def get_u(self, x, t):
        # Validate inputs
        if x < 0 or x > self.domain_length:
            raise ValueError(f"x={x} is out of domain bounds [0, {self.domain_length}]")
        if t < 0 or t > self.max_time:
            raise ValueError(f"t={t} is out of time bounds [0, {self.max_time}]")
        
        # Find surrounding time steps
        t_index_float = t / self.time_step_size
        t_index_lower = int(np.floor(t_index_float))
        t_index_upper = int(np.ceil(t_index_float))
        
        # Handle edge case where t is exactly at last time step
        if t_index_upper >= self.time_steps:
            t_index_upper = self.time_steps - 1
            t_index_lower = t_index_upper
        
        # Load data for both time steps
        x_array_lower, u_array_lower = self._get_time_step(t_index_lower)
        
        if t_index_lower == t_index_upper:
            # No temporal interpolation needed
            return self._interpolate_spatial(x, x_array_lower, u_array_lower)
        
        x_array_upper, u_array_upper = self._get_time_step(t_index_upper)
        
        u_lower = self._interpolate_spatial(x, x_array_lower, u_array_lower)
        u_upper = self._interpolate_spatial(x, x_array_upper, u_array_upper)
        
        t_lower = t_index_lower * self.time_step_size
        t_upper = t_index_upper * self.time_step_size
        
        if t_upper == t_lower:
            return u_lower
        
        weight = (t - t_lower) / (t_upper - t_lower)
        u_interpolated = u_lower * (1 - weight) + u_upper * weight
        
        return u_interpolated
    

    def clear_cache(self):
        self._cache.clear()
    

    def get_time_step_data(self, time_step_index):
        # basically gets like one csv file
        return self._get_time_step(time_step_index)
    

    def __repr__(self):
        return (f"BurgersSolution(sample='{self.sample_name}', "
                f"domain=[0, {self.domain_length:.2f}], "
                f"time=[0, {self.max_time:.4f}], "
                f"cached_steps={len(self._cache)})")