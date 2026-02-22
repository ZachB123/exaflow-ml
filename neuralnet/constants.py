from pathlib import Path

# make path relative to constants.py
FILE_DIR = Path(__file__).resolve().parent
REPO_ROOT = FILE_DIR.parent
DEFAULT_TRAINING_DATA_DIR = REPO_ROOT / "training_data"

METADATA_FILENAME = "metadata.json"
CSV_FILENAME_FORMAT = "timestep_{:05d}.csv"

CONFIG_KEY = "config"
SOLVER_KEY = "solver"

DOMAIN_LENGTH_KEY = "domain_length"

SPATIAL_STEP_SIZE_KEY = "spatial_step_size"
NUM_DOMAIN_POINTS_KEY = "num_domain_points"
TIME_STEPS_KEY = "time_steps"
TIME_STEP_SIZE_KEY = "time_step_size"

X_COLUMN = "x"
U_COLUMN = "u"

BIAS_KEY = "bias"
TERMS_KEY = "terms"
AMPLITUDE_KEY = "amplitude"
FREQUENCY_KEY = "frequency"
PHASE_SHIFT_KEY = "phase_shift"