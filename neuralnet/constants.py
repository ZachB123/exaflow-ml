from pathlib import Path

# make path relative to constants.py
FILE_DIR = Path(__file__).resolve().parent
REPO_ROOT = FILE_DIR.parent
DEFAULT_TRAINING_DATA_DIR = REPO_ROOT / "training_data"
NEURAL_NET_DIR = REPO_ROOT / "neuralnet"

SOLUTION_DATA_FILENAME = "solution.bin"
METADATA_FILENAME = "metadata.json"

CSV_FILENAME_FORMAT = "timestep_{:05d}.csv"

CONFIG_KEY = "config"
SOLVER_KEY = "solver"

DOMAIN_LENGTH_KEY = "domain_length"

SPATIAL_STEP_SIZE_KEY = "spatial_step_size"
NUM_DOMAIN_POINTS_KEY = "num_domain_points"
TIME_STEPS_KEY = "time_steps"
TIME_STEP_SIZE_KEY = "time_step_size"
KINEMATIC_VISCOSITY = "kinematic_viscosity"
MAX_U = "max_u"

X_COLUMN = "x"
U_COLUMN = "u"

FEATURE_MATRICES_PATH = NEURAL_NET_DIR / "feature_matrices.npz"

BIAS_KEY = "bias"
TERMS_KEY = "terms"
AMPLITUDE_KEY = "amplitude"
FREQUENCY_KEY = "frequency"
PHASE_SHIFT_KEY = "phase_shift"

# from solver for dt
ALPHA = 0.4
BETA = 0.2

# how much coarser our nn solve will be, for creating dt from the original solve dt
# TODO make this a range so we don't rely on a single multiplier
COARSENESS_MULTIPLIER = 2
CQ_DENOMINATOR_EPSILON = 1e-8
CQ_MAX_MAGNITUDE = 10.0
# how many extra features to have around the center
# so like if 3 we will have our center and 3 more in each direction
U_RADIUS = 3

EPOCHS = 1000

