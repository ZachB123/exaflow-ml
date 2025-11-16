# 1D Burgers' Equation Solver

A C++ implementation of a 1D Burgers' equation solver using finite difference methods. This solver can simulate both inviscid and viscous Burgers' equations with customizable initial conditions.

## What is Burgers' Equation?

Burgers' equation is a fundamental partial differential equation used in fluid mechanics that models how velocity fields evolve over time. It combines:
- **Convection**: How the wave shape moves and steepens
- **Diffusion**: How viscosity smooths out sharp features

This makes it useful for understanding shock waves, traffic flow, and basic fluid dynamics.

## Project Structure

```
.
├── include/
│   └── burgers.h          # Header file with BurgersSolver1d class
├── src/
│   ├── main.cpp           # Example simulations
│   └── burgers.cpp        # Solver implementation
├── visualizer/
│   └── visualize.py       # Display the burgers equations you have ran
├── CMakeLists.txt         # Build configuration
├── requirements.txt       # Python Libraries
└── data/                  # Output directory (created at runtime)

```

## Building the Project

### Prerequisites
- CMake 3.10 or higher
- C++17 compatible compiler (GCC, Clang, or MSVC)

### Build Instructions

1. Create a build directory and navigate to it:
```bash
mkdir build
cd build
```

2. Generate the build files with CMake:
```bash
cmake ..
```

3. Compile the project:
```bash
make
```

4. Run the executable:
```bash
./burgers
```

The program will generate output in the `../data` directory relative to the build folder.

## How to Use

### Creating a Simulation

To set up a simulation, you need two components:

#### 1. Solver Configuration

Define a `SolverConfig` struct with the following parameters:

```cpp
SolverConfig config = {
    .kinematic_viscosity = 0.01,    // Viscosity coefficient (ν)
    .num_domain_points = 201,       // Number of spatial grid points
    .time_steps = 2000,             // Number of time iterations
    .domain_length = 2.0,           // Length of spatial domain
    .time_step_size = 0.001         // Time step size (Δt)
};
```

**Parameter Explanation:**
- `kinematic_viscosity`: Controls diffusion/smoothing (higher = more smoothing)
- `num_domain_points`: Spatial resolution (more points = finer detail)
- `time_steps`: How long to simulate
- `domain_length`: The interval [0, L] over which the equation is solved
- `time_step_size`: Time between updates (smaller = more accurate but slower)

#### 2. Initial Condition Function

Define a lambda or function that specifies the initial velocity profile:

```cpp
// Example: Step function
std::function step_function = [](double x) -> double {
    if (x >= 0.5 && x <= 1.0) {
        return 2.0;
    } else {
        return 1.0;
    }
};

// Example: Sine wave
std::function sine_function = [](double x) -> double {
    return std::sin(x);
};
```

The function takes a position `x` and returns the initial velocity `u(x, t=0)` at that point.

#### 3. Create and Run the Solver

```cpp
// Create solver with configuration and initial conditions
BurgersSolver1d solver(config, step_function);

// Run the simulation
solver.solve();

// Save results (folder, run_name, gap)
solver.saveSolution("../data", "my_simulation", 10);
```

The `gap` parameter in `saveSolution` controls output frequency (e.g., `10` saves every 10th timestep).

## How the Solver Works

### The BurgersSolver1d Class

The solver uses an explicit finite difference scheme to discretize Burgers' equation:

**Burgers' Equation:**
```
∂u/∂t + u(∂u/∂x) = ν(∂²u/∂x²)
```

**Discretization:**
- **Time derivative**: Forward difference
- **Convection term**: Backward difference (upwind scheme)
- **Diffusion term**: Central difference

**Boundary Conditions:**
The solver implements periodic boundary conditions, meaning the solution wraps around at the domain edges (u[0] = u[N-1]).

### Key Methods

- `setInitialConditions()`: Initializes the velocity field u(x, 0)
- `solve()`: Iterates through time, computing u at each timestep
- `getSolution()`: Returns the complete solution history
- `saveSolution()`: Exports results to CSV files

## Output Structure

Results are saved in the following structure:

```
data/
└── [run_name]/
    ├── timestep_00000.csv
    ├── timestep_00010.csv
    ├── timestep_00020.csv
    └── ...
```

Each CSV file contains two columns:
```csv
x,u
0.0,1.0
0.01,1.0
0.02,1.0
...
```

- `x`: Spatial coordinate
- `u`: Velocity at that position

### Output Location

- The base output folder is `../data` relative to where you run the executable
- If running from the `build/` directory, data appears in `data/` at the project root
- Each simulation creates a subfolder named with the `run_name` parameter
- Existing folders with the same name are automatically overwritten

## Example Simulations

The included `main.cpp` demonstrates two classic test cases:

### 1. Step Function
- Initial condition: velocity jump from 1 to 2
- Shows shock formation and diffusion

### 2. Sine Wave
- Initial condition: sin(x) over [0, 2π]
- Demonstrates wave steepening and breaking

## Visualizer

The visualizer is written in Python using MatPlotLib. To run it we must first create a virtual environment, install the requirements and then execute visualize.py

```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
// speed is an optional parameter and will default to one
python visualize.py <folder name in data/> --speed 1
```


## Notes

- The solver uses periodic boundary conditions (suitable for wave problems)
- Stability depends on the CFL condition: ensure time steps are small enough
- For large viscosity or long time periods, you may need more timesteps
- The solver stores the complete history in memory (consider memory usage for large simulations)
