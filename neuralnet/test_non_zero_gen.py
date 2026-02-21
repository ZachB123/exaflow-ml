from burgers_solution import BurgersSolution
from non_zero_viscosity_generator import non_zero_viscosity_generator

# Initialize solver with folder and example name
solver = BurgersSolution("example", "./training_data")

# Print all non-zero artificial viscosity points
for t, x, av in non_zero_viscosity_generator(solver):
    print(f"Timestep {t}, X index {x}, Artificial viscosity = {av}")
