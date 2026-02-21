def non_zero_viscosity_generator(solver):
    """
    Yield timestep, x index, and artificial viscosity
    only where viscosity is non-zero.
    """
    for t, u in enumerate(solver.solution_history):  # solver stores timesteps
        for x_idx, av in enumerate(solver.calculateArtificialViscosity(u)):
            if av != 0:
                yield t, x_idx, av
