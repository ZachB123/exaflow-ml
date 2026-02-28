#include <cmath>
#include <fstream>
#include <functional>
#include <vector>

#include "burger_scheme.h"
#include "burgers.h"
#include "initial_condition_generator.h"

int main() {

  SolverConfig step_function_config = {.kinematic_viscosity = 0.01,
                                       .time_steps = 2000,
                                       .domain_length = 2.0,
                                       .time_step_size = 0.001};

  // 1 everywhere in domain
  // for 0.5 <= x <= 1 the function value is 2
  std::function<double(double)> step_function = [](double x) -> double {
    if (x >= 0.5 && x <= 1.0) {
      return 2.0;
    } else {
      return 1.0;
    }
  };

  BurgersSolver1d step_function_solver(std::make_unique<Godunov>(),
                                       step_function_config, step_function);

  step_function_solver.solve();
  step_function_solver.saveSolution("../data", "step_function", 1);
  std::cout << "step function was nan detected: "
            << step_function_solver.wasNanDetected() << std::endl;

  SolverConfig sine_wave_config = {.kinematic_viscosity = 0.01,
                                   .time_steps = 5000,
                                   .domain_length = 2.0 * M_PI,
                                   .time_step_size = 0.001};

  // initial condition: one full sine wave over [0, 2π]
  std::function<double(double)> sine_function = [](double x) -> double {
    return std::sin(x);
  };

  BurgersSolver1d solver(std::make_unique<Godunov>(), sine_wave_config,
                         sine_function);

  solver.solve();
  solver.saveSolution("../data", "sine_wave", 1);
  std::cout << "sine wave was nan detected: "
            << step_function_solver.wasNanDetected() << std::endl;

  RandomInitialConditionConfig functionConfig;
  RandomInitialCondition f(functionConfig, false, true);
  std::cout << f.toString() << std::endl;

  SolverConfig random_function_config = {
      .kinematic_viscosity = 0.01,
      .time_steps = 10000,
      .domain_length = 10.0,
      .time_step_size = 0.0001,
  };

  BurgersSolver1d random_function_solver(std::make_unique<Godunov>(),
                                         random_function_config, f);

  random_function_solver.solve();
  random_function_solver.saveSolution("../data", "random_function", 1);
  std::cout << "random function was nan detected: "
            << step_function_solver.wasNanDetected() << std::endl;

  // 1D cases that match the 2D Godunov ICs exactly
  // These are meant to be run alongside the 2D solver to compare
  // the mid-y cross-section slice directly against the 1D solution.

  // Matches 2D square-wave: u=1 on [0.3, 0.7], domain [0, 1], inviscid
  SolverConfig match_squarewave_config = {.kinematic_viscosity = 0.0,
                                          .time_steps = 400,
                                          .domain_length = 1.0,
                                          .time_step_size = 0.001};

  std::function<double(double)> match_squarewave = [](double x) -> double {
    return (x >= 0.3 && x <= 0.7) ? 1.0 : 0.0;
  };

  BurgersSolver1d sq_solver(std::make_unique<Godunov>(),
                            match_squarewave_config, match_squarewave);

  sq_solver.solve();
  // gap=10 -> saves every 10th step -> 41 files, matching the 2D snapshot count
  sq_solver.saveSolution("../data", "match_squarewave_Godunov", 10);
  std::cout << "match_squarewave NaN detected: " << sq_solver.wasNanDetected()
            << std::endl;

  // Matches 2D Gaussian: exp(-(x-0.5)^2 / 0.02), sigma=0.1, domain [0, 1],
  // inviscid
  SolverConfig match_gaussian_config = {.kinematic_viscosity = 0.0,
                                        .time_steps = 400,
                                        .domain_length = 1.0,
                                        .time_step_size = 0.001};

  std::function<double(double)> match_gaussian = [](double x) -> double {
    return std::exp(-((x - 0.5) * (x - 0.5)) / 0.02);
  };

  BurgersSolver1d gauss_solver(std::make_unique<Godunov>(),
                               match_gaussian_config, match_gaussian);

  gauss_solver.solve();
  gauss_solver.saveSolution("../data", "match_gaussian_Godunov", 10);
  std::cout << "match_gaussian NaN detected: " << gauss_solver.wasNanDetected()
            << std::endl;
}
