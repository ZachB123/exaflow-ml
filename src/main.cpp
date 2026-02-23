#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#include "burger_scheme.h"
#include "burgers.h"
#include "initial_condition_generator.h"

void run_all_schemes(const std::string &base_domain,
                     const std::string &base_name, const SolverConfig &config,
                     const std::function<double(double)> &func,
                     int save_gap = 1) {
  // 1. LaxWendroff Run
  std::cout << "\nRunning " << base_name << " on LaxWendroff..." << std::endl;
  BurgersSolver1d lax_solver(std::make_unique<LaxWendroff>(), config, func);
  lax_solver.solve();
  lax_solver.saveSolution(
      base_domain, base_name + "_" + lax_solver.getSchemeName(), save_gap);
  std::cout << base_name
            << " on LaxWendroff NaN detected: " << lax_solver.wasNanDetected()
            << std::endl;

  // 2. Godunov Run
  std::cout << "\nRunning " << base_name << " on Godunov..." << std::endl;
  BurgersSolver1d godunov_solver(std::make_unique<Godunov>(), config, func);
  godunov_solver.solve();
  godunov_solver.saveSolution(
      base_domain, base_name + "_" + godunov_solver.getSchemeName(), save_gap);
  std::cout << base_name
            << " on Godunov NaN detected: " << godunov_solver.wasNanDetected()
            << std::endl;
}

int main() {
  // --- PREVIOUS TEST CASES ---
  std::cout << "\n--- PREVIOUS TEST CASES ---\n";

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

  run_all_schemes("../data", "step_function", step_function_config,
                  step_function, 1);

  SolverConfig sine_wave_config = {.kinematic_viscosity = 0.01,
                                   .time_steps = 5000,
                                   .domain_length = 2.0 * M_PI,
                                   .time_step_size = 0.001};

  // initial condition: one full sine wave over [0, 2π]
  std::function<double(double)> sine_function = [](double x) -> double {
    return std::sin(x);
  };

  run_all_schemes("../data", "sine_wave", sine_wave_config, sine_function, 1);

  RandomInitialConditionConfig functionConfig;
  RandomInitialCondition f(functionConfig, false, true);
  std::cout << "\nGenerated Random Function String: " << f.toString()
            << std::endl;

  SolverConfig random_function_config = {
      .kinematic_viscosity = 0.01,
      .time_steps = 10000,
      .domain_length = 10.0,
      .time_step_size = 0.0001,
  };

  run_all_schemes("../data", "random_function", random_function_config, f, 1);
}
