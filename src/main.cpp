#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "burger_scheme.h"
#include "burgers.h"
#include "initial_condition_generator.h"

auto make_scheme(const std::string &name) -> std::unique_ptr<BurgerScheme> {
  if (name == "LaxWendroff")
    return std::make_unique<LaxWendroff>();
  return std::make_unique<Godunov>();
}

int main(int argc, char *argv[]) {
  std::string scheme = "Godunov";
  bool append = false;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-s" || arg == "--scheme") && i + 1 < argc)
      scheme = argv[++i];
    else if (arg == "-a" || arg == "--append")
      append = true;
  }

  if (!append) {
    for (auto &entry : std::filesystem::directory_iterator("../data"))
      if (entry.is_regular_file())
        std::filesystem::remove(entry.path());
  }
  // --- WIKIPEDIA TEST CASES ---
  std::cout << "\n--- WIKIPEDIA TEST CASES ---\n";

  SolverConfig wiki_config_base = {
      .kinematic_viscosity = 1.0, // Will be overwritten in loop
      .time_steps = 5000,
      .domain_length = 12.0, // Matches [-6, 6] range
      .time_step_size = 0.001};

  // Gaussian Initial Condition: u(x,0) = e^(-x^2/2)
  std::function<double(double)> gaussian_function = [](double x) -> double {
    double shifted_x = x - 6.0; // Center the plot at x=6 to mimic [-6, 6]
    return std::exp(-(shifted_x * shifted_x) / 2.0);
  };

  // N-wave type Initial Condition: u(x,0) = e^(-(x-1)^2/2) - e^(-(x+1)^2/2)
  std::function<double(double)> n_wave_function = [](double x) -> double {
    double shifted_x = x - 6.0; // Center the plot at x=6 to mimic [-6, 6]
    double part1 = std::exp(-(std::pow(shifted_x - 1.0, 2)) / 2.0);
    double part2 = std::exp(-(std::pow(shifted_x + 1.0, 2)) / 2.0);
    return part1 - part2;
  };

  std::vector<double> viscosities = {1.0, 0.1, 0.01};

  for (double v : viscosities) {
    wiki_config_base.kinematic_viscosity = v;
    std::string v_str = (v == 1.0) ? "1.0" : (v == 0.1) ? "0.1" : "0.01";

    std::cout << "Solving Gaussian for v=" << v << "...\n";
    BurgersSolver1d gaussian_solver(make_scheme(scheme), wiki_config_base,
                                    gaussian_function);
    gaussian_solver.solve();
    gaussian_solver.saveSolution("../data",
                                 "wiki_gaussian_" + v_str + "_" + scheme, 10);

    std::cout << "Solving N-wave for v=" << v << "...\n";
    BurgersSolver1d n_wave_solver(make_scheme(scheme), wiki_config_base,
                                  n_wave_function);
    n_wave_solver.solve();
    n_wave_solver.saveSolution("../data", "wiki_n_wave_" + v_str + "_" + scheme,
                               10);
  }

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
