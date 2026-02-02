#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>

#include "burgers.h"
#include "burger_stencil.h"
#include "initial_condition_generator.h"

int main() {

    SolverConfig step_function_config = {
        .kinematic_viscosity = 0.01,
        .time_steps = 2000,
        .domain_length = 2.0,
        .time_step_size = 0.001
    };

    // 1 everywhere in domain
    // for 0.5 <= x <= 1 the function value is 2
    std::function<double(double)> step_function = [](double x) -> double {
        if (x >= 0.5 && x <= 1.0) {
            return 2.0;
        } else {
            return 1.0;
        }
    };

    BurgersSolver1d step_function_solver(
        std::make_unique<LaxWendroff>(),
        step_function_config, 
        step_function);

    step_function_solver.solve();
    step_function_solver.saveSolution("../data", "step_function", 1);
    std::cout << "step function was nan detected: " << step_function_solver.wasNanDetected() << std::endl;

    SolverConfig sine_wave_config = {
        .kinematic_viscosity = 0.01,
        .time_steps = 5000,
        .domain_length = 2.0 * M_PI,
        .time_step_size = 0.001
    };

    // initial condition: one full sine wave over [0, 2π]
    std::function<double(double)> sine_function = [](double x) -> double {
        return std::sin(x);
    };

    BurgersSolver1d solver(
        std::make_unique<LaxWendroff>(), 
        sine_wave_config, 
        sine_function);

    solver.solve();
    solver.saveSolution("../data", "sine_wave", 1);
    std::cout << "sine wave was nan detected: " << step_function_solver.wasNanDetected() << std::endl;

    
    RandomInitialConditionConfig functionConfig;
    RandomInitialCondition f(functionConfig, false, true);
    std::cout << f.toString() << std::endl;

    SolverConfig random_function_config = {
        .kinematic_viscosity = 0.01,
        .time_steps = 10000,
        .domain_length = 10.0,
        .time_step_size = 0.0001,
    };

    BurgersSolver1d random_function_solver(
        std::make_unique<LaxWendroff>(),
        random_function_config, 
        f);

    random_function_solver.solve();
    random_function_solver.saveSolution("../data", "random_function", 1);
    std::cout << "random function was nan detected: " << step_function_solver.wasNanDetected() << std::endl;

}
