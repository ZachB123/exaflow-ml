#include <cmath>
#include <fstream>
#include <functional>
#include <vector>
#include <iostream>

#include "burger_scheme_2d.h"
#include "burgers_2d.h"

int main() {
    SolverConfig2D step_function_config = {.kinematic_viscosity = 0.01,
                                           .time_steps = 2000,
                                           .domain_length_x = 2.0,
                                           .domain_length_y = 2.0,
                                           .time_step_size = 0.001};

    // 1 everywhere in domain
    // for 0.5 <= x,y <= 1 the function value is 2
    std::function<double(double, double)> step_function = [](double x, double y) -> double {
        return (x >= 0.5 && x <= 1.0 && y >= 0.5 && y <= 1.0) ? 2.0 : 1.0;
    };

    BurgersSolver2d step_function_solver(std::make_unique<Godunov2D>(), step_function_config, step_function);

    step_function_solver.solve();
    step_function_solver.saveSolution("../data_2d", "step_function", 1);
    std::cout << "step function was nan detected: "
              << step_function_solver.wasNanDetected() << std::endl;
}