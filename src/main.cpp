#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>

#include "burgers.h"

int main() {
    SolverConfig config = {
        .kinematic_viscosity = 0.01,
        .num_domain_points = 101,
        .time_steps = 200,
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

    BurgersSolver1d solver(config, step_function);

    solver.solve();

    solver.saveSolution("../data", "step_function");
}

