#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <cmath>  

#include "burgers_2d.h"

BurgersSolver2d::BurgersSolver2d(std::unique_ptr<BurgerScheme2D> scheme, 
                                 const SolverConfig2D& config):
    scheme(std::move(scheme)),
    kinematic_viscosity(config.kinematic_viscosity),
    time_steps(config.time_steps),
    domain_length_x(config.domain_length_x),
    domain_length_y(config.domain_length_y),
    time_step_size(config.time_step_size),
    solution_history(),
    nan_detected(false)
{}

BurgersSolver2d::BurgersSolver2d(std::unique_ptr<BurgerScheme2D> scheme,
                                 const SolverConfig2D& config,
                                 const std::function<double(double, double)>& initialize_conditions):
    scheme(std::move(scheme)),
    kinematic_viscosity(config.kinematic_viscosity),
    time_steps(config.time_steps),
    domain_length_x(config.domain_length_x),
    domain_length_y(config.domain_length_y),
    time_step_size(config.time_step_size),
    solution_history(),
    initial_conditions(initial_conditions),
    nan_detected(false)
{}

void BurgersSolver2d::setInitialConditions(const std::function<double(double, double)>& initialize_conditions) {
    this->initial_conditions = initial_conditions;
}

double BurgersSolver2d::approximate_max_u() const {
    double max_u = 0.0;

    double approximate_step_size = time_step_size / ALPHA;
    int approximate_number_of_domain_points_x = std::floor(domain_length_x / approximate_step_size);
    // int approximate_number_of_domain_points_y = std::floor(domain_length_y / approximate_step_size);

    for (int i = 0; i < approximate_number_of_domain_points_x; i++) {
        double x = i * approximate_step_size;
        double y = i * approximate_step_size;
        max_u = std::max(max_u, std::abs(initial_conditions(x, y)));
    }

    return max_u;
}

