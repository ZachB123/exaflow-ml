#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <cmath>  
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "burgers_2d.h"

/*
    Each cell value lives at its centre: x_i = (i + 0.5)*dx, y_j = (j + 0.5)*dy.
    Data layout: flat vector, idx(i,j) = i*Ny + j (j is the fast index).
    All loops are  for(i) for(j)  so the inner loop walks through sequential memory addresses
 */

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
                                 const std::function<double(double, double)> initial_conditions_u,
                                 const std::function<double(double, double)> initial_conditions_v):
    scheme(std::move(scheme)),
    kinematic_viscosity(config.kinematic_viscosity),
    time_steps(config.time_steps),
    domain_length_x(config.domain_length_x),
    domain_length_y(config.domain_length_y),
    time_step_size(config.time_step_size),
    solution_history(),
    initial_conditions_u(initial_conditions_u),
    initial_conditions_v(initial_conditions_v),
    nan_detected(false)
{}

void BurgersSolver2d::setInitialConditions(const std::function<double(double, double)> &u_initial,
                                           const std::function<double(double, double)> &v_initial) {
    this->initial_conditions_u = u_initial;
    this->initial_conditions_v = v_initial;
}

double BurgersSolver2d::approximate_max_u() const {
    double max_u = 0.0;

    double approximate_step_size = time_step_size / ALPHA;
    int approximate_number_of_domain_points_x = std::floor(domain_length_x / approximate_step_size);

    for (int i = 0; i < approximate_number_of_domain_points_x; i++) {
        double x = i * approximate_step_size;
        double y = i * approximate_step_size;
        max_u = std::max(max_u, std::abs(initial_conditions_u(x, y)));
    }

    return max_u;
}

double BurgersSolver2d::approximate_max_v() const {
    double max_v = 0.0;

    double approximate_step_size = time_step_size / ALPHA;
    int approximate_number_of_domain_points_y = std::floor(domain_length_y / approximate_step_size);

    for (int i = 0; i < approximate_number_of_domain_points_y; i++) {
        double x = i * approximate_step_size;
        double y = i * approximate_step_size;
        max_v = std::max(max_v, std::abs(initial_conditions_v(x, y)));
    }

    return max_v;
}

void BurgersSolver2d::solve(double cq) {
    std::cout << "Solving...\n";

    if (!initial_conditions_u or !initial_conditions_v) {
        throw std::runtime_error("No Initial Condition Function.");
    }

    nan_detected = false;
    double max_u = approximate_max_u();
    double max_v = approximate_max_v();
    
    // cfl condition
    double spatial_step_size_x = time_step_size * max_u / ALPHA;
    double spatial_step_size_y = time_step_size * max_v / ALPHA;
    // This spatial stepsize is bullshit and not actually what we want to use
    // if we use the previous calculation it just blows everything up
    spatial_step_size_x = 0.01;
    spatial_step_size_y = 0.01;
    double num_domain_points_x = std::floor(domain_length_x / spatial_step_size_x);
    double num_domain_points_y = std::floor(domain_length_y / spatial_step_size_y);
    u.assign(num_domain_points_x, 0.0);
    v.assign(num_domain_points_y, 0.0);

    most_recent_spatial_step_size_x = spatial_step_size_x;
    most_recent_spatial_step_size_y = spatial_step_size_y;
    most_recent_num_domain_points_x = num_domain_points_x;
    most_recent_num_domain_points_y = num_domain_points_y;

    std::cout << "Computing with " << num_domain_points_x << "/" << num_domain_points_y << " in x/y.\n";

    std::cout << "Setting initial conditions...\n";
    // Cell centers: x = (i + 0.5)*dx,  y = (j + 0.5)*dy
    // i outer, j inner -> sequential memory access
    for (int i = 0; i < num_domain_points_x; ++i) {
        double x = (i + 0.5) * spatial_step_size_x;
        for (int j = 0; j < num_domain_points_y; ++j) {
            double y = (j + 0.5) * spatial_step_size_y;
            u[i * num_domain_points_x + j] = initial_conditions_u(x, y);
            v[i * num_domain_points_y + j] = initial_conditions_v(x, y);
        }
    }
}

