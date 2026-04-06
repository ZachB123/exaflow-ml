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

namespace {
    inline int idx(int i, int j, int Ny) {
        return i * Ny + j;
    }
}

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

    if (!initial_conditions_u || !initial_conditions_v) {
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
    u.assign(num_domain_points_x * num_domain_points_y, 0.0);
    v.assign(num_domain_points_y * num_domain_points_y, 0.0);

    most_recent_spatial_step_size_x = spatial_step_size_x;
    most_recent_spatial_step_size_y = spatial_step_size_y;
    most_recent_num_domain_points_x = num_domain_points_x;
    most_recent_num_domain_points_y = num_domain_points_y;

    std::cout << "Computing with " << num_domain_points_x << "x" << num_domain_points_y << " domain points.\n";

    std::cout << "Setting initial conditions...\n";
    // i outer, j inner -> sequential memory access
    for (int i = 0; i < num_domain_points_x; ++i) {
        double x = i * spatial_step_size_x;
        for (int j = 0; j < num_domain_points_y; ++j) {
            double y = j * spatial_step_size_y;
            u[idx(i, j, num_domain_points_y)] = initial_conditions_u(x, y);
            v[idx(i, j, num_domain_points_y)] = initial_conditions_v(x, y);
        }
    }

    solution_history.clear();
    solution_history.push_back({u, v});
    std::vector<double> u_next(num_domain_points_x * num_domain_points_y, 0.0);
    std::vector<double> v_next(num_domain_points_x * num_domain_points_y, 0.0);

    for (int time_step = 0; time_step < time_steps; ++time_step) {
        scheme->calculateNextU(
            u,
            v,
            u_next,
            v_next,
            cq,
            num_domain_points_x,
            num_domain_points_y,
            time_step_size,
            spatial_step_size_x,
            spatial_step_size_y,
            kinematic_viscosity
        );

        if (!nan_detected) {
            for (int i = 0; i < num_domain_points_x; i++) {
                double x = i * spatial_step_size_x;
                for (int j = 0; j < num_domain_points_y; ++j) {
                    double y = j * spatial_step_size_y;
                    if (std::isnan(u_next[idx(i, j, num_domain_points_y)])) {
                        nan_detected = true;
                        break;
                    }
                }
            }
        }

        std::swap(u, u_next);
        std::swap(v, v_next);

        solution_history.push_back({u, v});
    }
}

std::vector<BurgersSolver2d::Snapshot> BurgersSolver2d::getSolution() const {
    return solution_history;
}

void BurgersSolver2d::saveSolution(const std::string& base_folder, const std::string& run_name, int gap) const {
    if (!std::filesystem::exists(base_folder)) {
        std::filesystem::create_directory(base_folder);
    }

    // create (or recreate) run folder
    std::string run_folder = base_folder + "/" + run_name;
    if (std::filesystem::exists(run_folder)) {
        std::cout << "Overwriting existing folder: " << run_folder << std::endl;
        std::filesystem::remove_all(run_folder); // delete everything inside
    }
    std::filesystem::create_directory(run_folder);

    // Write one file per time step
    std::cout << "Writing results to: " << run_folder << std::endl;
    for (size_t t = 0; t < solution_history.size(); ++t) {
        if (t % gap != 0) {
            continue;
        }

        std::ostringstream filename;
        filename << run_folder << "/timestep_" << std::setw(5) << std::setfill('0') << t << ".csv";
        
        std::ofstream file(filename.str());
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename.str() << std::endl;
            continue;
        }

        file << "x,y,u,v\n";
        const auto &s_u = solution_history[t].u;
        const auto &s_v = solution_history[t].v;
        for (int i = 0; i < most_recent_num_domain_points_x; i++) {
            double x = i * most_recent_spatial_step_size_x;
            for (int j = 0; j < most_recent_num_domain_points_y; ++j) {
                double y = j * most_recent_spatial_step_size_y;
                int k = idx(i, j, most_recent_num_domain_points_y);
                file << x << "," << y << "," << s_u[k] << "," << s_v[k] << "\n";
            }
        }   
    }

    std::cout << "All timesteps written successfully.\n";
}

bool BurgersSolver2d::wasNanDetected() const {
    return nan_detected;
}

int BurgersSolver2d::getNumDomainPointsX() const {
    return most_recent_num_domain_points_x;
}

double BurgersSolver2d::getSpatialStepSizeX() const {
    return most_recent_spatial_step_size_x;
}

int BurgersSolver2d::getNumDomainPointsY() const {
    return most_recent_num_domain_points_y;
}

double BurgersSolver2d::getSpatialStepSizeY() const {
    return most_recent_spatial_step_size_y;
}

std::string BurgersSolver2d::getSchemeName() const {
    return scheme->getName();
}
