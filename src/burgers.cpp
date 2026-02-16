#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <cmath>  

#include "burgers.h"


BurgersSolver1d::BurgersSolver1d(std::unique_ptr<BurgerStencil> stencil, const SolverConfig& config)
    :   
        stencil(std::move(stencil)),
        kinematic_viscosity(config.kinematic_viscosity),
        time_steps(config.time_steps),
        domain_length(config.domain_length),
        time_step_size(config.time_step_size),
        solution_history(),
        nan_detected(false)
{}

BurgersSolver1d::BurgersSolver1d(
    std::unique_ptr<BurgerStencil> stencil,
    const SolverConfig& config,
    const std::function<double(double)>& initial_conditions
)
    :   
        stencil(std::move(stencil)),
        kinematic_viscosity(config.kinematic_viscosity),
        time_steps(config.time_steps),
        domain_length(config.domain_length),
        time_step_size(config.time_step_size),
        solution_history(),
        initial_conditions(initial_conditions),
        nan_detected(false)
{}

void BurgersSolver1d::setInitialConditions(const std::function<double(double)>& initial_conditions) {
    this->initial_conditions = initial_conditions;
}

double BurgersSolver1d::approximate_max_u() const {
    double max_u = 0.0;
    // this is a good approximation since the actual step size will be (time_step_size * max_u) / ALPHA
    double approximate_step_size = time_step_size / ALPHA;
    int approximate_number_of_domain_points = std::floor(domain_length / approximate_step_size);

    for (int i = 0; i < approximate_number_of_domain_points; ++i) {
        double x = i * approximate_step_size;
        max_u = std::max(max_u, std::abs(initial_conditions(x)));
    }

    return max_u;
}

void BurgersSolver1d::solve(double cq) {
    std::cout << "Solving...\n";

    if (!initial_conditions) {
        throw std::runtime_error("No Initial Condition Function.");
    }
    
    nan_detected = false;
    double max_u = approximate_max_u();
    // cfl condition
    double spatial_step_size = time_step_size * max_u * ALPHA;
    // This spatial stepsize is bullshit and not actually what we want to use
    // if we use the previous calculation it just blows everything up
    spatial_step_size = 0.01;
    // spatial_step_size = domain_length / 1000.0;
    double num_domain_points = std::floor(domain_length / spatial_step_size);
    u.assign(num_domain_points, 0.0);

    most_recent_spatial_step_size = spatial_step_size;
    most_recent_num_domain_points = num_domain_points;

    std::cout << "Computing with " << num_domain_points << " domain points.\n";

    std::cout << "Setting initial conditions...\n";
    for (int i = 0; i < num_domain_points; ++i) {
        double x = i * spatial_step_size;
        u[i] = initial_conditions(x);
    }

    solution_history.clear();
    std::vector<double> u_next(num_domain_points, 0.0);
    solution_history.push_back(u);

    for (int time_step = 0; time_step < time_steps; ++time_step) {
        stencil->calculateNextU(
            u, 
            u_next, 
            cq, 
            num_domain_points, 
            time_step_size, 
            spatial_step_size, 
            kinematic_viscosity
        );       
        
        // figure out if the solver crashed anywhere
        if (!nan_detected) {
            for (int i = 0; i < num_domain_points; i++) {
                if (std::isnan(u_next[i])) {
                    nan_detected = true;
                    break;
                }
            }
        }

        std::swap(u, u_next);
        solution_history.push_back(u);
    }
}

std::vector<std::vector<double>> BurgersSolver1d::getSolution() const {
    return solution_history;
}

void BurgersSolver1d::saveSolution(const std::string& base_folder, const std::string& run_name, int gap) const {

    if (!std::filesystem::exists(base_folder)) {
        std::filesystem::create_directory(base_folder);
    }

    // Create (or recreate) run folder
    std::string run_folder = base_folder + "/" + run_name;
    if (std::filesystem::exists(run_folder)) {
        std::cout << "Overwriting existing folder: " << run_folder << std::endl;
        std::filesystem::remove_all(run_folder); // delete everything inside
    }
    std::filesystem::create_directory(run_folder);

    std::cout << "Writing results to: " << run_folder << std::endl;

    // Write one file per time step
    for (size_t t = 0; t < solution_history.size(); ++t) {

        if (t % gap != 0) {
            continue;
        }

        std::ostringstream filename;
        filename << run_folder << "/timestep_"
                 << std::setw(5) << std::setfill('0') << t << ".csv";

        std::ofstream file(filename.str());
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename.str() << std::endl;
            continue;
        }

        file << "x,u\n";

        const auto& u_t = solution_history[t];
        for (int i = 0; i < most_recent_num_domain_points; ++i) {
            double x = i * most_recent_spatial_step_size;  // compute x-coordinate
            file << x << "," << u_t[i] << "\n";
        }
    }

    std::cout << "All timesteps written successfully.\n";
}

bool BurgersSolver1d::wasNanDetected() const {
    return nan_detected;
}

int BurgersSolver1d::getNumDomainPoints() const {
    return most_recent_num_domain_points;
}

double BurgersSolver1d::getSpatialStepSize() const {
    return most_recent_spatial_step_size;
}

std::string BurgersSolver1d::getStencilName() const {
    return stencil->getName();
}