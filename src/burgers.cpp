#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <sstream>  

#include "burgers.h"


BurgersSolver1d::BurgersSolver1d(const SolverConfig& config)
    :   kinematic_viscosity(config.kinematic_viscosity),
        num_domain_points(config.num_domain_points),
        time_steps(config.time_steps),
        domain_length(config.domain_length),
        spatial_step_size(config.domain_length / (config.num_domain_points - 1)),
        time_step_size(config.time_step_size),
        u(num_domain_points, 0.0),
        solution_history()
{}

BurgersSolver1d::BurgersSolver1d(
    const SolverConfig& config,
    const std::function<double(double)>& initialize_conditions
)
    :   kinematic_viscosity(config.kinematic_viscosity),
        num_domain_points(config.num_domain_points),
        time_steps(config.time_steps),
        domain_length(config.domain_length),
        spatial_step_size(config.domain_length / (config.num_domain_points - 1)),
        time_step_size(config.time_step_size),
        u(num_domain_points, 0.0),
        solution_history()
{
    setInitialConditions(initialize_conditions);
}

void BurgersSolver1d::setInitialConditions(const std::function<double(double)>& initialize_conditions) {
    std::cout << "Setting initial conditions...\n";
    for (int i = 0; i < num_domain_points; ++i) {
        double x = i * spatial_step_size;
        u[i] = initialize_conditions(x);
    }
}

void BurgersSolver1d::solve() {
    std::cout << "Solving...\n";

    solution_history.clear();
    std::vector<double> u_next(num_domain_points, 0.0);
    solution_history.push_back(u);

    for (int time_step = 0; time_step < time_steps; ++time_step) {

        for (int i = 1; i < num_domain_points - 1; ++i) {
            u_next[i] = u[i]
                - u[i] * time_step_size / spatial_step_size * (u[i] - u[i - 1])
                + kinematic_viscosity * time_step_size / (spatial_step_size * spatial_step_size)
                  * (u[i + 1] - 2 * u[i] + u[i - 1]);
        }

        // wrap around
        u_next[0] = u[0]
            - u[0] * time_step_size / spatial_step_size * (u[0] - u[num_domain_points - 2])
            + kinematic_viscosity * time_step_size / (spatial_step_size * spatial_step_size)
              * (u[1] - 2 * u[0] + u[num_domain_points - 2]);

        u_next[num_domain_points - 1] = u_next[0];

        std::swap(u, u_next);
        solution_history.push_back(u);
    }
}

std::vector<std::vector<double>> BurgersSolver1d::getSolution() const {
    return solution_history;
}


void BurgersSolver1d::saveSolution(const std::string& base_folder, const std::string& run_name) const {

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
        std::ostringstream filename;
        filename << run_folder << "/timestep_"
                 << std::setw(5) << std::setfill('0') << t << ".txt";

        std::ofstream file(filename.str());
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename.str() << std::endl;
            continue;
        }

        const auto& u_t = solution_history[t];
        for (int i = 0; i < num_domain_points; ++i) {
            double x = i * spatial_step_size;  // compute x-coordinate
            file << i << "," << x << "," << u_t[i] << "\n";
        }
    }

    std::cout << "All timesteps written successfully.\n";
}




