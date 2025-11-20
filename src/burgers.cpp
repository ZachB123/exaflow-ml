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

void BurgersSolver1d::solve(double cq, Scheme scheme) {
    std::cout << "Solving...\n";

    solution_history.clear();
    std::vector<double> u_next(num_domain_points, 0.0);
    solution_history.push_back(u);

    for (int time_step = 0; time_step < time_steps; ++time_step) {
        switch (scheme) {
            case Scheme::FTCS:
            {
                for (int i = 1; i < num_domain_points - 1; ++i) {
            u_next[i] = u[i]
                - u[i] * time_step_size / spatial_step_size * (u[i + 1] - u[i - 1])
                + (kinematic_viscosity + calculateArtificialViscosity(cq, Scheme::FTCS, i)) * time_step_size / (spatial_step_size * spatial_step_size)
                  * (u[i + 1] - 2 * u[i] + u[i - 1]);
            }

            // wrap around
            u_next[0] = u[0]
                - u[0] * time_step_size / spatial_step_size * (u[0] - u[num_domain_points - 2])
                + kinematic_viscosity * time_step_size / (spatial_step_size * spatial_step_size)
                * (u[1] - 2 * u[0] + u[num_domain_points - 2]);

            u_next[num_domain_points - 1] = u_next[0];
            }
            break;
            case Scheme::LAX_WENDROFF:
            {
                double dt = time_step_size;
                double dx = spatial_step_size;

                // interior
                for (int i = 1; i < num_domain_points - 1; ++i) {

                    double f_ip = 0.5 * u[i+1] * u[i+1];
                    double f_i  = 0.5 * u[i]   * u[i];
                    double f_im = 0.5 * u[i-1] * u[i-1];

                    double a_ip = 0.5 * (u[i] + u[i+1]);
                    double a_im = 0.5 * (u[i-1] + u[i]);

                    double convective =
                        -(dt/(2*dx)) * (f_ip - f_im)
                        + (dt*dt/(2*dx*dx)) * ( a_ip * (f_ip - f_i) - a_im * (f_i - f_im) );

                    u_next[i] =
                        u[i] + convective +
                        (kinematic_viscosity + calculateArtificialViscosity(cq, Scheme::LAX_WENDROFF, i))
                        * dt/(dx*dx)
                        * (u[i+1] - 2*u[i] + u[i-1]);
                }

                // boundary i=0 (periodic)
                int i = 0;
                int ip = 1;
                int im = num_domain_points - 2;

                double f_ip = 0.5 * u[ip] * u[ip];
                double f_i  = 0.5 * u[i]  * u[i];
                double f_im = 0.5 * u[im] * u[im];

                double a_ip = 0.5 * (u[i] + u[ip]);
                double a_im = 0.5 * (u[im] + u[i]);

                double convective =
                    -(dt/(2*dx)) * (f_ip - f_im)
                    + (dt*dt/(2*dx*dx)) * ( a_ip * (f_ip - f_i) - a_im * (f_i - f_im) );

                u_next[i] =
                    u[i] + convective +
                    (kinematic_viscosity + calculateArtificialViscosity(cq, Scheme::LAX_WENDROFF, i))
                    * dt/(dx*dx)
                    * (u[ip] - 2*u[i] + u[im]);

                // enforce periodicity
                u_next[num_domain_points - 1] = u_next[0];
            }
            break;
            default:
                {
                    std::cerr << "Error: unknown scheme selected!\n";
                    return;   // or throw, depending on your design
                }
        }            

        std::swap(u, u_next);
        solution_history.push_back(u);
    }
}

double BurgersSolver1d::calculateArtificialViscosity(double cq, Scheme s, int i) {
    double artvis;
    double ux;

    switch (s) {
    case Scheme::FTCS:
        // linear artificial viscosity model, not quadratic RVN
        ux = (u[i + 1] - u[i - 1]) / (2.0 * spatial_step_size);
        // Only take artificial viscosity when in compression
        artvis = (ux < 0) ? cq * spatial_step_size * spatial_step_size * std::abs(ux) : 0.0;
        break;
    case Scheme::LAX_WENDROFF:
        ux = (u[i + 1] - u[i - 1]) / (2.0 * spatial_step_size);
        artvis = (ux < 0) ? cq * spatial_step_size * spatial_step_size * std::abs(ux) : 0.0;
        break;
    default:
        artvis = 0.0;
        break;
    }

    return artvis; 
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
        for (int i = 0; i < num_domain_points; ++i) {
            double x = i * spatial_step_size;  // compute x-coordinate
            file << x << "," << u_t[i] << "\n";
        }
    }

    std::cout << "All timesteps written successfully.\n";
}
