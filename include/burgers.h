#pragma once
#include <functional>
#include <vector>

struct SolverConfig {
    // viscosity coefficient for viscous burgers
    double kinematic_viscosity;
    // how many points our domain is discretized to
    int num_domain_points;
    // how many iterations we will simulate
    int time_steps;
    // how long is the x axis are equation is defined on
    double domain_length;
    // length of time each time step is
    double time_step_size;
};

class BurgersSolver1d {
    const double kinematic_viscosity;

    const int num_domain_points;
    
    const int time_steps;

    [[maybe_unused]] const double domain_length;

    // how big of jumps we take in the domain duch that we have num_domain_points points
    const double spatial_step_size;

    const double time_step_size;

    // the current velocity field
    std::vector<double> u;

    // stores all u's as we solve
    std::vector<std::vector<double>> solution_history;

public:

    BurgersSolver1d(const SolverConfig& config);

    BurgersSolver1d(const SolverConfig& config, const std::function<double(double)>& initialize_conditions);

    void setInitialConditions(const std::function<double(double)>& initialize_conditions);

    void solve();

    std::vector<std::vector<double>> getSolution() const;

    void saveSolution(const std::string& base_folder, const std::string& run_name, int gap) const;
};