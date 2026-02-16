#pragma once
#include <functional>
#include <vector>

#include "burger_stencil.h"

struct SolverConfig {
    // viscosity coefficient for viscous burgers
    double kinematic_viscosity;
    // how many iterations we will simulate
    int time_steps;
    // how long is the x axis are equation is defined on
    double domain_length;
    // length of time each time step is
    double time_step_size;
};

class BurgersSolver1d {

private:
    // used for calculating spatial step size and number of domain points
    inline static const double ALPHA = 0.9;

    // stencil for solving
    std::unique_ptr<BurgerStencil> stencil;

    // kinematic viscosity as shown in the burgers equation
    const double kinematic_viscosity;
    
    // how many times will we run our solver
    const int time_steps;

    // how long is our initial condition domain (x axis)
    const double domain_length;

    // how much time do we step for each solver iteration
    const double time_step_size;

    // the current velocity field
    std::vector<double> u;

    // stores all u's as we solve
    std::vector<std::vector<double>> solution_history;

    // the function that we are using to define our initial conditions
    std::function<double(double)> initial_conditions;

    // how many points we calculate the value for in our domain
    int most_recent_num_domain_points;

    // how big of jumps we take in the domain duch that we have num_domain_points points
    double most_recent_spatial_step_size;

    double approximate_max_u() const;

    // tracks if a NaN was found during the solve
    bool nan_detected;

public:
    
    BurgersSolver1d(
        const std::unique_ptr<BurgerStencil> stencil, 
        const SolverConfig& config);

    BurgersSolver1d(
        const std::unique_ptr<BurgerStencil> stencil, 
        const SolverConfig& config, 
        const std::function<double(double)>& initialize_conditions);

    void setInitialConditions(const std::function<double(double)>& initialize_conditions);

    void solve(double cq = 2.0);

    std::vector<std::vector<double>> getSolution() const;

    void saveSolution(const std::string& base_folder, const std::string& run_name, int gap) const;

    bool wasNanDetected() const;

    int getNumDomainPoints() const;

    double getSpatialStepSize() const;

    std::string getStencilName() const;
};
