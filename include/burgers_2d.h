#pragma once

#include <functional>
#include <vector>
#include <memory>

#include "burger_scheme_2d.h"

struct SolverConfig2D {
    double kinematic_viscosity;
    int time_steps;
    double domain_length_x;
    double domain_length_y;
    double time_step_size;
};

class BurgersSolver2d {
    private:
        // used for calculating spatial step size and number of domain points
        inline static const double ALPHA = 0.9;

        // scheme for solving
        std::unique_ptr<BurgerScheme2D> scheme;

        const double kinematic_viscosity;
        const int time_steps;
        const double domain_length_x;
        const double domain_length_y;
        const double time_step_size;

        // the current velocity fields
        std::vector<double> u;
        std::vector<double> v;

        // stores all u's as we solve
        struct Snapshot {
            std::vector<double> u;
            std::vector<double> v;
        };
        std::vector<Snapshot> solution_history;

        // the functions that we are using to define our initial conditions
        std::function<double(double, double)> initial_conditions_u;
        std::function<double(double, double)> initial_conditions_v;

        // how many points we calculate the value for in our domain
        int most_recent_num_domain_points_x;
        int most_recent_num_domain_points_y;

        // how big of jumps we take in the domain duch that we have num_domain_points points
        double most_recent_spatial_step_size_x;
        double most_recent_spatial_step_size_y;

        double approximate_max_u() const;
        double approximate_max_v() const;

        // tracks if a NaN was found during the solve
        bool nan_detected;

    public:
        BurgersSolver2d(
            const std::unique_ptr<BurgerScheme2D> scheme, 
            const SolverConfig2D& config);

        BurgersSolver2d(
            const std::unique_ptr<BurgerScheme2D> scheme, 
            const SolverConfig2D& config, 
            const std::function<double(double, double)> initial_conditions_u,
            const std::function<double(double, double)> initial_conditions_v);

        void setInitialConditions(const std::function<double(double, double)> &u_initial,
                                  const std::function<double(double, double)> &v_initial);

        void solve(double cq = 2.0);

        std::vector<Snapshot> getSolution() const;

        void saveSolution(const std::string& base_folder, const std::string& run_name, int gap) const;

        bool wasNanDetected() const;

        int getNumDomainPointsX() const;

        int getNumDomainPointsY() const;

        double getSpatialStepSizeX() const;

        double getSpatialStepSizeY() const;

        std::string getSchemeName() const;
};