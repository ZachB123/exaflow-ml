#include <iostream>

#include "burger_scheme.h"

void FTCS::calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double cq, int num_domain_points, double time_step_size, double spatial_step_size, double kinematic_viscosity) {
    for (int i = 1; i < num_domain_points - 1; ++i) {
        u_next[i] = u[i]
            - u[i] * time_step_size / spatial_step_size * (u[i + 1] - u[i - 1])
            + (kinematic_viscosity + calculateArtificialViscosity(u, cq, spatial_step_size, i, num_domain_points)) * time_step_size / (spatial_step_size * spatial_step_size)
                * (u[i + 1] - 2 * u[i] + u[i - 1]);
    }

    // wrap around
    u_next[0] = u[0]
        - u[0] * time_step_size / spatial_step_size * (u[0] - u[num_domain_points - 2])
        + kinematic_viscosity * time_step_size / (spatial_step_size * spatial_step_size)
        * (u[1] - 2 * u[0] + u[num_domain_points - 2]);

    u_next[num_domain_points - 1] = u_next[0];
}

double FTCS::calculateArtificialViscosity(const std::vector<double>& u, double cq, double spatial_step_size, int i, int num_domain_points) const {
    // linear artificial viscosity model, not quadratic RVN
    double ux = (i == 0) ? (u[i + 1] - u[num_domain_points - 2]) / (2.0 * spatial_step_size) : (u[i + 1] - u[i - 1]) / (2.0 * spatial_step_size);
    // Only take artificial viscosity when in compression
    double artificial_viscosity = (ux < 0) ? cq * spatial_step_size * spatial_step_size * std::abs(ux) : 0.0;

    return artificial_viscosity;
}

std::string FTCS::getName() const {
    return "FTCS";
}