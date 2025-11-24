#include "burger_stencil.h"

void LaxWendroff::calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double cq, int num_domain_points, double time_step_size, double spatial_step_size, double kinematic_viscosity) {
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
            - (dt / (2 * dx)) * (f_ip - f_im)
            + (dt * dt / (2 * dx * dx)) * (a_ip * (f_ip - f_i) - a_im * (f_i - f_im));

        u_next[i] =
            u[i] + convective +
            (kinematic_viscosity + calculateArtificialViscosity(u, cq, spatial_step_size, i))
            * dt / (dx * dx)
            * (u[i + 1] - 2 * u[i] + u[i - 1]);
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
        - (dt / (2 * dx)) * (f_ip - f_im)
        + (dt * dt / (2 * dx * dx)) * ( a_ip * (f_ip - f_i) - a_im * (f_i - f_im));

    u_next[i] =
        u[i] + convective +
        (kinematic_viscosity + calculateArtificialViscosity(u, cq, spatial_step_size, i))
        * dt / (dx * dx)
        * (u[ip] - 2 * u[i] + u[im]);

    // enforce periodicity
    u_next[num_domain_points - 1] = u_next[0];
}

double LaxWendroff::calculateArtificialViscosity(const std::vector<double>& u, double cq, double spatial_step_size, int i) const {
    // linear artificial viscosity model, not quadratic RVN
    double ux = (u[i + 1] - u[i - 1]) / (2.0 * spatial_step_size);
    // Only take artificial viscosity when in compression
    double artvis = (ux < 0) ? cq * spatial_step_size * spatial_step_size * std::abs(ux) : 0.0;

    return artvis;
}

