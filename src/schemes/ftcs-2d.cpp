#include <iostream>
#include <vector>
#include "burger_scheme_2d.h"
#include <cmath>

void FTCS2D::calculateNextUandV(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& u_next,
    std::vector<std::vector<double>>& v_next,
    double cq, 
    int num_domain_points_x, 
    int num_domain_points_y, 
    double time_step_size, 
    double spatial_step_size_x, 
    double spatial_step_size_y, 
    double kinematic_viscosity)
    {
    for (int i = 1; i < num_domain_points_x - 1; ++i) {
        for (int j = 1; j < num_domain_points_y - 1; ++j) {

            double du_dx = (u[i + 1][j] - u[i - 1][j]) / (2.0 * spatial_step_size_x);
            double du_dy = (u[i][j + 1] - u[i][j - 1]) / (2.0 * spatial_step_size_y);

            double d2u_dx2 = (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / (spatial_step_size_x * spatial_step_size_x);
            double d2u_dy2 = (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / (spatial_step_size_y * spatial_step_size_y);

            double dv_dx = (v[i + 1][j] - v[i - 1][j]) / (2.0 * spatial_step_size_x);
            double dv_dy = (v[i][j + 1] - v[i][j - 1]) / (2.0 * spatial_step_size_y);

            double d2v_dx2 = (v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / (spatial_step_size_x * spatial_step_size_x);
            double d2v_dy2 = (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / (spatial_step_size_y * spatial_step_size_y);

            double artificial_viscosity = calculateArtificialViscosity(u, v, cq, spatial_step_size_x, spatial_step_size_y, i, j, num_domain_points_x, num_domain_points_y);

            u_next[i][j] = u[i][j] - time_step_size * (u[i][j] * du_dx + v[i][j] * du_dy) + (kinematic_viscosity + artificial_viscosity) * time_step_size * (d2u_dx2 + d2u_dy2);

            v_next[i][j] = v[i][j] - time_step_size * (u[i][j] * dv_dx + v[i][j] * dv_dy) + (kinematic_viscosity + artificial_viscosity) * time_step_size * (d2v_dx2 + d2v_dy2);
        }
    }
}


double FTCS2D::calculateArtificialViscosity(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v, 
    double cq, 
    double spatial_step_size_x, 
    double spatial_step_size_y, 
    int i, 
    int j, 
    int num_domain_points_x, 
    int num_domain_points_y) 
    const {

    double ux = (u[i + 1][j] - u[i - 1][j]) / (2.0 * spatial_step_size_x);
    double vy = (v[i][j + 1] - v[i][j - 1]) / (2.0 * spatial_step_size_y);

    double compression = ux + vy;

    double artificial_viscosity = (compression < 0) ? cq * (spatial_step_size_x * spatial_step_size_x + spatial_step_size_y * spatial_step_size_y) * std::abs(compression) : 0.0;

    return artificial_viscosity;    
}