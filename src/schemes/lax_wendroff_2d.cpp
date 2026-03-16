#include "burger_scheme_2d.h"
#include <cmath>

void LaxWendroff2D::calculateNextU(const std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& u_next, double cq, int num_domain_points_x, int num_domain_points_y, double time_step_size, double spatial_step_size_x, double spatial_step_size_y, double kinematic_viscosity) {
    double dt = time_step_size;
    double dx = spatial_step_size_x;
    double dy = spatial_step_size_y;

    // interior
    for (int i = 1; i < num_domain_points_x - 1; ++i) {
        for (int j = 1; j < num_domain_points_y - 1; ++j) {

            double f_ip = 0.5 * u[i + 1][j] * u[i + 1][j];
            double f_i  = 0.5 * u[i][j]     * u[i][j];
            double f_im = 0.5 * u[i - 1][j] * u[i - 1][j];

            double g_jp = 0.5 * u[i][j + 1] * u[i][j + 1];
            double g_j  = 0.5 * u[i][j]     * u[i][j];
            double g_jm = 0.5 * u[i][j - 1] * u[i][j - 1];

            double a_ip = 0.5 * (u[i][j] + u[i + 1][j]);
            double a_im = 0.5 * (u[i - 1][j] + u[i][j]);

            double b_jp = 0.5 * (u[i][j] + u[i][j + 1]);
            double b_jm = 0.5 * (u[i][j - 1] + u[i][j]);

            double convective =
                - (dt / (2 * dx)) * (f_ip - f_im)
                + (dt * dt / (2 * dx * dx)) * (a_ip * (f_ip - f_i) - a_im * (f_i - f_im))
                - (dt / (2 * dy)) * (g_jp - g_jm)
                + (dt * dt / (2 * dy * dy)) * (b_jp * (g_jp - g_j) - b_jm * (g_j - g_jm));

            u_next[i][j] =
                u[i][j] + convective +
                (kinematic_viscosity + calculateArtificialViscosity(u, cq, dx, dy, i, j, num_domain_points_x, num_domain_points_y))
                * dt
                * ((u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / (dx * dx)
                 + (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / (dy * dy));
        }
    }

    // boundary i=0 (periodic in x), interior j
    int i = 0;
    int ip = 1;
    int im = num_domain_points_x - 2;

    for (int j = 1; j < num_domain_points_y - 1; ++j) {
        double f_ip = 0.5 * u[ip][j] * u[ip][j];
        double f_i  = 0.5 * u[i][j]  * u[i][j];
        double f_im = 0.5 * u[im][j] * u[im][j];

        double g_jp = 0.5 * u[i][j + 1] * u[i][j + 1];
        double g_j  = 0.5 * u[i][j]     * u[i][j];
        double g_jm = 0.5 * u[i][j - 1] * u[i][j - 1];

        double a_ip = 0.5 * (u[i][j] + u[ip][j]);
        double a_im = 0.5 * (u[im][j] + u[i][j]);

        double b_jp = 0.5 * (u[i][j] + u[i][j + 1]);
        double b_jm = 0.5 * (u[i][j - 1] + u[i][j]);

        double convective =
            - (dt / (2 * dx)) * (f_ip - f_im)
            + (dt * dt / (2 * dx * dx)) * (a_ip * (f_ip - f_i) - a_im * (f_i - f_im))
            - (dt / (2 * dy)) * (g_jp - g_jm)
            + (dt * dt / (2 * dy * dy)) * (b_jp * (g_jp - g_j) - b_jm * (g_j - g_jm));

        u_next[i][j] =
            u[i][j] + convective +
            (kinematic_viscosity + calculateArtificialViscosity(u, cq, dx, dy, i, j, num_domain_points_x, num_domain_points_y))
            * dt
            * ((u[ip][j] - 2 * u[i][j] + u[im][j]) / (dx * dx)
             + (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / (dy * dy));
    }

    // boundary j=0 (periodic in y), interior i
    int j = 0;
    int jp = 1;
    int jm = num_domain_points_y - 2;

    for (int i = 1; i < num_domain_points_x - 1; ++i) {
        double f_ip = 0.5 * u[i + 1][j] * u[i + 1][j];
        double f_i  = 0.5 * u[i][j]     * u[i][j];
        double f_im = 0.5 * u[i - 1][j] * u[i - 1][j];

        double g_jp = 0.5 * u[i][jp] * u[i][jp];
        double g_j  = 0.5 * u[i][j]  * u[i][j];
        double g_jm = 0.5 * u[i][jm] * u[i][jm];

        double a_ip = 0.5 * (u[i][j] + u[i + 1][j]);
        double a_im = 0.5 * (u[i - 1][j] + u[i][j]);

        double b_jp = 0.5 * (u[i][j] + u[i][jp]);
        double b_jm = 0.5 * (u[i][jm] + u[i][j]);

        double convective =
            - (dt / (2 * dx)) * (f_ip - f_im)
            + (dt * dt / (2 * dx * dx)) * (a_ip * (f_ip - f_i) - a_im * (f_i - f_im))
            - (dt / (2 * dy)) * (g_jp - g_jm)
            + (dt * dt / (2 * dy * dy)) * (b_jp * (g_jp - g_j) - b_jm * (g_j - g_jm));

        u_next[i][j] =
            u[i][j] + convective +
            (kinematic_viscosity + calculateArtificialViscosity(u, cq, dx, dy, i, j, num_domain_points_x, num_domain_points_y))
            * dt
            * ((u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / (dx * dx)
             + (u[i][jp] - 2 * u[i][j] + u[i][jm]) / (dy * dy));
    }

    // corner i=0, j=0 (periodic in both x and y)
    i = 0;
    j = 0;
    ip = 1;
    im = num_domain_points_x - 2;
    jp = 1;
    jm = num_domain_points_y - 2;

    double f_ip = 0.5 * u[ip][j] * u[ip][j];
    double f_i  = 0.5 * u[i][j]  * u[i][j];
    double f_im = 0.5 * u[im][j] * u[im][j];

    double g_jp = 0.5 * u[i][jp] * u[i][jp];
    double g_j  = 0.5 * u[i][j]  * u[i][j];
    double g_jm = 0.5 * u[i][jm] * u[i][jm];

    double a_ip = 0.5 * (u[i][j] + u[ip][j]);
    double a_im = 0.5 * (u[im][j] + u[i][j]);

    double b_jp = 0.5 * (u[i][j] + u[i][jp]);
    double b_jm = 0.5 * (u[i][jm] + u[i][j]);

    double convective =
        - (dt / (2 * dx)) * (f_ip - f_im)
        + (dt * dt / (2 * dx * dx)) * (a_ip * (f_ip - f_i) - a_im * (f_i - f_im))
        - (dt / (2 * dy)) * (g_jp - g_jm)
        + (dt * dt / (2 * dy * dy)) * (b_jp * (g_jp - g_j) - b_jm * (g_j - g_jm));

    u_next[i][j] =
        u[i][j] + convective +
        (kinematic_viscosity + calculateArtificialViscosity(u, cq, dx, dy, i, j, num_domain_points_x, num_domain_points_y))
        * dt
        * ((u[ip][j] - 2 * u[i][j] + u[im][j]) / (dx * dx)
         + (u[i][jp] - 2 * u[i][j] + u[i][jm]) / (dy * dy));

    // enforce periodicity
    for (int j = 0; j < num_domain_points_y - 1; ++j) {
        u_next[num_domain_points_x - 1][j] = u_next[0][j];
    }

    for (int i = 0; i < num_domain_points_x - 1; ++i) {
        u_next[i][num_domain_points_y - 1] = u_next[i][0];
    }

    u_next[num_domain_points_x - 1][num_domain_points_y - 1] = u_next[0][0];
}

double LaxWendroff2D::calculateArtificialViscosity(const std::vector<std::vector<double>>& u, double cq, double spatial_step_size_x, double spatial_step_size_y, int i, int j, int num_domain_points_x, int num_domain_points_y) const {
    double ux = (i == 0) ? (u[i + 1][j] - u[num_domain_points_x - 2][j]) / (2.0 * spatial_step_size_x) : (u[i + 1][j] - u[i - 1][j]) / (2.0 * spatial_step_size_x);
    double uy = (j == 0) ? (u[i][j + 1] - u[i][num_domain_points_y - 2]) / (2.0 * spatial_step_size_y) : (u[i][j + 1] - u[i][j - 1]) / (2.0 * spatial_step_size_y);
    // Only take artificial viscosity when in compression
    double artvis = (ux + uy < 0) ? cq * 0.5 * (spatial_step_size_x * spatial_step_size_x + spatial_step_size_y * spatial_step_size_y) * std::abs(ux + uy) : 0.0;

    return artvis;
}

std::string LaxWendroff2D::getName() const {
    return "LaxWendroff2D";
}