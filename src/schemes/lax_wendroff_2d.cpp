#include "burger_scheme_2d.h"
#include <cmath>

void LaxWendroff2D::calculateNextUandV(const std::vector<double>& u, const std::vector<double>& v, std::vector<double>& u_next, std::vector<double>& v_next, double cq, int num_domain_points_x, int num_domain_points_y, double time_step_size, double spatial_step_size_x, double spatial_step_size_y, double kinematic_viscosity) {
    double dt = time_step_size;
    double dx = spatial_step_size_x;
    double dy = spatial_step_size_y;

    auto wrap_x = [num_domain_points_x](int i) { return (i + num_domain_points_x - 1) % (num_domain_points_x - 1); };
    auto wrap_y = [num_domain_points_y](int j) { return (j + num_domain_points_y - 1) % (num_domain_points_y - 1); };

    auto idx = [num_domain_points_y](int i, int j) { return i * num_domain_points_y + j; };

    for (int i = 0; i < num_domain_points_x - 1; ++i) {
        for (int j = 0; j < num_domain_points_y - 1; ++j) {

            double artvis = calculateArtificialViscosity(u, v, cq, dx, dy, i, j, num_domain_points_x, num_domain_points_y);
            int ip = wrap_x(i + 1);
            int im = wrap_x(i - 1);
            int jp = wrap_y(j + 1);
            int jm = wrap_y(j - 1);

            //calculating u_next[i][j]
            double f_ip = 0.5 * u[idx(ip, j)] * u[idx(ip, j)];
            double f_i  = 0.5 * u[idx(i, j)]  * u[idx(i, j)];
            double f_im = 0.5 * u[idx(im, j)] * u[idx(im, j)];

            double g_jp = u[idx(i, jp)] * v[idx(i, jp)];
            double g_i  = u[idx(i, j)]  * v[idx(i, j)];
            double g_jm = u[idx(i, jm)] * v[idx(i, jm)];

            double a_ip = 0.5 * (u[idx(i, j)] + u[idx(ip, j)]);
            double a_im = 0.5 * (u[idx(im, j)] + u[idx(i, j)]);

            double b_jp = 0.5 * (v[idx(i, j)] + v[idx(i, jp)]);
            double b_jm = 0.5 * (v[idx(i, jm)] + v[idx(i, j)]);


            double convective_u =
                - (dt / (2 * dx)) * (f_ip - f_im)
                + (dt * dt / (2 * dx * dx)) * (a_ip * (f_ip - f_i) - a_im * (f_i - f_im))
                - (dt / (2 * dy)) * (g_jp - g_jm)
                + (dt * dt / (2 * dy * dy)) * (b_jp * (g_jp - g_i) - b_jm * (g_i - g_jm));

            u_next[idx(i, j)] =
                u[idx(i, j)] + convective_u +
                (kinematic_viscosity + artvis) * dt *
                ((u[idx(ip, j)] - 2 * u[idx(i, j)] + u[idx(im, j)]) / (dx * dx)
               + (u[idx(i, jp)] - 2 * u[idx(i, j)] + u[idx(i, jm)]) / (dy * dy));

            //calculating v_next[i][j]
            double p_ip = u[idx(ip, j)] * v[idx(ip, j)];
            double p_i  = u[idx(i, j)]  * v[idx(i, j)];
            double p_im = u[idx(im, j)] * v[idx(im, j)];

            double q_jp = 0.5 * v[idx(i, jp)] * v[idx(i, jp)];
            double q_i  = 0.5 * v[idx(i, j)]  * v[idx(i, j)];
            double q_jm = 0.5 * v[idx(i, jm)] * v[idx(i, jm)];

            double c_ip = 0.5 * (u[idx(i, j)] + u[idx(ip, j)]);
            double c_im = 0.5 * (u[idx(im, j)] + u[idx(i, j)]);

            double d_jp = 0.5 * (v[idx(i, j)] + v[idx(i, jp)]);
            double d_jm = 0.5 * (v[idx(i, jm)] + v[idx(i, j)]);

            double convective_v =
                - (dt / (2 * dx)) * (p_ip - p_im)
                + (dt * dt / (2 * dx * dx)) * (c_ip * (p_ip - p_i) - c_im * (p_i - p_im))
                - (dt / (2 * dy)) * (q_jp - q_jm)
                + (dt * dt / (2 * dy * dy)) * (d_jp * (q_jp - q_i) - d_jm * (q_i - q_jm)); 

            v_next[idx(i, j)] =
                v[idx(i, j)] + convective_v +
                (kinematic_viscosity + artvis) * dt *
                ((v[idx(ip, j)] - 2 * v[idx(i, j)] + v[idx(im, j)]) / (dx * dx)
               + (v[idx(i, jp)] - 2 * v[idx(i, j)] + v[idx(i, jm)]) / (dy * dy));
        }
    }

    // enforce periodicity
    for (int j = 0; j < num_domain_points_y - 1; ++j) {
        u_next[idx(num_domain_points_x - 1, j)] = u_next[idx(0, j)];
        v_next[idx(num_domain_points_x - 1, j)] = v_next[idx(0, j)];
    }

    for (int i = 0; i < num_domain_points_x - 1; ++i) {
        u_next[idx(i, num_domain_points_y - 1)] = u_next[idx(i, 0)];
        v_next[idx(i, num_domain_points_y - 1)] = v_next[idx(i, 0)];
    }

    u_next[idx(num_domain_points_x - 1, num_domain_points_y - 1)] = u_next[idx(0, 0)];
    v_next[idx(num_domain_points_x - 1, num_domain_points_y - 1)] = v_next[idx(0, 0)];
}

double LaxWendroff2D::calculateArtificialViscosity(const std::vector<double>& u, const std::vector<double>& v, double cq, double spatial_step_size_x, double spatial_step_size_y, int i, int j, int num_domain_points_x, int num_domain_points_y) const {
    int ip = (i + 1) % (num_domain_points_x - 1);
    int im = (i - 1 + (num_domain_points_x - 1)) % (num_domain_points_x - 1);
    int jp = (j + 1) % (num_domain_points_y - 1);
    int jm = (j - 1 + (num_domain_points_y - 1)) % (num_domain_points_y - 1);

    auto idx = [num_domain_points_y](int i, int j) { return i * num_domain_points_y + j; };

    double ux = (u[idx(ip, j)] - u[idx(im, j)]) / (2.0 * spatial_step_size_x);
    double vy = (v[idx(i, jp)] - v[idx(i, jm)]) / (2.0 * spatial_step_size_y);
    // Only take artificial viscosity when in compression
    double artvis = (ux + vy < 0) ? cq * 0.5 * (spatial_step_size_x * spatial_step_size_x + spatial_step_size_y * spatial_step_size_y) * std::abs(ux + vy) : 0.0;

    return artvis;
}

std::string LaxWendroff2D::getName() const {
    return "LaxWendroff2D";
}
