#include "burger_scheme_2d.h"
#include <cmath>

void LaxWendroff2D::calculateNextUandV(const std::vector<double>& u, const std::vector<double>& v, std::vector<double>& u_next, std::vector<double>& v_next, double cq, int num_domain_points_x, int num_domain_points_y, double time_step_size, double spatial_step_size_x, double spatial_step_size_y, double kinematic_viscosity) {
    double dt = time_step_size;
    double dx = spatial_step_size_x;
    double dy = spatial_step_size_y;

    auto wrap_x = [num_domain_points_x](int i) { return (i + num_domain_points_x - 1) % (num_domain_points_x - 1); };
    auto wrap_y = [num_domain_points_y](int j) { return (j + num_domain_points_y - 1) % (num_domain_points_y - 1); };

    for (int i = 0; i < num_domain_points_x - 1; ++i) {
        for (int j = 0; j < num_domain_points_y - 1; ++j) {

            double artvis = calculateArtificialViscosity(u, v, cq, dx, dy, i, j, num_domain_points_x, num_domain_points_y);
            int ip = wrap_x(i + 1);
            int im = wrap_x(i - 1);
            int jp = wrap_y(j + 1);
            int jm = wrap_y(j - 1);

            //calculating u_next[i][j]
            double f_ip = 0.5 * u[ip][j] * u[ip][j];
            double f_i  = 0.5 * u[i][j]  * u[i][j];
            double f_im = 0.5 * u[im][j] * u[im][j];

            double g_jp = u[i][jp] * v[i][jp];
            double g_i  = u[i][j]  * v[i][j];
            double g_jm = u[i][jm] * v[i][jm];

            double a_ip = 0.5 * (u[i][j] + u[ip][j]);
            double a_im = 0.5 * (u[im][j] + u[i][j]);

            double b_jp = 0.5 * (v[i][j] + v[i][jp]);
            double b_jm = 0.5 * (v[i][jm] + v[i][j]);

            double convective_u =
                - (dt / (2 * dx)) * (f_ip - f_im)
                + (dt * dt / (2 * dx * dx)) * (a_ip * (f_ip - f_i) - a_im * (f_i - f_im))
                - (dt / (2 * dy)) * (g_jp - g_jm)
                + (dt * dt / (2 * dy * dy)) * (b_jp * (g_jp - g_i) - b_jm * (g_i - g_jm));

            u_next[i][j] =
                u[i][j] + convective_u +
                (kinematic_viscosity + artvis) * dt *
                ((u[ip][j] - 2 * u[i][j] + u[im][j]) / (dx * dx)
               + (u[i][jp] - 2 * u[i][j] + u[i][jm]) / (dy * dy));

            //calculating v_next[i][j]
            double p_ip = u[ip][j] * v[ip][j];
            double p_i  = u[i][j]  * v[i][j]; 
            double p_im = u[im][j] * v[im][j]; 

            double q_jp = 0.5 * v[i][jp] * v[i][jp];
            double q_i  = 0.5 * v[i][j]  * v[i][j]; 
            double q_jm = 0.5 * v[i][jm] * v[i][jm]; 

            double c_ip = 0.5 * (u[i][j] + u[ip][j]); 
            double c_im = 0.5 * (u[im][j] + u[i][j]); 

            double d_jp = 0.5 * (v[i][j] + v[i][jp]); 
            double d_jm = 0.5 * (v[i][jm] + v[i][j]); 

            double convective_v =
                - (dt / (2 * dx)) * (p_ip - p_im)
                + (dt * dt / (2 * dx * dx)) * (c_ip * (p_ip - p_i) - c_im * (p_i - p_im))
                - (dt / (2 * dy)) * (q_jp - q_jm)
                + (dt * dt / (2 * dy * dy)) * (d_jp * (q_jp - q_i) - d_jm * (q_i - q_jm)); 

            v_next[i][j] =
                v[i][j] + convective_v +
                (kinematic_viscosity + artvis) * dt *
                ((v[ip][j] - 2 * v[i][j] + v[im][j]) / (dx * dx)
               + (v[i][jp] - 2 * v[i][j] + v[i][jm]) / (dy * dy)); 
        }
    }

    // enforce periodicity
    for (int j = 0; j < num_domain_points_y - 1; ++j) {
        u_next[num_domain_points_x - 1][j] = u_next[0][j];
        v_next[num_domain_points_x - 1][j] = u_next[0][j];
    }

    for (int i = 0; i < num_domain_points_x - 1; ++i) {
        u_next[i][num_domain_points_y - 1] = u_next[i][0];
        v_next[i][num_domain_points_y - 1] = u_next[i][0];
    }

    u_next[num_domain_points_x - 1][num_domain_points_y - 1] = u_next[0][0];
    v_next[num_domain_points_x - 1][num_domain_points_y - 1] = v_next[0][0];
}

double LaxWendroff2D::calculateArtificialViscosity(const std::vector<double>& u, const std::vector<double>& v, double cq, double spatial_step_size_x, double spatial_step_size_y, int i, int j, int num_domain_points_x, int num_domain_points_y) const {
    int ip = (i + 1) % (num_domain_points_x - 1);
    int im = (i - 1 + (num_domain_points_x - 1)) % (num_domain_points_x - 1);
    int jp = (j + 1) % (num_domain_points_y - 1);
    int jm = (j - 1 + (num_domain_points_y - 1)) % (num_domain_points_y - 1);
    
    double ux = (u[ip][j] - u[im][j]) / (2.0 * spatial_step_size_x);
    double vy = (v[i][jp] - v[i][jm]) / (2.0 * spatial_step_size_y);
    // Only take artificial viscosity when in compression
    double artvis = (ux + vy < 0) ? cq * 0.5 * (spatial_step_size_x * spatial_step_size_x + spatial_step_size_y * spatial_step_size_y) * std::abs(ux + vy) : 0.0;

    return artvis;
}

std::string LaxWendroff2D::getName() const {
    return "LaxWendroff2D";
}
