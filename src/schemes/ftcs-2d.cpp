#include <iostream>
#include <vector>
#include "burger_scheme_2d.h"
#include <cmath>

//'regular function instead for now, can redefine when FTCS2D is a class'
void calculateNextU(const std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& u_next, double cq, int nx, int ny, double dt, double dx, double dy, double nu){
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {

            double du_dx = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dx);
            double du_dy = (u[i][j + 1] - u[i][j - 1]) / (2.0 * dy);

            double d2u_dx2 = (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / (dx * dx);
            double d2u_dy2 = (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / (dy * dy);

            double artificial_viscosity = calculateArtificialViscosity(u, cq, dx, dy, i, j, nx, ny);

            u_next[i][j] = u[i][j] - u[i][j] * dt * (du_dx + du_dy) + (nu + artificial_viscosity) * dt * (d2u_dx2 + d2u_dy2);
        }
    }
}

double calculateArtificialViscosity(
    const std::vector<std::vector<double>>& u, double cq, double dx, double dy, int i, int j, int nx, int ny) {
    double ux = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dx);
    double uy = (u[i][j + 1] - u[i][j - 1]) / (2.0 * dy);

    double compression = ux + uy;
    double artificial_viscosity = (compression < 0) ? cq * (dx * dx + dy * dy) * std::abs(compression) : 0.0;

    return artificial_viscosity;    
}