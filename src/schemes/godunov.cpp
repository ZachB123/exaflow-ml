#include "burger_scheme.h"
#include <vector>

/*
    Godunov scheme for 1D Burgers equation:

        u_t + (u^2/2)_x = nu u_xx

    - First-order finite volume
    - Exact Riemann solver for Burgers
    - Conservative flux-difference form
    - Periodic boundary conditions
    - No artificial viscosity
*/

void Godunov::calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double /*cq (unused)*/, int n /* size data array */, double dt, double dx, double nu /* kinematic_viscosity */)
{
    // Interfaces: between i and i+1 for i=0..n-3, plus periodic interface (n-2, 0)
    std::vector<double> F(n - 1);

    for (int i = 0; i < n - 2; ++i) {
        F[i] = godunovFlux(u[i], u[i + 1]);
    }
    F[n - 2] = godunovFlux(u[n - 2], u[0]);  // periodic interface

    // Update ALL physical cells: i = 0..n-2
    for (int i = 0; i <= n - 2; ++i) {
        int ip = (i == n - 2) ? 0     : i + 1;
        int im = (i == 0)     ? n - 2 : i - 1;

        // Right interface flux is F[i]
        double F_right = (i == n - 2) ? F[n - 2] : F[i];
        // Left interface flux is F[i-1], except for i=0 where it’s the periodic interface
        double F_left  = (i == 0)     ? F[n - 2] : F[i - 1];

        double convective = -(dt / dx) * (F_right - F_left);

        double diffusion  = nu * dt / (dx * dx) * (u[ip] - 2.0 * u[i] + u[im]);

        u_next[i] = u[i] + convective + diffusion;
    }

    // Enforce periodic duplicate point for convenience
    u_next[n - 1] = u_next[0];
}

// Exact Godunov flux for Burgers equation
// Solves Riemann problem between states uLeft and uRight, returns physical flux f(u*) at interface

double Godunov::godunovFlux(double uLeft, double uRight) const
{
    // Burgers flux: f(u) = 0.5 u^2

    if (uLeft > uRight) {
        // Shock case [Shock speed from Rankine-Hugoniot: s = (uLeft + uRight)/2]
        double s = 0.5 * (uLeft + uRight);

        if (s > 0.0)
            return 0.5 * uLeft * uLeft;
        else
            return 0.5 * uRight * uRight;
    }
    else {
        // Rarefaction case
        if (uLeft >= 0.0)
            return 0.5 * uLeft * uLeft;
        else if (uRight <= 0.0)
            return 0.5 * uRight * uRight;
        else
            return 0.0;  // fan contains u = 0
    }
}


std::string Godunov::getName() const {
    return "Godunov";
}
