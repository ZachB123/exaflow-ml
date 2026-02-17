#include "burger_stencil.h"
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

void Godunov::calculateNextU(const std::vector<double>& u, std::vector<double>& u_next, double /*cq (unused)*/, int N, double dt, double dx, double kinematic_viscosity) {
    // ---------------------------------------------------------
    // 1) Compute all interface fluxes F[i] = F_{i+1/2}
    // ---------------------------------------------------------
    // There are N-1 physical cells (last point duplicates first
    // for periodicity), so N-1 interfaces.
    // F[i] corresponds to interface between cell i and i+1.
    // ---------------------------------------------------------

    std::vector<double> F(N - 1);

    for (int i = 0; i < N - 2; ++i) {
        F[i] = godunovFlux(u[i], u[i + 1]);
    }

    // Periodic interface between last physical cell and first
    F[N - 2] = godunovFlux(u[N - 2], u[0]);

    // ---------------------------------------------------------
    // 2) Update interior cells using conservative flux difference
    //
    //     u_i^{n+1} =
    //         u_i^n
    //         - (dt/dx)(F_{i+1/2} - F_{i-1/2})
    //         + nu dt/dx^2 (u_{i+1} - 2u_i + u_{i-1})
    // ---------------------------------------------------------

    for (int i = 1; i < N - 1; ++i) {

        int ip = (i == N - 2) ? 0 : i + 1;  // periodic forward
        int im = (i == 0)     ? N - 2 : i - 1;

        double convective =
            - (dt / dx) * (F[i] - F[i - 1]);

        double diffusion =
            kinematic_viscosity * dt / (dx * dx)
            * (u[ip] - 2.0 * u[i] + u[im]);

        u_next[i] = u[i] + convective + diffusion;
    }

    // ---------------------------------------------------------
    // 3) Periodic boundary enforcement
    // ---------------------------------------------------------
    // Last entry duplicates first
    // ---------------------------------------------------------

    u_next[0] = u_next[N - 2];
    u_next[N - 1] = u_next[0];
}


// ------------------------------------------------------------------
// Exact Godunov flux for Burgers equation
//
// Solves Riemann problem between states uL and uR
// Returns physical flux f(u*) at interface
// ------------------------------------------------------------------

double Godunov::godunovFlux(double uL, double uR) const
{
    // Burgers flux: f(u) = 0.5 u^2

    if (uL > uR) {
        // ---------------------
        // Shock case
        // ---------------------
        // Shock speed from Rankine-Hugoniot:
        // s = (uL + uR)/2
        double s = 0.5 * (uL + uR);

        if (s > 0.0)
            return 0.5 * uL * uL;
        else
            return 0.5 * uR * uR;
    }
    else {
        // ---------------------
        // Rarefaction case
        // ---------------------
        if (uL >= 0.0)
            return 0.5 * uL * uL;
        else if (uR <= 0.0)
            return 0.5 * uR * uR;
        else
            return 0.0;  // fan contains u = 0
    }
}


std::string Godunov::getName() const {
    return "Godunov";
}
