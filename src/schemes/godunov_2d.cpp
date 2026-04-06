#include "burger_scheme_2d.h"
#include <algorithm>
#include <cmath>
#include <vector>

/*
    Godunov scheme for 2D inviscid Burgers system:

        u_t + (u^2/2)_x + (u v)_y = 0
        v_t + (u v)_x   + (v^2/2)_y = 0

    - First-order finite volume
    - Conservative flux-difference form
    - Periodic boundary conditions
    - No artificial viscosity in this Godunov implementation
*/

namespace {
    inline int idx(int i, int j, int Ny) {
        return i * Ny + j;
    }

    inline int wrap(int a, int n) {
        return (a % n + n) % n;
    }
}

void Godunov2D::calculateNextUandV(
    const std::vector<double>& u,
    const std::vector<double>& v,
    std::vector<double>& u_next,
    std::vector<double>& v_next,
    double /* cq unused */,
    int Nx,
    int Ny,
    double dt,
    double dx,
    double dy,
    double /* kinematic_viscosity unused */
) {
    const double dtdx = dt / dx;
    const double dtdy = dt / dy;

    // Fluxes on x-faces: interface between (i,j) and (i+1,j)
    std::vector<double> Fu(Nx * Ny, 0.0);
    std::vector<double> Fv(Nx * Ny, 0.0);

    for (int i = 0; i < Nx; ++i) {
        int ip = wrap(i + 1, Nx);
        for (int j = 0; j < Ny; ++j) {
            double uL = u[idx(i, j, Ny)];
            double uR = u[idx(ip, j, Ny)];
            double vL = v[idx(i, j, Ny)];
            double vR = v[idx(ip, j, Ny)];

            Fu[idx(i, j, Ny)] = godunovFlux(uL, uR);
            Fv[idx(i, j, Ny)] = transportFlux(uL, uR, vL, vR);
        }
    }

    // Fluxes on y-faces: interface between (i,j) and (i,j+1)
    std::vector<double> Gu(Nx * Ny, 0.0);
    std::vector<double> Gv(Nx * Ny, 0.0);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int jp = wrap(j + 1, Ny);

            double uB = u[idx(i, j, Ny)];
            double uT = u[idx(i, jp, Ny)];
            double vB = v[idx(i, j, Ny)];
            double vT = v[idx(i, jp, Ny)];

            Gu[idx(i, j, Ny)] = transportFlux(vB, vT, uB, uT);
            Gv[idx(i, j, Ny)] = godunovFlux(vB, vT);
        }
    }

    // Conservative unsplit update
    for (int i = 0; i < Nx; ++i) {
        int im = wrap(i - 1, Nx);
        for (int j = 0; j < Ny; ++j) {
            int jm = wrap(j - 1, Ny);

            u_next[idx(i, j, Ny)] =
                u[idx(i, j, Ny)]
                - dtdx * (Fu[idx(i, j, Ny)] - Fu[idx(im, j, Ny)])
                - dtdy * (Gu[idx(i, j, Ny)] - Gu[idx(i, jm, Ny)]);

            v_next[idx(i, j, Ny)] =
                v[idx(i, j, Ny)]
                - dtdx * (Fv[idx(i, j, Ny)] - Fv[idx(im, j, Ny)])
                - dtdy * (Gv[idx(i, j, Ny)] - Gv[idx(i, jm, Ny)]);
        }
    }
}

// Exact Godunov flux for scalar Burgers flux f(q)=0.5*q^2
double Godunov2D::godunovFlux(double qLeft, double qRight) const {
    if (qLeft > qRight) {
        // shock
        double s = 0.5 * (qLeft + qRight);
        if (s > 0.0) {
            return 0.5 * qLeft * qLeft;
        } else {
            return 0.5 * qRight * qRight;
        }
    } else {
        // rarefaction
        if (qLeft >= 0.0) {
            return 0.5 * qLeft * qLeft;
        } else if (qRight <= 0.0) {
            return 0.5 * qRight * qRight;
        } else {
            return 0.0;
        }
    }
}

// Upwind transport flux for flux = w * q
double Godunov2D::transportFlux(double wLeft, double wRight,
                                double qLeft, double qRight) const {
    if (wLeft > wRight) {
        // shock in transport speed
        double s = 0.5 * (wLeft + wRight);
        return (s > 0.0) ? wLeft * qLeft : wRight * qRight;
    } else {
        // rarefaction in transport speed
        if (wLeft >= 0.0) {
            return wLeft * qLeft;
        } else if (wRight <= 0.0) {
            return wRight * qRight;
        } else {
            return 0.0;
        }
    }
}

std::string Godunov2D::getName() const {
    return "Godunov2D";
}