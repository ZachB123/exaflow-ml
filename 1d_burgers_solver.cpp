#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Define constants
const int nx = 101;  // Number of spatial points
const int nt = 2000; // Number of time steps
const double dx = 2.0 / (nx - 1);
const double nu = 0.01;  // Viscosity
const double dt = 0.001; // Time step

int main()
{
    std::vector<double> u(nx, 0.0);  // Current time step
    std::vector<double> un(nx, 0.0); // Previous time step

    // Initial condition: u = 2 between 0.5 and 1.0
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        if (x >= 0.5 && x <= 1.0)
            u[i] = 2.0;
        else
            u[i] = 1.0;
    }

    // Time-stepping loop
    for (int n = 0; n < nt; ++n)
    {
        un = u;
        for (int i = 1; i < nx - 1; ++i)
        {
            u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1])                // Convective term (upwind)
                   + nu * dt / (dx * dx) * (un[i + 1] - 2 * un[i] + un[i - 1]); // Diffusion
        }
    }

    // Write results to a file
    std::ofstream fout("burgers_output.dat");
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        fout << x << " " << u[i] << "\n";
    }
    fout.close();

    std::cout << "Simulation complete. Output written to burgers_output.dat\n";
    return 0;
}