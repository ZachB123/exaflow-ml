#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Defining constants

// Number of spatial points that the domain is discretized into
const int nx = 101;
// Number of time steps to simulate the evolution of the system
const int nt = 2000;
// Spatial step size, calculated based on the domain length (2.0) and number of points
const double dx = 2.0 / (nx - 1);
// Viscosity coefficient, representing the diffusion effect in the Burgers' equation
const double nu = 0.01;
// Time step size, determining how much time advances in each iteration
const double dt = 0.001;

int main()
{
    // Initialize the velocity field arrays
    // u: current time step, un: previous time step
    std::vector<double> u(nx, 0.0);
    std::vector<double> un(nx, 0.0);

    // Loop over spatial points to set initial conditions
    // Sets a step function as the initial velocity profile
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        if (x >= 0.5 && x <= 1.0)
            u[i] = 2.0;
        else
            u[i] = 1.0;
    }

    // Time-stepping loop
    // Updates the velocity field using the finite difference method
    for (int n = 0; n < nt; ++n)
    {
        un = u;
        for (int i = 1; i < nx - 1; ++i)
        {
            u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + nu * dt / (dx * dx) * (un[i + 1] - 2 * un[i] + un[i - 1]);
        }
    }

    // Writing results to a simple data file
    // The output file contains the spatial position and corresponding velocity
    std::ofstream fout("burgers_output.dat");
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        fout << x << " " << u[i] << "\n";
    }
    fout.close();

    // Indicates completion
    std::cout << "Simulation complete. Output written to burgers_output.dat\n";
    return 0;
}