#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Defining constants

// Number of spatial points that the domain is discretized into
// This number could change based on the number of times we want to calculate the velocity
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
    // Updates the velocity for each time step using the burger's equation
    for (int n = 0; n < nt; ++n)
    {
        un = u;
        for (int i = 1; i < nx - 1; ++i)
        {
            u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + nu * dt / (dx * dx) * (un[i + 1] - 2 * un[i] + un[i - 1]);
        }
    }

    // Writing results to a simple data file
    // The output file contains the spatial position (x) and corresponding velocity u(x, t)
    // For instance, if the output file has a row that reads "1.0 1.7", that means that the position is 1.0 and the velocity is 1.7
    std::ofstream fout("burgers_output_cpp.dat");
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        fout << x << " " << u[i] << "\n";
    }
    fout.close();

    // Prints that the process has been completed
    std::cout << "Simulation complete. Output written to burgers_output.dat\n";
    return 0;
}
