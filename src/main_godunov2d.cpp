#include <cmath>
#include <functional>
#include <iostream>

#include "burgers_godunov2d.h"

/*
    Driver for the 2D Unsplit Godunov Burgers solver.

    Runs two test cases back to back:
      1. Square-wave (top-hat) IC: u = v = 1 inside a centred patch, 0 outside.
         Produces a 2D shock visible as a sharp discontinuity propagating
         diagonally toward the domain corner.

      2. Gaussian bump IC: u = v = exp(-r^2 / (2*sigma^2)) centred at (Lx/2,
   Ly/2). Produces smooth steepening into a circular shock.

    Output data is written to ../data_2d/<run_name>/timestep_XXXXX.csv
    Each CSV: columns x, y, u, v.

    Usage:
        burgers_godunov2d.exe
*/

int main() {
  // ---- Shared grid / time setup ----
  Solver2dConfig cfg;
  cfg.Nx = 200; // cells in x
  cfg.Ny = 200; // cells in y
  cfg.Lx = 1.0; // domain [0, 1]
  cfg.Ly = 1.0;
  cfg.CFL = 0.45;           // safe 2D unsplit CFL (<0.5)
  cfg.t_end = 0.4;          // enough time to develop a clear shock
  cfg.save_interval = 0.01; // save every 0.01 time units -> 40 snapshots

  // ================================================================
  // Test Case 1: Square-wave (top-hat) initial condition
  // ================================================================
  std::cout << "========================================\n";
  std::cout << " Test Case 1: Square-wave IC\n";
  std::cout << "========================================\n";
  {
    BurgersGodunovSolver2d solver(cfg);

    // u = v = 1 inside [0.3, 0.7] x [0.3, 0.7], 0 elsewhere
    auto u0 = [&](double x, double y) -> double {
      if (x >= 0.3 && x <= 0.7 && y >= 0.3 && y <= 0.7)
        return 1.0;
      return 0.0;
    };
    auto v0 = u0; // same IC for v

    solver.setInitialConditions(u0, v0);
    solver.solve();
    solver.saveSolution("../data_2d", "godunov_squarewave");
  }

  // ================================================================
  // Test Case 2: Gaussian bump initial condition
  // ================================================================
  std::cout << "========================================\n";
  std::cout << " Test Case 2: Gaussian bump IC\n";
  std::cout << "========================================\n";
  {
    BurgersGodunovSolver2d solver(cfg);

    const double sigma = 0.1;
    const double cx = cfg.Lx / 2.0;
    const double cy = cfg.Ly / 2.0;

    auto u0 = [=](double x, double y) -> double {
      double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
      return std::exp(-r2 / (2.0 * sigma * sigma));
    };
    auto v0 = u0; // same IC for v

    solver.setInitialConditions(u0, v0);
    solver.solve();
    solver.saveSolution("../data_2d", "godunov_gaussian");
  }

  std::cout << "All runs complete.\n";
  return 0;
}
