#include <cmath>
#include <functional>
#include <iostream>

#include "burgers_godunov2d.h"

/*
    Driver for the 2D Godunov Burgers solver.

    Runs two test cases:

    1. Square-wave (top-hat): u = v = 1 inside the centre 40% of the domain,
       0 everywhere else. This produces sharp 2D shocks on all four edges of
       the patch that propagate diagonally outward.

    2. Gaussian bump: u = v = exp(-r^2 / (2*sigma^2)) centred at (Lx/2, Ly/2).
       Smooth at t=0, but the centre moves faster and eventually forms a
       circular shock front.

    Output goes to ../data_2d/<run_name>/timestep_XXXXX.csv
    Columns: x, y, u, v
*/

int main() {
  Solver2dConfig cfg;
  cfg.Nx = 200;
  cfg.Ny = 200;
  cfg.Lx = 1.0;
  cfg.Ly = 1.0;
  cfg.CFL = 0.45; // well below the 0.5 stability limit for 2D unsplit
  cfg.t_end = 0.4;
  cfg.save_interval = 0.01; // 40 snapshots + IC = 41 files total

  // ---- Test case 1: square-wave IC ----
  std::cout << "=== Test Case 1: Square-wave IC ===\n";
  {
    BurgersGodunovSolver2d solver(cfg);

    // u = v = 1 inside [0.3, 0.7] x [0.3, 0.7], 0 outside
    auto sq = [](double x, double y) -> double {
      return (x >= 0.3 && x <= 0.7 && y >= 0.3 && y <= 0.7) ? 1.0 : 0.0;
    };

    solver.setInitialConditions(sq, sq);
    solver.solve();
    solver.saveSolution("../data_2d", "godunov_squarewave");
  }

  // ---- Test case 2: Gaussian bump IC ----
  std::cout << "=== Test Case 2: Gaussian bump IC ===\n";
  {
    BurgersGodunovSolver2d solver(cfg);

    const double sigma = 0.1;
    const double cx = cfg.Lx / 2.0;
    const double cy = cfg.Ly / 2.0;

    // [=] captures sigma, cx, cy by value so the lambda stays valid
    auto gauss = [=](double x, double y) -> double {
      double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
      return std::exp(-r2 / (2.0 * sigma * sigma));
    };

    solver.setInitialConditions(gauss, gauss);
    solver.solve();
    solver.saveSolution("../data_2d", "godunov_gaussian");
  }

  std::cout << "All runs complete.\n";
  return 0;
}
