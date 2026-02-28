#include "burgers_godunov2d.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

/*
    2D Unsplit Godunov Scheme for the inviscid Burgers system:

        u_t + (u^2/2)_x + (u*v)_y = 0
        v_t + (u*v)_x  + (v^2/2)_y = 0

    Conservative finite-volume, first-order accurate, periodic BCs.

    Each cell value lives at its centre: x_i = (i + 0.5)*dx, y_j = (j + 0.5)*dy.
    Data layout: flat vector, idx(i,j) = i*Ny + j (j is the fast index).
    All loops are  for(i) for(j)  so the inner loop walks through
    sequential memory addresses -- much better for cache performance.

    The update each step is:
        q^{n+1}_{i,j} = q^n_{i,j}
            - (dt/dx) * (F_right - F_left)
            - (dt/dy) * (G_top   - G_bottom)
*/

BurgersGodunovSolver2d::BurgersGodunovSolver2d(const Solver2dConfig &cfg)
    : Nx_(cfg.Nx), Ny_(cfg.Ny), Lx_(cfg.Lx), Ly_(cfg.Ly), dx_(cfg.Lx / cfg.Nx),
      dy_(cfg.Ly / cfg.Ny), CFL_(cfg.CFL), t_end_(cfg.t_end),
      save_interval_(cfg.save_interval), u_(cfg.Nx * cfg.Ny, 0.0),
      v_(cfg.Nx * cfg.Ny, 0.0), u_next_(cfg.Nx * cfg.Ny, 0.0),
      v_next_(cfg.Nx * cfg.Ny, 0.0) {}

// Fill the initial u and v fields by evaluating the user's IC lambdas
// at each cell centre.
void BurgersGodunovSolver2d::setInitialConditions(
    const std::function<double(double, double)> &u0,
    const std::function<double(double, double)> &v0) {
  ic_u_ = u0;
  ic_v_ = v0;

  // Cell centres: x = (i + 0.5)*dx,  y = (j + 0.5)*dy
  // i outer, j inner -> sequential memory access
  for (int i = 0; i < Nx_; ++i) {
    double x = (i + 0.5) * dx_;
    for (int j = 0; j < Ny_; ++j) {
      double y = (j + 0.5) * dy_;
      u_[idx(i, j)] = u0(x, y);
      v_[idx(i, j)] = v0(x, y);
    }
  }
}

// ---------------------------------------------------------------------------
// Riemann solvers
// ---------------------------------------------------------------------------

// Godunov flux for Burgers self-advection: f(q) = q^2/2, wave speed = q.
// This is the same logic used in the 1D solver -- just called from two
// directions in the 2D update.
double BurgersGodunovSolver2d::compute_flux(double qL, double qR) const {
  if (qL > qR) {
    // Shock. Rankine-Hugoniot gives shock speed s = (qL + qR) / 2.
    // Pick the upwind flux based on which way the shock is moving.
    double s = 0.5 * (qL + qR);
    return (s > 0.0) ? 0.5 * qL * qL : 0.5 * qR * qR;
  } else {
    // Rarefaction. Upwind based on sign of the wave speed.
    if (qL >= 0.0)
      return 0.5 * qL * qL;
    else if (qR <= 0.0)
      return 0.5 * qR * qR;
    else
      return 0.0; // sonic point, flux is zero
  }
}

// Upwind flux for a cross-advection term: flux = w * q.
// w is the wave speed (from the other velocity component),
// q is what's being carried. We Riemann-solve w to find which
// cell is the donor, then return w_donor * q_donor.
double BurgersGodunovSolver2d::compute_transport_flux(double wL, double wR,
                                                      double qL,
                                                      double qR) const {
  if (wL > wR) {
    // Shock in the w field.
    double s = 0.5 * (wL + wR);
    return (s > 0.0) ? wL * qL : wR * qR;
  } else {
    // Rarefaction in the w field.
    if (wL >= 0.0)
      return wL * qL;
    else if (wR <= 0.0)
      return wR * qR;
    else
      return 0.0;
  }
}

// ---------------------------------------------------------------------------
// Single time step
// ---------------------------------------------------------------------------
void BurgersGodunovSolver2d::step(double dt) {
  const double dtdx = dt / dx_;
  const double dtdy = dt / dy_;

  // --- x-direction fluxes at the right face of each cell (i+1/2, j) ---
  // Fu[i,j] is the flux of u across the right face of cell (i,j).
  // Fv[i,j] is the flux of v across the same face.
  std::vector<double> Fu(Nx_ * Ny_), Fv(Nx_ * Ny_);

  for (int i = 0; i < Nx_; ++i) {
    int ip1 = wrapI(i + 1);
    for (int j = 0; j < Ny_; ++j) {
      double uL = u_[idx(i, j)];
      double uR = u_[idx(ip1, j)];
      double vL = v_[idx(i, j)];
      double vR = v_[idx(ip1, j)];

      Fu[idx(i, j)] = compute_flux(uL, uR); // (u^2/2)_x for u-eqn
      Fv[idx(i, j)] =
          compute_transport_flux(uL, uR, vL, vR); // (u*v)_x for v-eqn
    }
  }

  // --- y-direction fluxes at the top face of each cell (i, j+1/2) ---
  std::vector<double> Gu(Nx_ * Ny_), Gv(Nx_ * Ny_);

  for (int i = 0; i < Nx_; ++i) {
    for (int j = 0; j < Ny_; ++j) {
      int jp1 = wrapJ(j + 1);

      double uB = u_[idx(i, j)];
      double uT = u_[idx(i, jp1)];
      double vB = v_[idx(i, j)];
      double vT = v_[idx(i, jp1)];

      Gu[idx(i, j)] =
          compute_transport_flux(vB, vT, uB, uT); // (u*v)_y for u-eqn
      Gv[idx(i, j)] = compute_flux(vB, vT);       // (v^2/2)_y for v-eqn
    }
  }

  // --- Conservative update (unsplit) ---
  // F_{i-1/2,j} is the right-face flux stored on the cell to the left: Fu[i-1,
  // j]. G_{i,j-1/2} is the top-face flux stored on the cell below:  Gu[i, j-1].
  for (int i = 0; i < Nx_; ++i) {
    int im1 = wrapI(i - 1);
    for (int j = 0; j < Ny_; ++j) {
      int jm1 = wrapJ(j - 1);

      u_next_[idx(i, j)] = u_[idx(i, j)] -
                           dtdx * (Fu[idx(i, j)] - Fu[idx(im1, j)]) -
                           dtdy * (Gu[idx(i, j)] - Gu[idx(i, jm1)]);

      v_next_[idx(i, j)] = v_[idx(i, j)] -
                           dtdx * (Fv[idx(i, j)] - Fv[idx(im1, j)]) -
                           dtdy * (Gv[idx(i, j)] - Gv[idx(i, jm1)]);
    }
  }

  // Swap buffers -- no copying, just pointer swap.
  std::swap(u_, u_next_);
  std::swap(v_, v_next_);
}

// ---------------------------------------------------------------------------
// Main time loop
// ---------------------------------------------------------------------------
void BurgersGodunovSolver2d::solve() {
  if (!ic_u_ || !ic_v_) {
    throw std::runtime_error("Call setInitialConditions() before solve().");
  }

  std::cout << "BurgersGodunovSolver2d: starting solve\n"
            << "  Grid: " << Nx_ << " x " << Ny_ << "  dx=" << dx_
            << "  dy=" << dy_ << "\n"
            << "  t_end=" << t_end_ << "  CFL=" << CFL_ << "\n";

  history_.clear();
  history_.push_back({0.0, u_, v_});

  double t = 0.0;
  double t_next_save = save_interval_;
  int step_count = 0;

  // Record the initial conservation totals so we can track drift.
  double sum_u0 = 0.0, sum_v0 = 0.0;
  for (int k = 0; k < Nx_ * Ny_; ++k) {
    sum_u0 += u_[k];
    sum_v0 += v_[k];
  }
  std::cout << std::fixed << std::setprecision(8);
  std::cout << "  Initial:  sum(u)=" << sum_u0 << "  sum(v)=" << sum_v0
            << "\n\n";

  while (t < t_end_) {
    // Find the fastest wave speed anywhere on the grid.
    double max_speed = 0.0;
    for (int k = 0; k < Nx_ * Ny_; ++k) {
      max_speed =
          std::max(max_speed, std::max(std::abs(u_[k]), std::abs(v_[k])));
    }

    if (max_speed < 1e-14) {
      std::cout << "  WARNING: flow has stalled, stopping early.\n";
      break;
    }

    // 2D CFL condition: dt = CFL * min(dx,dy) / max_speed
    double min_h = std::min(dx_, dy_);
    double dt = CFL_ * min_h / max_speed;

    // Trim dt so we land exactly on save times and t_end.
    dt = std::min(dt, t_end_ - t);
    dt = std::min(dt, t_next_save - t);

    // Sanity check -- this should never fire if CFL is set correctly.
    double cfl_x = dt * max_speed / dx_;
    double cfl_y = dt * max_speed / dy_;
    if (cfl_x > 1.0 || cfl_y > 1.0) {
      std::cerr << "  CFL VIOLATION at t=" << t << "  cfl_x=" << cfl_x
                << "  cfl_y=" << cfl_y << "\n";
    }

    step(dt);
    t += dt;
    ++step_count;

    // Print conservation check every 100 steps.
    // sum(u) and sum(v) should stay very close to their initial values.
    if (step_count % 100 == 0) {
      double sum_u = 0.0, sum_v = 0.0;
      for (int k = 0; k < Nx_ * Ny_; ++k) {
        sum_u += u_[k];
        sum_v += v_[k];
      }
      std::cout << "  t=" << std::setw(10) << t << "  step=" << std::setw(6)
                << step_count << "  sum(u)=" << sum_u
                << "  err_u=" << std::abs(sum_u - sum_u0)
                << "  max_speed=" << max_speed << "\n";
    }

    // Save snapshot if we've reached the next save time.
    if (t >= t_next_save - 1e-14) {
      history_.push_back({t, u_, v_});
      t_next_save += save_interval_;
    }
  }

  // Always capture the final state.
  if (history_.empty() || history_.back().t < t - 1e-14) {
    history_.push_back({t, u_, v_});
  }

  // Final report.
  double sum_u = 0.0, sum_v = 0.0;
  for (int k = 0; k < Nx_ * Ny_; ++k) {
    sum_u += u_[k];
    sum_v += v_[k];
  }
  std::cout << "\nDone.  steps=" << step_count << "  t_final=" << t << "\n"
            << "  Final sum(u)=" << sum_u
            << "  err=" << std::abs(sum_u - sum_u0) << "\n"
            << "  Final sum(v)=" << sum_v
            << "  err=" << std::abs(sum_v - sum_v0) << "\n"
            << "  Snapshots: " << history_.size() << "\n\n";
}

// ---------------------------------------------------------------------------
// Write all snapshots to disk
// ---------------------------------------------------------------------------
void BurgersGodunovSolver2d::saveSolution(const std::string &base_folder,
                                          const std::string &run_name) const {
  namespace fs = std::filesystem;

  if (!fs::exists(base_folder)) {
    fs::create_directories(base_folder);
  }

  std::string run_folder = base_folder + "/" + run_name;
  if (fs::exists(run_folder)) {
    std::cout << "Overwriting: " << run_folder << "\n";
    fs::remove_all(run_folder);
  }
  fs::create_directory(run_folder);

  std::cout << "Writing " << history_.size() << " snapshots to: " << run_folder
            << "\n";

  for (size_t snap = 0; snap < history_.size(); ++snap) {
    std::ostringstream fname;
    fname << run_folder << "/timestep_" << std::setw(5) << std::setfill('0')
          << snap << ".csv";

    std::ofstream file(fname.str());
    if (!file.is_open()) {
      std::cerr << "Failed to open: " << fname.str() << "\n";
      continue;
    }

    file << "x,y,u,v\n";
    const auto &su = history_[snap].u;
    const auto &sv = history_[snap].v;

    // i outer, j inner matches the flat layout, so this is a sequential read.
    for (int i = 0; i < Nx_; ++i) {
      double x = (i + 0.5) * dx_;
      for (int j = 0; j < Ny_; ++j) {
        double y = (j + 0.5) * dy_;
        int k = i * Ny_ + j;
        file << x << "," << y << "," << su[k] << "," << sv[k] << "\n";
      }
    }
  }

  std::cout << "Done.\n";
}
