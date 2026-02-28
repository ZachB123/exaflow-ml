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

    Finite-volume, first-order, periodic BCs.

    Data layout: flat vector, idx(i,j) = i*Ny + j  (row-major, j fast).
    All nested loops are  for(i) for(j)  to keep j as the fast-varying
    index and maximise cache-line utilisation.

    Update formula (unsplit):
        q^{n+1}_{i,j} = q^n_{i,j}
            - (dt/dx) * (F_{i+1/2,j} - F_{i-1/2,j})   <- x fluxes
            - (dt/dy) * (G_{i,j+1/2} - G_{i,j-1/2})   <- y fluxes
*/

// ============================================================
//  Constructor
// ============================================================
BurgersGodunovSolver2d::BurgersGodunovSolver2d(const Solver2dConfig &cfg)
    : Nx_(cfg.Nx), Ny_(cfg.Ny), Lx_(cfg.Lx), Ly_(cfg.Ly), dx_(cfg.Lx / cfg.Nx),
      dy_(cfg.Ly / cfg.Ny), CFL_(cfg.CFL), t_end_(cfg.t_end),
      save_interval_(cfg.save_interval), u_(cfg.Nx * cfg.Ny, 0.0),
      v_(cfg.Nx * cfg.Ny, 0.0), u_next_(cfg.Nx * cfg.Ny, 0.0),
      v_next_(cfg.Nx * cfg.Ny, 0.0) {}

// ============================================================
//  Initial conditions
// ============================================================
void BurgersGodunovSolver2d::setInitialConditions(
    const std::function<double(double, double)> &u0,
    const std::function<double(double, double)> &v0) {
  ic_u_ = u0;
  ic_v_ = v0;

  // Cell centres: x_i = (i + 0.5)*dx,  y_j = (j + 0.5)*dy
  // Loop order: i outer, j inner  -> sequential memory access
  for (int i = 0; i < Nx_; ++i) {
    double x = (i + 0.5) * dx_;
    for (int j = 0; j < Ny_; ++j) {
      double y = (j + 0.5) * dy_;
      u_[idx(i, j)] = u0(x, y);
      v_[idx(i, j)] = v0(x, y);
    }
  }
}

// ============================================================
//  Riemann solvers
// ============================================================

// Godunov flux for self-advection f(q) = q^2/2, wave speed = q.
// Exact entropy-satisfying solution to the Burgers Riemann problem.
double BurgersGodunovSolver2d::compute_flux(double qL, double qR) const {
  if (qL > qR) {
    // Shock: Rankine-Hugoniot speed s = (qL + qR)/2
    double s = 0.5 * (qL + qR);
    return (s > 0.0) ? 0.5 * qL * qL : 0.5 * qR * qR;
  } else {
    // Rarefaction (entropy fan)
    if (qL >= 0.0)
      return 0.5 * qL * qL;
    else if (qR <= 0.0)
      return 0.5 * qR * qR;
    else
      return 0.0; // sonic point: u* = 0
  }
}

// Upwind donor-cell flux for cross-term: flux = w * q
// Wave speed is w (from the v or u field); transported quantity is q.
// Uses the same shock/rarefaction logic on the w Riemann problem
// to decide which cell donates its q value.
double BurgersGodunovSolver2d::compute_transport_flux(double wL, double wR,
                                                      double qL,
                                                      double qR) const {
  if (wL > wR) {
    // Shock in w field
    double s = 0.5 * (wL + wR);
    return (s > 0.0) ? wL * qL : wR * qR;
  } else {
    // Rarefaction in w field
    if (wL >= 0.0)
      return wL * qL;
    else if (wR <= 0.0)
      return wR * qR;
    else
      return 0.0; // sonic: w* = 0
  }
}

// ============================================================
//  Single time step
// ============================================================
void BurgersGodunovSolver2d::step(double dt) {
  const double dtdx = dt / dx_;
  const double dtdy = dt / dy_;

  // ---- x-direction fluxes at interfaces (i+1/2, j) ----
  // F_u[i,j] = flux of u across right face of cell (i,j)
  // F_v[i,j] = flux of v across right face of cell (i,j)
  // Size: Nx * Ny, layout identical to u_/v_
  std::vector<double> Fu(Nx_ * Ny_), Fv(Nx_ * Ny_);

  // i outer, j inner  ->  sequential memory (cache-friendly)
  for (int i = 0; i < Nx_; ++i) {
    int ip1 = wrapI(i + 1); // periodic right neighbour
    for (int j = 0; j < Ny_; ++j) {
      double uL = u_[idx(i, j)];
      double uR = u_[idx(ip1, j)];
      double vL = v_[idx(i, j)];
      double vR = v_[idx(ip1, j)];

      // u-equation x-flux: (u^2/2)_x  -> self-advection
      Fu[idx(i, j)] = compute_flux(uL, uR);
      // v-equation x-flux: (u*v)_x  -> transport of v by u-wave
      Fv[idx(i, j)] = compute_transport_flux(uL, uR, vL, vR);
    }
  }

  // ---- y-direction fluxes at interfaces (i, j+1/2) ----
  // G_u[i,j] = flux of u across top face of cell (i,j)
  // G_v[i,j] = flux of v across top face of cell (i,j)
  std::vector<double> Gu(Nx_ * Ny_), Gv(Nx_ * Ny_);

  for (int i = 0; i < Nx_; ++i) {
    for (int j = 0; j < Ny_; ++j) {
      int jp1 = wrapJ(j + 1); // periodic top neighbour

      double uB = u_[idx(i, j)];   // u below interface
      double uT = u_[idx(i, jp1)]; // u above interface
      double vB = v_[idx(i, j)];
      double vT = v_[idx(i, jp1)];

      // u-equation y-flux: (u*v)_y  -> transport of u by v-wave
      Gu[idx(i, j)] = compute_transport_flux(vB, vT, uB, uT);
      // v-equation y-flux: (v^2/2)_y -> self-advection
      Gv[idx(i, j)] = compute_flux(vB, vT);
    }
  }

  // ---- Unsplit conservative update ----
  //   q^{n+1}_{i,j} = q^n_{i,j}
  //       - (dt/dx)(F_{i+1/2,j} - F_{i-1/2,j})
  //       - (dt/dy)(G_{i,j+1/2} - G_{i,j-1/2})
  //
  // F_{i-1/2,j} = Fu[idx(im1, j)]   (stored on left cell)
  // G_{i,j-1/2} = Gu[idx(i, jm1)]   (stored on bottom cell)
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

  std::swap(u_, u_next_);
  std::swap(v_, v_next_);
}

// ============================================================
//  Solve
// ============================================================
void BurgersGodunovSolver2d::solve() {
  if (!ic_u_ || !ic_v_) {
    throw std::runtime_error(
        "Initial conditions not set. Call setInitialConditions() first.");
  }

  std::cout << "BurgersGodunovSolver2d: starting solve\n"
            << "  Grid: " << Nx_ << " x " << Ny_ << "  dx=" << dx_
            << "  dy=" << dy_ << "\n"
            << "  t_end=" << t_end_ << "  CFL=" << CFL_ << "\n";

  // Save snapshot at t=0
  history_.clear();
  history_.push_back({0.0, u_, v_});

  double t = 0.0;
  double t_next_save = save_interval_;
  int step_count = 0;

  // Precompute initial conservation reference
  double sum_u0 = 0.0, sum_v0 = 0.0;
  for (int k = 0; k < Nx_ * Ny_; ++k) {
    sum_u0 += u_[k];
    sum_v0 += v_[k];
  }
  std::cout << std::fixed << std::setprecision(8);
  std::cout << "  Initial conservation:  sum(u)=" << sum_u0
            << "  sum(v)=" << sum_v0 << "\n\n";

  while (t < t_end_) {
    // ---- CFL-based dt ----
    double max_speed = 0.0;
    for (int k = 0; k < Nx_ * Ny_; ++k) {
      max_speed =
          std::max(max_speed, std::max(std::abs(u_[k]), std::abs(v_[k])));
    }

    if (max_speed < 1e-14) {
      std::cout << "  WARNING: max wave speed is effectively zero. Stopping.\n";
      break;
    }

    // 2D unsplit CFL: dt <= CFL * min(dx,dy) / max_speed
    double min_h = std::min(dx_, dy_);
    double dt = CFL_ * min_h / max_speed;

    // Don't overshoot t_end or the next save point
    dt = std::min(dt, t_end_ - t);
    dt = std::min(dt, t_next_save - t);

    // ---- CFL Monitor ----
    double cfl_x = dt * max_speed / dx_;
    double cfl_y = dt * max_speed / dy_;
    if (cfl_x > 1.0 || cfl_y > 1.0) {
      std::cerr << "  CFL VIOLATION at t=" << t << "  cfl_x=" << cfl_x
                << "  cfl_y=" << cfl_y << "\n";
    }

    // ---- Advance ----
    step(dt);
    t += dt;
    ++step_count;

    // ---- Conservation Monitor (every 100 steps) ----
    if (step_count % 100 == 0) {
      double sum_u = 0.0, sum_v = 0.0;
      for (int k = 0; k < Nx_ * Ny_; ++k) {
        sum_u += u_[k];
        sum_v += v_[k];
      }
      std::cout << "  t=" << std::setw(10) << t << "  step=" << std::setw(6)
                << step_count << "  sum(u)=" << sum_u
                << "  |err_u|=" << std::abs(sum_u - sum_u0)
                << "  sum(v)=" << sum_v
                << "  |err_v|=" << std::abs(sum_v - sum_v0)
                << "  max_speed=" << max_speed << "\n";
    }

    // ---- Save snapshot ----
    if (t >= t_next_save - 1e-14) {
      history_.push_back({t, u_, v_});
      t_next_save += save_interval_;
    }
  }

  // Always save the final state
  if (history_.empty() || history_.back().t < t - 1e-14) {
    history_.push_back({t, u_, v_});
  }

  // Final conservation check
  double sum_u = 0.0, sum_v = 0.0;
  for (int k = 0; k < Nx_ * Ny_; ++k) {
    sum_u += u_[k];
    sum_v += v_[k];
  }
  std::cout << "\nSolve complete.  Steps=" << step_count << "  t_final=" << t
            << "\n"
            << "  Final conservation:  sum(u)=" << sum_u
            << "  |err_u|=" << std::abs(sum_u - sum_u0) << "\n"
            << "                       sum(v)=" << sum_v
            << "  |err_v|=" << std::abs(sum_v - sum_v0) << "\n"
            << "  Snapshots saved: " << history_.size() << "\n\n";
}

// ============================================================
//  Save to CSV
// ============================================================
void BurgersGodunovSolver2d::saveSolution(const std::string &base_folder,
                                          const std::string &run_name) const {
  namespace fs = std::filesystem;

  if (!fs::exists(base_folder)) {
    fs::create_directories(base_folder);
  }

  std::string run_folder = base_folder + "/" + run_name;
  if (fs::exists(run_folder)) {
    std::cout << "Overwriting existing folder: " << run_folder << "\n";
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
    const auto &snap_u = history_[snap].u;
    const auto &snap_v = history_[snap].v;

    // i outer, j inner -> sequential read (cache-friendly)
    for (int i = 0; i < Nx_; ++i) {
      double x = (i + 0.5) * dx_;
      for (int j = 0; j < Ny_; ++j) {
        double y = (j + 0.5) * dy_;
        int k = i * Ny_ + j;
        file << x << "," << y << "," << snap_u[k] << "," << snap_v[k] << "\n";
      }
    }
  }

  std::cout << "Done.\n";
}
