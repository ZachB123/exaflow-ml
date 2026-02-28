#pragma once

#include <functional>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Configuration for the 2D Godunov Burgers solver
// ---------------------------------------------------------------------------
struct Solver2dConfig {
  int Nx;               // number of cells in x
  int Ny;               // number of cells in y
  double Lx;            // domain length in x  [0, Lx]
  double Ly;            // domain length in y  [0, Ly]
  double CFL;           // CFL number (<=0.45 recommended for 2D unsplit)
  double t_end;         // final simulation time
  double save_interval; // physical time between saved snapshots
};

// ---------------------------------------------------------------------------
// 2D Unsplit Godunov Solver for the inviscid Burgers system
//
//   u_t + (u^2/2)_x + (uv)_y = 0
//   v_t + (uv)_x  + (v^2/2)_y = 0
//
// Finite-volume, first-order, periodic BCs.
// Data layout: flat vector, idx(i,j) = i*Ny + j  (row-major, cache-friendly).
// All loops are structured as  for(i ...) for(j ...)  to keep j as the
// fast-varying index and maximise cache-line utilisation.
// ---------------------------------------------------------------------------
class BurgersGodunovSolver2d {
public:
  BurgersGodunovSolver2d(const Solver2dConfig &config);

  // Set initial conditions via two separate lambdas: u0(x,y) and v0(x,y).
  void setInitialConditions(const std::function<double(double, double)> &u0,
                            const std::function<double(double, double)> &v0);

  // Run the solver to t_end.  Prints conservation + CFL diagnostics each step.
  void solve();

  // Write all saved snapshots to <base_folder>/<run_name>/timestep_XXXXX.csv
  // Each CSV has columns: x,y,u,v
  void saveSolution(const std::string &base_folder,
                    const std::string &run_name) const;

private:
  // ------- config -------
  int Nx_, Ny_;
  double Lx_, Ly_;
  double dx_, dy_;
  double CFL_;
  double t_end_;
  double save_interval_;

  // ------- state -------
  std::vector<double> u_, v_; // current fields
  std::vector<double> u_next_, v_next_;

  // ------- history -------
  // Each snapshot is {t, u[Nx*Ny], v[Nx*Ny]}
  struct Snapshot {
    double t;
    std::vector<double> u;
    std::vector<double> v;
  };
  std::vector<Snapshot> history_;

  // ------- IC store -------
  std::function<double(double, double)> ic_u_;
  std::function<double(double, double)> ic_v_;

  // ------- helpers -------
  // Flat index (row-major, j fast)
  inline int idx(int i, int j) const { return i * Ny_ + j; }

  // Periodic wrap
  inline int wrapI(int i) const { return (i + Nx_) % Nx_; }
  inline int wrapJ(int j) const { return (j + Ny_) % Ny_; }

  // Godunov flux for self-advection: f(q) = q^2/2  with wave speed = q
  double compute_flux(double qL, double qR) const;

  // Upwind transport flux: transports quantity q with wave speed w
  // Flux = w * q  using donor-cell upwinding on the w Riemann problem
  double compute_transport_flux(double wL, double wR, double qL,
                                double qR) const;

  // One time-step update (dt already computed)
  void step(double dt);
};
