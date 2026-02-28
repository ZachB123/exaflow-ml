#pragma once

#include <functional>
#include <string>
#include <vector>

// All the settings you need to configure a 2D run.
struct Solver2dConfig {
  int Nx;               // grid cells in x
  int Ny;               // grid cells in y
  double Lx;            // domain width  [0, Lx]
  double Ly;            // domain height [0, Ly]
  double CFL;           // time step safety factor, keep <= 0.45 for 2D
  double t_end;         // stop time
  double save_interval; // how often (in sim-time) to write a snapshot to disk
};

// Solves the inviscid 2D Burgers system:
//   u_t + (u^2/2)_x + (u*v)_y = 0
//   v_t + (u*v)_x  + (v^2/2)_y = 0
//
// First-order Godunov finite-volume scheme, unsplit, periodic boundaries.
// Grid data is stored flat: index(i,j) = i*Ny + j (row-major).
// All loops run i-outer / j-inner so the inner index is the fast one,
// which keeps memory accesses sequential and cache-friendly.
class BurgersGodunovSolver2d {
public:
  BurgersGodunovSolver2d(const Solver2dConfig &config);

  // Set u(x,y,0) and v(x,y,0) via two lambdas.
  void setInitialConditions(const std::function<double(double, double)> &u0,
                            const std::function<double(double, double)> &v0);

  // Run solver to t_end. Prints conservation and CFL info along the way.
  void solve();

  // Write snapshots to <base_folder>/<run_name>/timestep_XXXXX.csv
  // Columns: x, y, u, v
  void saveSolution(const std::string &base_folder,
                    const std::string &run_name) const;

private:
  int Nx_, Ny_;
  double Lx_, Ly_;
  double dx_, dy_;
  double CFL_;
  double t_end_;
  double save_interval_;

  std::vector<double> u_, v_;           // current fields, size Nx*Ny
  std::vector<double> u_next_, v_next_; // scratch buffers for the next step

  // One saved snapshot: time value + full u and v arrays at that time.
  struct Snapshot {
    double t;
    std::vector<double> u;
    std::vector<double> v;
  };
  std::vector<Snapshot> history_;

  std::function<double(double, double)> ic_u_;
  std::function<double(double, double)> ic_v_;

  // Flat 2D index, j is the fast (inner) dimension.
  inline int idx(int i, int j) const { return i * Ny_ + j; }

  // Periodic wrap so we don't need special boundary code.
  inline int wrapI(int i) const { return (i + Nx_) % Nx_; }
  inline int wrapJ(int j) const { return (j + Ny_) % Ny_; }

  // Godunov flux for a scalar self-advection term: f(q) = q^2/2.
  double compute_flux(double qL, double qR) const;

  // Upwind flux for a transport term: flux = w*q.
  // w is the wave speed, q is the thing being transported.
  double compute_transport_flux(double wL, double wR, double qL,
                                double qR) const;

  void step(double dt);
};
