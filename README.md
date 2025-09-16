# VIP ExaFlow – Burgers + Libtorch Setup

## ✅ Completed (as of [today’s date])
- **Libtorch (C++ API) installed and tested**
  - Verified on **Mac** and **ICE cluster** with `torch_cpp_test` (random tensor + matmul).
  - Environment variables (`LIBTORCH`, `CMAKE_PREFIX_PATH`, `LD_LIBRARY_PATH`) set and working.

- **Baseline 1D Burgers solver (C++)**
  - Periodic BCs, Rusanov flux, adaptive `dt`.
  - Runs stably on ICE (CFL < 0.4, mean ≈ 0, L2 norm decays).
  - Snapshots written as CSV for plotting.
  - Slurm script (`run_burgers.slurm`) to submit jobs on ICE.

- **Neural Network demo (C++)**
  - MLP (`mlp_demo`) implemented in LibTorch.
  - Successfully trains to approximate `sin(x)` (loss decreases, predictions reasonable).
  - Confirms NN forward/backward + optimizer all working in C++.

- **Coupled NN + Burgers (`burgers_nn`)**
  - Added `WeightPredictor` stub (MLP, 5-point stencil → 3 weights).
  - Integrated into solver: per-interface dissipation scaling modulated by NN output.
  - Compiles + runs cleanly on ICE.
  - Produces stable results with NN-coupled dissipation (near baseline, as weights ~uniform).
  - Ready for TorchScript model loading when a trained model is available.

## 📂 Project Structure

vip/
├─ libtorch_test/   # initial torch demo + libtorch install
│   ├─ test.cpp
│   ├─ CMakeLists.txt
│   └─ run_torch.slurm
├─ burgers/         # baseline + NN-coupled solvers
│   ├─ burgers.hpp / burgers.cpp   # solver core
│   ├─ main.cpp     # baseline Burgers
│   ├─ main_nn.cpp  # NN-coupled Burgers
│   ├─ mlp.hpp      # tiny MLP definition
│   ├─ nn_weights.hpp # weight predictor stub
│   ├─ CMakeLists.txt
│   ├─ run_burgers.slurm
│   └─ run_burgers_nn.slurm
└─ nn_demo/         # standalone NN training demo
├─ mlp.hpp
├─ mlp_demo.cpp
├─ CMakeLists.txt
└─ run_mlp.slurm

## 🔜 Next Steps
1. **Diagnostics**
   - Add matching CFL/mean/L2 prints in `burgers_nn` (like baseline).
   - Optionally save CSV snapshots (`snap_XXXX.csv`) from both runs for direct comparisons.

2. **TorchScript Integration**
   - Train WENO-weight predictor in Python.
   - Export via `torch.jit.script` → `weno_weights.pt`.
   - Load in C++ (`torch::jit::load`) inside `WeightPredictor`.

3. **Training Data**
   - Generate stencils/labels from smooth + shock test cases.
   - Phase 1: sanity (learn uniform weights).
   - Phase 2: learn WENO-5 weights or flux limiter proxy.

4. **Safety**
   - Add ENO fallback + renormalization of weights at inference.
   - Blend NN effect with baseline (`tau` knob) for stable rollout.

5. **Presentation Prep**
   - Plot baseline vs. NN-coupled results.
   - Show alpha scaling stats (min/mean/max).
   - Include README + Slurm scripts for reproducibility.

---

## 🚀 Usage

**Build (baseline + NN)**
```bash
cd burgers
export LIBTORCH="$HOME/vip/libtorch_test/libtorch"
export CMAKE_PREFIX_PATH="$LIBTORCH"
export LD_LIBRARY_PATH="$LIBTORCH/lib:$LD_LIBRARY_PATH"

cmake -S . -B build
cmake --build build -j


Run

./build/burgers     # baseline
./build/burgers_nn  # NN-coupled

Run with Slurm

sbatch run_burgers.slurm
sbatch run_burgers_nn.slurm