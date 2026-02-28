"""
visualize_godunov2d.py
======================
Interactive animated visualizer for the 2D Godunov Burgers solver output.

Usage:
    python visualize_godunov2d.py <run_folder> [--speed N]

Examples:
    python visualize_godunov2d.py data_2d/godunov_squarewave
    python visualize_godunov2d.py data_2d/godunov_gaussian --speed 3

<run_folder> can be:
  - an absolute path
  - a path relative to the project root
  - a folder name under  data_2d/

Each CSV (timestep_XXXXX.csv) must have columns: x, y, u, v

Layout (three panels):
  Left:   pcolormesh heatmap of  u  field (colormap locked to IC range)
          with a quiver overlay showing the (u, v) vector field
  Right:  cross-section line plot of  u  along the mid-y slice vs. x

Controls (matching visualize.py style):
  Play / Pause / Prev / Next buttons
  Speed slider  (frames per timer tick)
  Lock Y checkbox  (cross-section panel)
"""

import argparse
import glob
import os

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets

# -----------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)


# -----------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------
def locate_run_folder(folder_name: str) -> str:
    """Resolve run_folder: absolute, relative, or under data_2d/."""
    candidates = [
        folder_name,
        os.path.join(PROJECT_ROOT, folder_name),
        os.path.join(PROJECT_ROOT, "data_2d", folder_name),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return os.path.abspath(c)
    raise ValueError(f"Could not find run folder: '{folder_name}'\n"
                     f"Tried: {candidates}")


def load_frames(data_dir: str):
    """
    Load all timestep_XXXXX.csv files.

    Returns
    -------
    xs, ys  : 1-D sorted unique coordinate arrays
    U_frames, V_frames : list of 2-D arrays shaped (Nx, Ny)
    file_labels : list of strings (basename for each frame)
    """
    files = sorted(glob.glob(os.path.join(data_dir, "timestep_*.csv")))
    if not files:
        raise ValueError(f"No timestep_*.csv files found in {data_dir}")

    print(f"Loading {len(files)} snapshots from {data_dir} …", flush=True)

    # Read first file to determine grid dimensions
    first = np.loadtxt(files[0], delimiter=",", skiprows=1)
    xs_raw = first[:, 0]
    ys_raw = first[:, 1]

    xs = np.unique(xs_raw)
    ys = np.unique(ys_raw)
    Nx, Ny = len(xs), len(ys)

    U_frames, V_frames = [], []

    for f in files:
        data = np.loadtxt(f, delimiter=",", skiprows=1)
        # data rows are in idx(i,j) = i*Ny+j order (row-major, j fast)
        U = data[:, 2].reshape(Nx, Ny)
        V = data[:, 3].reshape(Nx, Ny)
        U_frames.append(U)
        V_frames.append(V)

    labels = [os.path.basename(f) for f in files]
    return xs, ys, U_frames, V_frames, labels


# -----------------------------------------------------------------------
# Main visualizer
# -----------------------------------------------------------------------
def run_visualizer(folder_name: str, initial_speed: int):
    data_dir = locate_run_folder(folder_name)
    xs, ys, U_frames, V_frames, labels = load_frames(data_dir)

    Nx, Ny   = len(xs), len(ys)
    mid_j    = Ny // 2          # index of mid-y slice
    n_frames = len(U_frames)

    # Quiver subsampling: show every K-th cell so arrows are readable
    K   = max(1, min(Nx, Ny) // 20)
    qi  = np.arange(0, Nx, K)
    qj  = np.arange(0, Ny, K)
    qX, qY = np.meshgrid(xs[qi], ys[qj], indexing="ij")

    # Colormap range locked to initial condition
    u_min = float(U_frames[0].min())
    u_max = float(U_frames[0].max())
    if u_max == u_min:
        u_min -= 0.5
        u_max += 0.5

    # ---- Figure layout ----
    fig = plt.figure(figsize=(13, 6))
    fig.canvas.manager.set_window_title(
        f"Godunov 2D Visualizer – {os.path.basename(data_dir)}")
    plt.subplots_adjust(bottom=0.28, left=0.06, right=0.98, wspace=0.35)

    ax_heat  = fig.add_subplot(1, 2, 1)
    ax_slice = fig.add_subplot(1, 2, 2)

    # -- Heatmap (pcolormesh) --
    # pcolormesh expects (Ny+1) x (Nx+1) edge arrays; we use cell centres here.
    dx = xs[1] - xs[0] if Nx > 1 else 1.0
    dy = ys[1] - ys[0] if Ny > 1 else 1.0
    x_edges = np.append(xs - dx / 2, xs[-1] + dx / 2)
    y_edges = np.append(ys - dy / 2, ys[-1] + dy / 2)
    Xe, Ye  = np.meshgrid(x_edges, y_edges, indexing="ij")

    pcm = ax_heat.pcolormesh(
        Xe, Ye, U_frames[0],
        cmap="RdBu_r", vmin=u_min, vmax=u_max, shading="flat")
    fig.colorbar(pcm, ax=ax_heat, label="u")
    ax_heat.set_xlabel("x")
    ax_heat.set_ylabel("y")
    ax_heat.set_aspect("equal")

    # -- Quiver overlay --
    quiv = ax_heat.quiver(
        qX, qY,
        U_frames[0][np.ix_(qi, qj)],
        V_frames[0][np.ix_(qi, qj)],
        color="k", alpha=0.55, scale=None,
        width=0.003, headwidth=4)
    heat_title = ax_heat.set_title(labels[0])

    # -- Cross-section slice --
    (slice_line,) = ax_slice.plot(xs, U_frames[0][:, mid_j], color="steelblue", lw=1.5)
    ax_slice.set_xlabel("x")
    ax_slice.set_ylabel(f"u at y = {ys[mid_j]:.3f}  (mid-y slice)")
    ax_slice.set_title("Cross-section slice")
    ax_slice.set_xlim(xs[0], xs[-1])
    slice_title = ax_slice.set_title(f"Mid-y slice  (frame 0)")

    # Lock y-limits to initial slice
    y_slice_min = float(U_frames[0][:, mid_j].min())
    y_slice_max = float(U_frames[0][:, mid_j].max())
    margin = 0.1 * max(abs(y_slice_max - y_slice_min), 1e-6)
    ax_slice.set_ylim(y_slice_min - margin, y_slice_max + margin)

    # ---- Shared state ----
    state = dict(frame_idx=0, frame_pos=0.0, playing=False,
                 speed=float(initial_speed), y_lock=True)

    # ---- Update function ----
    def update_plot():
        fi = state["frame_idx"]
        U = U_frames[fi]
        V = V_frames[fi]

        # heatmap
        pcm.set_array(U.ravel())
        heat_title.set_text(labels[fi])

        # quiver
        quiv.set_UVC(
            U[np.ix_(qi, qj)].ravel(),
            V[np.ix_(qi, qj)].ravel())

        # slice
        slice_data = U[:, mid_j]
        slice_line.set_ydata(slice_data)
        slice_title.set_text(f"Mid-y slice  (frame {fi})")

        if not state["y_lock"]:
            s_min = float(slice_data.min())
            s_max = float(slice_data.max())
            if s_max == s_min:
                s_min -= 0.5; s_max += 0.5
            mg = 0.05 * (s_max - s_min)
            ax_slice.set_ylim(s_min - mg, s_max + mg)

        fig.canvas.draw_idle()

    # ---- Button callbacks ----
    def next_frame(_):
        state["frame_idx"] = (state["frame_idx"] + 1) % n_frames
        state["frame_pos"] = float(state["frame_idx"])
        update_plot()

    def prev_frame(_):
        state["frame_idx"] = (state["frame_idx"] - 1) % n_frames
        state["frame_pos"] = float(state["frame_idx"])
        update_plot()

    def play_pause(_):
        if state["frame_idx"] >= n_frames - 1:
            state["frame_idx"] = 0
            state["frame_pos"] = 0.0
            update_plot()
        state["playing"] = not state["playing"]
        play_btn.label.set_text("Pause" if state["playing"] else "Play")
        fig.canvas.draw_idle()

    def change_speed(val):
        state["speed"] = float(val)

    def toggle_lock(_):
        state["y_lock"] = not state["y_lock"]
        if state["y_lock"]:
            ax_slice.set_ylim(ax_slice.get_ylim())
        else:
            update_plot()
        fig.canvas.draw_idle()

    # ---- Timer ----
    def on_timer():
        if not state["playing"]:
            return
        state["frame_pos"] += state["speed"]
        if state["frame_pos"] >= n_frames - 1:
            state["frame_pos"] = n_frames - 1
            state["frame_idx"] = int(state["frame_pos"])
            update_plot()
            state["playing"] = False
            play_btn.label.set_text("Play")
            fig.canvas.draw_idle()
            return
        state["frame_idx"] = int(state["frame_pos"])
        update_plot()

    timer = fig.canvas.new_timer(interval=30)
    timer.add_callback(on_timer)
    timer.start()

    # ---- Widget layout ----
    ax_prev   = plt.axes([0.10, 0.12, 0.09, 0.06])
    ax_play   = plt.axes([0.21, 0.12, 0.12, 0.06])
    ax_next   = plt.axes([0.35, 0.12, 0.09, 0.06])
    ax_speed  = plt.axes([0.50, 0.14, 0.25, 0.04])
    ax_lock   = plt.axes([0.80, 0.10, 0.12, 0.08])

    prev_btn    = widgets.Button(ax_prev,  "Prev")
    play_btn    = widgets.Button(ax_play,  "Play")
    next_btn    = widgets.Button(ax_next,  "Next")
    speed_slider= widgets.Slider(ax_speed, "Speed", 0.1, 20.0,
                                 valinit=float(initial_speed))
    lock_cb     = widgets.CheckButtons(ax_lock, ["Lock Y"], [state["y_lock"]])

    prev_btn.on_clicked(prev_frame)
    next_btn.on_clicked(next_frame)
    play_btn.on_clicked(play_pause)
    speed_slider.on_changed(change_speed)
    lock_cb.on_clicked(toggle_lock)

    plt.show()


# -----------------------------------------------------------------------
# CLI entry point
# -----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Interactive 2D Godunov Burgers visualizer")
    parser.add_argument(
        "folder",
        type=str,
        help="Run folder (absolute path, relative, or name under data_2d/)",
    )
    parser.add_argument(
        "--speed",
        type=int,
        default=1,
        help="Initial playback speed (frames per tick, default: 1)",
    )
    args = parser.parse_args()
    run_visualizer(args.folder, args.speed)


if __name__ == "__main__":
    main()
