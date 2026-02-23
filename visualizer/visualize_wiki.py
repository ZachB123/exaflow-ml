import argparse
import glob
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

def load_frames(data_dir):
    """
    Loads all timestep_XXXXX.csv files inside data_dir.
    Returns:
        x: np.array of shape (N,)
        frames: list of np.array, each shape (N,)
    """
    files = sorted(glob.glob(os.path.join(data_dir, "timestep_*.csv")))

    if not files:
        raise ValueError(f"No timestep_*.csv files found in {data_dir}")

    frames = []
    x = None

    for f in files:
        data = np.loadtxt(f, delimiter=",", skiprows=1)
        if x is None:
            x = data[:, 0] - 6.0  # Shift to [-6, 6]
        u = data[:, 1]
        # Replace |u| > 5.0 with NaN to prevent axis issues
        u = np.where(np.abs(u) > 5.0, np.nan, u)
        frames.append(u)

    return x, frames, files

def run_visualizer(case, initial_speed, scheme):
    viscosities = ["1.0", "0.1", "0.01"]
    colors = {"1.0": "black", "0.1": "red", "0.01": "blue"}
    labels = {"1.0": "v=1.0", "0.1": "v=0.1", "0.01": "v=0.01"}

    # Load data for each viscosity
    data = {}
    max_frames = 0
    for v in viscosities:
        folder_name = f"wiki_{case}_{v}_{scheme}"
        data_dir = os.path.join(PROJECT_ROOT, "data_v2", folder_name)
        x, frames, files = load_frames(data_dir)
        data[v] = {"x": x, "frames": frames, "files": files}
        max_frames = max(max_frames, len(frames))

    # Ensure all have the same number of frames (pad with last frame if needed)
    for v in viscosities:
        while len(data[v]["frames"]) < max_frames:
            data[v]["frames"].append(data[v]["frames"][-1])

    fig, ax = plt.subplots()
    fig.canvas.manager.set_window_title(f"Wikipedia Burgers Visualizer – {case}")
    plt.subplots_adjust(bottom=0.25)

    # Plot initial frame for all viscosities
    lines = {}
    for v in viscosities:
        line, = ax.plot(data[v]["x"], data[v]["frames"][0], color=colors[v], label=labels[v])
        lines[v] = line

    ax.set_xlim([-7, 7])
    ax.set_ylim([-1.1, 1.1])
    ax.set_title(f"{case.capitalize()} Initial Condition (t=0)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    frame_pos = 0.0
    frame_idx = 0
    playing = False
    speed = initial_speed

    def update_plot():
        """
        Update the plotted lines for all viscosities.
        """
        for v in viscosities:
            lines[v].set_ydata(data[v]["frames"][frame_idx])
        ax.set_title(f"{case.capitalize()} (frame {frame_idx})")
        fig.canvas.draw_idle()

    # Button functions
    def next_frame(event):
        nonlocal frame_pos, frame_idx
        frame_idx = (frame_idx + 1) % max_frames
        frame_pos = float(frame_idx)
        update_plot()

    def prev_frame(event):
        nonlocal frame_pos, frame_idx
        frame_idx = (frame_idx - 1) % max_frames
        frame_pos = float(frame_idx)
        update_plot()

    def play_pause(event):
        nonlocal playing, frame_pos, frame_idx

        # If we are at the last frame, reset to 0
        if frame_idx >= max_frames - 1:
            frame_pos = 0.0
            frame_idx = 0
            update_plot()

        playing = not playing
        play_button.label.set_text("Pause" if playing else "Play")
        fig.canvas.draw_idle()

    # Slider: speed override
    def change_speed(val):
        nonlocal speed
        speed = float(val)

    # --- UI Layout ---
    axprev = plt.axes([0.1, 0.1, 0.1, 0.075])
    axplay = plt.axes([0.25, 0.1, 0.15, 0.075])
    axnext = plt.axes([0.45, 0.1, 0.1, 0.075])
    axspeed = plt.axes([0.65, 0.1, 0.25, 0.05])

    prev_button = widgets.Button(axprev, "Prev")
    play_button = widgets.Button(axplay, "Play")
    next_button = widgets.Button(axnext, "Next")
    speed_slider = widgets.Slider(axspeed, "Speed", 0.1, 50.0, valinit=speed)

    prev_button.on_clicked(prev_frame)
    next_button.on_clicked(next_frame)
    play_button.on_clicked(play_pause)
    speed_slider.on_changed(change_speed)

    timer = fig.canvas.new_timer(interval=30)  # milliseconds

    def on_timer():
        nonlocal frame_pos, frame_idx, playing
        if not playing:
            return

        frame_pos += speed

        # Reached or passed the last frame
        if frame_pos >= max_frames - 1:
            frame_pos = max_frames - 1
            frame_idx = int(frame_pos)
            update_plot()

            # Stop playback
            playing = False
            play_button.label.set_text("Play")
            fig.canvas.draw_idle()
            return

        # Normal update
        frame_idx = int(frame_pos)
        update_plot()

    timer.add_callback(on_timer)
    timer.start()

    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Visualize Wikipedia Burgers' equation data")
    parser.add_argument(
        "case",
        type=str,
        choices=["gaussian", "n_wave"],
        help="Initial condition case: gaussian or n_wave",
    )
    parser.add_argument(
        "--speed",
        type=float,
        default=1.0,
        help="Initial playback speed (frames per tick)",
    )
    parser.add_argument(
        "--scheme",
        type=str,
        default="LaxWendroff",
        choices=["LaxWendroff", "Godunov"],
        help="Numerical scheme used for the simulation",
    )

    args = parser.parse_args()
    run_visualizer(args.case, args.speed, args.scheme)

if __name__ == "__main__":
    main()