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
            x = data[:, 0]
        frames.append(data[:, 1])

    return x, frames, files


def run_visualizer(folder_name, initial_speed):
    data_dir = os.path.join(PROJECT_ROOT, "data", folder_name)
    x, frames, files = load_frames(data_dir)

    fig, ax = plt.subplots()
    fig.canvas.manager.set_window_title(f"Burgers Visualizer – {folder_name}")
    plt.subplots_adjust(bottom=0.25)

    line, = ax.plot(x, frames[0])
    ax.set_title(os.path.basename(files[0]))

    frame_pos = 0.0
    frame_idx = 0
    playing = False
    speed = initial_speed

    y_lock = True            # False = autoscale every frame; True = keep current y-limits
    y_margin_frac = 0.05      # 5% margin above/below data range when autoscaling

    def update_plot():
        """
        Update the plotted line and optionally autoscale the y-axis.
        Uses frame_idx (int) and frames (list of arrays) from outer scope.
        """
        line.set_ydata(frames[frame_idx])
        ax.set_title(f"{os.path.basename(files[frame_idx])}   (frame {frame_idx})")

        if not y_lock:
            # compute min/max of the current frame
            y_min = float(np.min(frames[frame_idx]))
            y_max = float(np.max(frames[frame_idx]))

            # avoid zero-height limits
            if y_max == y_min:
                # expand by +-0.5 or 10% if value is large
                delta = max(abs(y_min) * 0.1, 0.5)
                y_min -= delta
                y_max += delta

            # add a small margin so plot doesn't touch the edges
            margin = y_margin_frac * (y_max - y_min)
            ax.set_ylim(y_min - margin, y_max + margin)

        # if y_lock is True we do nothing (y-limits stay as they are)
        fig.canvas.draw_idle()

    # --- Checkbox to lock/unlock y-axis ---
    # Place it near the right side of the figure (adjust coords if overlap)
    ax_lock = plt.axes([0.82, 0.15, 0.12, 0.08])  # [left, bottom, width, height] in figure coords
    lock_cb = widgets.CheckButtons(ax_lock, ["Lock Y"], [y_lock])

    def on_toggle_lock(label):
        nonlocal y_lock
        # CheckButtons toggles label's state; simply flip y_lock
        y_lock = not y_lock

        # If we just locked the y-axis, capture current limits as the locked limits
        if y_lock:
            # read current limits (use what is visible now)
            current_ylim = ax.get_ylim()
            ax.set_ylim(current_ylim)  # enforce them explicitly
        else:
            # if unlocking, do an immediate autoscale to current frame so view updates
            # (this uses update_plot's autoscale branch)
            update_plot()

        fig.canvas.draw_idle()

    lock_cb.on_clicked(on_toggle_lock)

    # button: next frame
    def next_frame(event):
        nonlocal frame_pos, frame_idx
        frame_idx = (frame_idx + 1) % len(frames)
        frame_pos = float(frame_idx)
        update_plot()

    def prev_frame(event):
        nonlocal frame_pos, frame_idx
        frame_idx = (frame_idx - 1) % len(frames)
        frame_pos = float(frame_idx)
        update_plot()

    def play_pause(event):
        nonlocal playing, frame_pos, frame_idx

        # If we are at the last frame, reset to 0
        if frame_idx >= len(frames) - 1:
            frame_pos = 0.0
            frame_idx = 0
            update_plot()

        playing = not playing
        play_button.label.set_text("Pause" if playing else "Play")
        fig.canvas.draw_idle()


    # slider: speed override
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

        # reached or passed the last frame
        if frame_pos >= len(frames) - 1:
            frame_pos = len(frames) - 1
            frame_idx = int(frame_pos)
            update_plot()

            # stop playback
            playing = False
            play_button.label.set_text("Play")
            fig.canvas.draw_idle()
            return

        # normal update
        frame_idx = int(frame_pos)
        update_plot()

    timer.add_callback(on_timer)
    timer.start()

    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Visualize Burgers' equation data")
    parser.add_argument(
        "folder",
        type=str,
        help="Folder name inside data/, e.g. step_function or sine_wave",
    )
    parser.add_argument(
        "--speed",
        type=int,
        default=1,
        help="Initial playback speed (frames per tick)",
    )

    args = parser.parse_args()
    run_visualizer(args.folder, args.speed)


if __name__ == "__main__":
    main()
