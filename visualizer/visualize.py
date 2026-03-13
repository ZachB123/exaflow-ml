import argparse
import json
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

SEARCH_DIRS = ["training_data", "data"]

# maximum number of spatial points passed to matplotlib per frame;
# grids larger than this are downsampled by a uniform stride before display
MAX_DISPLAY_POINTS = 4000


def load_frames(data_dir):
    """
    Loads Burgers solution from metadata.json and solution.bin

    Returns:
        x: np.array of shape (N,)
        frames: np.memmap of shape (T, N)
        frame_labels: list[str] of length T
    """
    metadata_path = os.path.join(data_dir, "metadata.json")
    solution_path = os.path.join(data_dir, "solution.bin")

    if not os.path.exists(metadata_path):
        raise ValueError(f"metadata.json not found in {data_dir}")

    if not os.path.exists(solution_path):
        raise ValueError(f"solution.bin not found in {data_dir}")

    with open(metadata_path, "r") as f:
        metadata = json.load(f)

    solver = metadata["solver"]
    time_steps = int(solver["time_steps"])
    num_domain_points = int(solver["num_domain_points"])
    spatial_step_size = float(solver["spatial_step_size"])

    x = spatial_step_size * np.arange(num_domain_points)

    frames = np.memmap(
        solution_path,
        dtype=np.float64,
        mode="r",
        shape=(time_steps, num_domain_points),
    )

    frame_labels = [f"timestep_{i:05d}" for i in range(time_steps)]

    return x, frames, frame_labels


def resolve_folder(folder_name):
    # Resolves the folder path for a given sample name or relative path.
    # If folder_name is an integer, treat as sample_XXXXXX in training_data
    if folder_name.isdigit():
        sample_name = f"sample_{int(folder_name):06d}"
        candidate = os.path.join(PROJECT_ROOT, "training_data", sample_name)
        if os.path.isdir(candidate):
            return candidate
        raise FileNotFoundError(f"Could not find folder 'sample_{int(folder_name):06d}' in training_data.")
    if os.path.sep in folder_name:
        candidate = os.path.join(PROJECT_ROOT, folder_name)
        if os.path.isdir(candidate):
            return candidate
        raise FileNotFoundError(f"Folder not found: {folder_name}")
    # abstract search locations into a list in case we add more locations later
    candidates = []
    for search_dir in SEARCH_DIRS:
        candidate = os.path.join(PROJECT_ROOT, search_dir, folder_name)
        if os.path.isdir(candidate):
            candidates.append(candidate)
    if len(candidates) == 1:
        return candidates[0]
    elif len(candidates) > 1:
        print(f"Error: Multiple folders found for '{folder_name}':")
        for idx, c in enumerate(candidates):
            print(f"  [{idx+1}] {c}")
        print("Please specify the full path to the folder you want to visualize.")
        return None
    else:
        raise FileNotFoundError(f"Could not find folder '{folder_name}' in any of: {', '.join(SEARCH_DIRS)}.")


def run_visualizer(folder_name, initial_speed):
    data_dir = resolve_folder(folder_name)
    if data_dir is None:
        return
    x, frames, frame_labels = load_frames(data_dir)

    display_stride = max(1, len(x) // MAX_DISPLAY_POINTS)
    x_display = x[::display_stride]

    fig, ax = plt.subplots()
    fig.canvas.manager.set_window_title(f"Burgers Visualizer – {folder_name}")
    plt.subplots_adjust(bottom=0.25)
    line, = ax.plot(x_display, frames[0, ::display_stride])
    ax.set_title(f"{frame_labels[0]}   (frame 0)")

    frame_pos = 0.0
    frame_idx = 0
    playing = False
    speed = float(initial_speed)

    y_lock = True
    y_margin_frac = 0.05

    def update_plot():
        line.set_ydata(frames[frame_idx, ::display_stride])
        ax.set_title(f"{frame_labels[frame_idx]}   (frame {frame_idx})")

        if not y_lock:
            y_min = float(np.min(frames[frame_idx]))
            y_max = float(np.max(frames[frame_idx]))

            if y_max == y_min:
                delta = max(abs(y_min) * 0.1, 0.5)
                y_min -= delta
                y_max += delta

            margin = y_margin_frac * (y_max - y_min)
            ax.set_ylim(y_min - margin, y_max + margin)

        fig.canvas.draw_idle()

    ax_lock = plt.axes([0.82, 0.15, 0.12, 0.08])
    lock_cb = widgets.CheckButtons(ax_lock, ["Lock Y"], [y_lock])

    def on_toggle_lock(label):
        nonlocal y_lock
        y_lock = not y_lock

        if y_lock:
            current_ylim = ax.get_ylim()
            ax.set_ylim(current_ylim)
        else:
            update_plot()
        fig.canvas.draw_idle()
    lock_cb.on_clicked(on_toggle_lock)

    def next_frame(event):
        nonlocal frame_pos, frame_idx
        frame_idx = min(frame_idx + 1, len(frames) - 1)
        frame_pos = float(frame_idx)
        update_plot()

    def prev_frame(event):
        nonlocal frame_pos, frame_idx
        frame_idx = max(frame_idx - 1, 0)
        frame_pos = float(frame_idx)
        update_plot()

    def play_pause(event):
        nonlocal playing, frame_pos, frame_idx

        if frame_idx >= len(frames) - 1:
            frame_pos = 0.0
            frame_idx = 0
            update_plot()
        playing = not playing
        play_button.label.set_text("Pause" if playing else "Play")
        fig.canvas.draw_idle()

    def reset(event):
        nonlocal frame_pos, frame_idx, playing
        playing = False
        play_button.label.set_text("Play")
        frame_pos = 0.0
        frame_idx = 0
        update_plot()
        fig.canvas.draw_idle()

    def change_speed(val):
        nonlocal speed
        speed = float(val)

    axprev = plt.axes([0.1, 0.1, 0.1, 0.075])
    axplay = plt.axes([0.23, 0.1, 0.15, 0.075])
    axreset = plt.axes([0.40, 0.1, 0.12, 0.075])
    axnext = plt.axes([0.55, 0.1, 0.1, 0.075])
    axspeed = plt.axes([0.74, 0.1, 0.15, 0.04])

    prev_button = widgets.Button(axprev, "Prev")
    play_button = widgets.Button(axplay, "Play")
    reset_button = widgets.Button(axreset, "Reset")
    next_button = widgets.Button(axnext, "Next")
    speed_slider = widgets.Slider(axspeed, "Speed", 0.1, 100.0, valinit=speed)
    prev_button.on_clicked(prev_frame)
    next_button.on_clicked(next_frame)
    play_button.on_clicked(play_pause)
    reset_button.on_clicked(reset)
    speed_slider.on_changed(change_speed)

    timer = fig.canvas.new_timer(interval=30)

    def on_timer():
        nonlocal frame_pos, frame_idx, playing
        if not playing:
            return
        frame_pos += speed

        if frame_pos >= len(frames) - 1:
            frame_pos = len(frames) - 1
            frame_idx = int(frame_pos)
            update_plot()

            playing = False
            play_button.label.set_text("Play")
            fig.canvas.draw_idle()
            return

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
        help="Sample name (e.g. sample_000000), integer (e.g. 5), or folder path. Searches training_data/ and data/",
    )
    parser.add_argument(
        "--speed",
        type=float,
        default=1,
        help="Initial playback speed (frames per tick)",
    )
    args = parser.parse_args()
    try:
        run_visualizer(args.folder, args.speed)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        exit(1)


if __name__ == "__main__":
    main()
