import argparse
import glob
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

SEARCH_DIRS = ["training_data", "data"]

# GLOBAL CONTROLS
# << and >> for playing all cells at the same speed (all cells will sync to the same time)
# < and > for stepping all cells one timestep based on the largest dt among all samples
# R for resetting all cells and global time 0
# LOCAL CONTROLS
# << and >> for playing the cell at the local speed
# < and > for stepping the cell one frame (same as one timestep based on its local dt)
# R for resetting the cell to time 0 (this does not affect global time so using the global controls will snap it back to global time)

def load_frames(data_dir):
    """
    Loads solution data from either:
    - Binary format: solution.bin + metadata.json
    - CSV format: timestep_*.csv files
    
    Returns:
        x: np.array of shape (N,)
        frames: list of np.array, each shape (N,)
        files: list of strings (file names or synthetic names)
        dt: float, time step size from metadata (or 1.0 for CSV fallback)
    """
    bin_path = os.path.join(data_dir, "solution.bin")
    meta_path = os.path.join(data_dir, "metadata.json")
    csv_files = sorted(glob.glob(os.path.join(data_dir, "timestep_*.csv")))
    
    dt = 1.0  # Default value
    
    # Try binary format first
    if os.path.exists(bin_path) and os.path.exists(meta_path):
        # Parse solution metadata (JSON)
        with open(meta_path, "r") as f:
            meta = json.load(f)

        num_domain_points = int(meta["solver"]["num_domain_points"])
        num_timesteps = int(meta["solver"]["time_steps"])
        dx = float(meta["solver"]["spatial_step_size"])
        dt = float(meta["solver"]["time_step_size"])

        # Load binary data
        data = np.fromfile(bin_path, dtype=np.float64)

        expected = num_domain_points * num_timesteps
        if data.size != expected:
            raise ValueError(
                f"Binary size mismatch for {bin_path}. "
                f"Expected {expected} float64 values (num_timesteps={num_timesteps}, num_domain_points={num_domain_points}), got {data.size}."
            )

        u = data.reshape((num_timesteps, num_domain_points))

        # Reconstruct x grid
        x = dx * np.arange(num_domain_points)

        # Return frames as a list of 1D arrays
        frames = [u[i, :] for i in range(num_timesteps)]

        # Return synthetic file names
        files = [f"timestep_{i:05d}.csv" for i in range(num_timesteps)]

        return x, frames, files, dt
    
    # Fall back to CSV format (this is legacy now but i've kept it for compatibility with the old format)
    elif csv_files:
        frames = []
        x = None

        for f in csv_files:
            data = np.loadtxt(f, delimiter=",", skiprows=1)
            if x is None:
                x = data[:, 0]
            frames.append(data[:, 1])

        return x, frames, csv_files, dt
    
    else:
        raise ValueError(f"No solution data found in {data_dir}. Expected either solution.bin+metadata.json or timestep_*.csv files.")

def resolve_folder(folder_name):
    """
    Resolves a folder path for a given sample name, integer, or relative path.
    Supports:
    - Integer (e.g., 5 → sample_000005 in training_data)
    - Folder name (e.g., sample_000000 → searches SEARCH_DIRS)
    - Relative path (e.g., training_data/sample_000000)
    
    Returns the resolved folder path or raises FileNotFoundError/ValueError.
    """
    # If folder_name is an integer string, treat as sample_XXXXXX in training_data
    if folder_name.isdigit():
        sample_name = f"sample_{int(folder_name):06d}"
        candidate = os.path.join(PROJECT_ROOT, "training_data", sample_name)
        if os.path.isdir(candidate):
            return candidate
        raise FileNotFoundError(f"Could not find folder 'sample_{int(folder_name):06d}' in training_data.")
    
    # If folder_name contains path separators, try a few things (I extended this because I was having some issues when I was in certain CWDs)
    # Try in order: absolute path, path relative to current working directory, path relative to PROJECT_ROOT
    if os.path.sep in folder_name:
        # Absolute path
        if os.path.isabs(folder_name):
            if os.path.isdir(folder_name):
                return folder_name
            raise FileNotFoundError(f"Folder not found: {folder_name}")

        # Path relative to current working directory
        cwd_candidate = os.path.abspath(folder_name)
        if os.path.isdir(cwd_candidate):
            return cwd_candidate

        # Also try relative to the script directory (useful when user runs from inside visualizer/)
        script_candidate = os.path.normpath(os.path.join(SCRIPT_DIR, folder_name))
        if os.path.isdir(script_candidate):
            return script_candidate

        # Path relative to project root (original behavior)
        candidate = os.path.join(PROJECT_ROOT, folder_name)
        if os.path.isdir(candidate):
            return candidate

        raise FileNotFoundError(f"Folder not found: {folder_name}")
    
    # Otherwise search in SEARCH_DIRS
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


class VisualizerCell:
    """
    Represents one complete visualizer window (plot + all controls).
    Each cell manages its own state independently.
    Controls are placed in a separate axis below the plot.
    """
    # Class variable to track locked y-axis state (global when multiple cells)
    y_locked_min = None
    y_locked_max = None
    all_cells = []  # Registry of all cells for coordinated updates
    
    def __init__(self, ax_plot, ax_controls, folder_path, initial_speed=1, preloaded_data=None, main_visualizer=None):
        self.ax = ax_plot
        self.ax_controls = ax_controls
        self.folder_path = folder_path
        self.folder_name = os.path.basename(folder_path)
        self.main_visualizer = main_visualizer  # Reference to MultiGridVisualizer
        
        # Load data (or use preloaded data)
        if preloaded_data is not None:
            self.x, self.frames, self.files, self.dt = preloaded_data
        else:
            self.x, self.frames, self.files, self.dt = load_frames(folder_path)
        
        # State variables
        self.time_pos = 0.0          # Current continuous time in seconds
        self.frame_idx = 0           # Current discrete frame index (round(time_pos / dt))
        self.playing = False
        self.speed = initial_speed
        self.y_lock = True
        self.y_margin_frac = 0.05
        
        # Global sync state: whether this cell is synchronized to global time
        self.synced_to_global = True
        
        # Create plot line
        self.line, = self.ax.plot(self.x, self.frames[0])
        # Set initial title with folder name (sample identifier) and frame/time info
        self._update_plot()
        
        # Register this cell
        VisualizerCell.all_cells.append(self)
        
        # Create control widgets in the controls axis
        self._create_controls()
        
        # Create and start timer for this cell
        self.timer = None
        self._start_timer()
    
    def _create_controls(self):
        """
        Create custom controls: < << >> R > (play forward/backward/reset) and pause block.
        Speed slider below buttons, and y-lock checkbox with global lock checkbox.
        """
        # Create button/slider axes within the control area using relative positioning
        fig = self.ax.figure
        pos = self.ax_controls.get_position()
        
        # Control area dimensions
        ctrl_left = pos.x0
        ctrl_bottom = pos.y0
        ctrl_width = pos.width
        ctrl_height = pos.height
        
        # Button height (40% of control area height)
        button_height = ctrl_height * 0.4
        
        # Convert button_height to figure coordinates to calculate actual width
        # Calculate the width needed to be square based on figure aspect ratio
        fig_width_inches = fig.get_figwidth()
        fig_height_inches = fig.get_figheight()
        aspect = fig_height_inches / fig_width_inches
        
        # Button width in figure coords to be square
        button_width_fig = button_height * aspect
        
        # Positions within the control area: buttons at top (55%-95%), slider below (35%-50%)
        button_y = ctrl_bottom + ctrl_height * 0.55
        slider_y = ctrl_bottom + ctrl_height * 0.35
        
        # Button row (top): 5 square buttons with consistent spacing
        button_spacing = button_width_fig * 0.1
        total_button_width = button_width_fig * 5 + button_spacing * 4
        buttons_left = ctrl_left + ctrl_width * 0.02
        
        ax_prev = fig.add_axes([buttons_left, button_y, button_width_fig, button_height])
        ax_rewind = fig.add_axes([buttons_left + (button_width_fig + button_spacing), button_y, button_width_fig, button_height])
        ax_play = fig.add_axes([buttons_left + (button_width_fig + button_spacing) * 2, button_y, button_width_fig, button_height])
        ax_next = fig.add_axes([buttons_left + (button_width_fig + button_spacing) * 3, button_y, button_width_fig, button_height])
        ax_reset = fig.add_axes([buttons_left + (button_width_fig + button_spacing) * 4, button_y, button_width_fig, button_height])
        
        # Speed slider (below buttons, same width as button group)
        slider_width = total_button_width
        ax_speed = fig.add_axes([buttons_left, slider_y, slider_width, button_height * 0.6])
        
        # Lock Y checkbox: anchor to right edge of control area, spaced away from playback buttons
        lock_width = ctrl_width * 0.15  # reasonable size
        # place it a small margin from right edge
        lock_x = ctrl_left + ctrl_width - lock_width - ctrl_width * 0.02
        ax_lock = fig.add_axes([lock_x, button_y, 
                                lock_width, button_height])
        
        self.prev_button = widgets.Button(ax_prev, "<")
        self.rewind_button = widgets.Button(ax_rewind, "<<")
        self.play_button = widgets.Button(ax_play, ">>")
        self.next_button = widgets.Button(ax_next, ">")
        self.reset_button = widgets.Button(ax_reset, "R")
        self.speed_slider = widgets.Slider(ax_speed, "", 1.0, 100.0, valinit=self.speed, valstep=0.25)
        # Adjust label position to avoid overlap with slider thumb at max
        self.speed_slider.label.set_x(0.2)
        self.lock_cb = widgets.CheckButtons(ax_lock, ["Lock Y"], [self.y_lock])
        
        # scale font size with button dimensions
        button_height_inches = button_height * fig_height_inches
        fontsize = max(8, button_height_inches * 14)  # scale factor 14, min 8pt
        
        self.prev_button.label.set_fontsize(fontsize)
        self.rewind_button.label.set_fontsize(fontsize)
        self.play_button.label.set_fontsize(fontsize)
        self.next_button.label.set_fontsize(fontsize)
        self.reset_button.label.set_fontsize(fontsize)
        self.speed_slider.label.set_fontsize(fontsize)
        
        # Track playback direction: 1 = forward, -1 = backward, 0 = stopped
        self.playback_direction = 0
        
        # Connect callbacks
        self.prev_button.on_clicked(self._prev_frame)
        self.rewind_button.on_clicked(self._play_backward)
        self.play_button.on_clicked(self._play_forward)
        self.next_button.on_clicked(self._next_frame)
        self.reset_button.on_clicked(self._on_reset)
        self.speed_slider.on_changed(self._change_speed)
        self.lock_cb.on_clicked(self._on_toggle_lock)
    
    def _update_plot(self):
        """Update the plotted line and optionally autoscale the y-axis."""
        # Calculate current frame index from time_pos
        max_frame = len(self.frames) - 1
        self.frame_idx = int(round(self.time_pos / self.dt))
        self.frame_idx = max(0, min(self.frame_idx, max_frame))
        
        self.line.set_ydata(self.frames[self.frame_idx])
        
        # Display time and frame (time is continuous, frame is discrete)
        name = self.folder_name
        title_text = (
            f"{name} | frame = {self.frame_idx}\n"
            f"time = {self.time_pos:.8f}"
        )
        self.ax.set_title(title_text)
        
        # Scale title fontsize based on figure dimensions
        fig = self.ax.figure
        fig_height_inches = fig.get_figheight()
        ax_bbox = self.ax.get_position()
        plot_height_inches = ax_bbox.height * fig_height_inches
        title_fontsize = max(10, plot_height_inches * 2.2)
        self.ax.title.set_fontsize(title_fontsize)
        
        # If lock is enabled, use stored lock limits
        if self.y_lock:
            if VisualizerCell.y_locked_min is not None and VisualizerCell.y_locked_max is not None:
                margin = self.y_margin_frac * (VisualizerCell.y_locked_max - VisualizerCell.y_locked_min)
                self.ax.set_ylim(VisualizerCell.y_locked_min - margin, VisualizerCell.y_locked_max + margin)
        # If lock is disabled, autoscale to current frame
        else:
            y_min = float(np.min(self.frames[self.frame_idx]))
            y_max = float(np.max(self.frames[self.frame_idx]))
            
            if y_max == y_min:
                delta = max(abs(y_min) * 0.1, 0.5)
                y_min -= delta
                y_max += delta
            
            margin = self.y_margin_frac * (y_max - y_min)
            self.ax.set_ylim(y_min - margin, y_max + margin)
        
        self.ax.figure.canvas.draw_idle()
    
    def _prev_frame(self, event):
        """Go to previous discrete frame."""
        self.desync_from_global()
        # Move back one frame relative to current frame index
        self.frame_idx = (self.frame_idx - 1) % len(self.frames)
        self.time_pos = float(self.frame_idx) * self.dt
        self.playing = False
        self.playback_direction = 0
        self._update_button_states()
        self._update_plot()
    
    def _next_frame(self, event):
        """Go to next discrete frame."""
        self.desync_from_global()
        # Move forward one frame relative to current frame index
        self.frame_idx = (self.frame_idx + 1) % len(self.frames)
        self.time_pos = float(self.frame_idx) * self.dt
        self.playing = False
        self.playback_direction = 0
        self._update_button_states()
        self._update_plot()
    
    def _play_forward(self, event):
        """Play forward (continuous time) or pause."""
        self.desync_from_global()
        if self.playing and self.playback_direction == 1:
            self.playing = False
            self.playback_direction = 0
        else:
            max_time = (len(self.frames) - 1) * self.dt
            if self.time_pos >= max_time - 1e-9:
                self.time_pos = 0.0
                self._update_plot()
            self.playing = True
            self.playback_direction = 1
        self._update_button_states()
    
    def _play_backward(self, event):
        """Play backward (continuous time) or pause."""
        self.desync_from_global()
        if self.playing and self.playback_direction == -1:
            self.playing = False
            self.playback_direction = 0
        else:
            max_time = (len(self.frames) - 1) * self.dt
            if self.time_pos <= 1e-9:
                self.time_pos = max_time
                self._update_plot()
            self.playing = True
            self.playback_direction = -1
        self._update_button_states()
    
    def _update_button_states(self):
        """Update button text to show pause state (■) when playing."""
        if self.playing:
            if self.playback_direction == 1:
                self.play_button.label.set_text("■")
                self.rewind_button.label.set_text("<<")
            elif self.playback_direction == -1:
                self.rewind_button.label.set_text("■")
                self.play_button.label.set_text(">>")
        else:
            self.play_button.label.set_text(">>")
            self.rewind_button.label.set_text("<<")
        self.ax.figure.canvas.draw_idle()
    
    def _change_speed(self, val):
        """Update speed from slider."""
        self.speed = float(val)
    
    def _on_reset(self, event):
        """Reset cell to initial state: time 0, stopped playback."""
        self.desync_from_global()
        self.reset()
    
    def _on_toggle_lock(self, label):
        """Toggle y-axis lock."""
        self.y_lock = not self.y_lock
        if self.y_lock:
            # When locking, calculate min/max across all frames
            # If multiple cells, use global lock; otherwise use just this cell's data
            if len(VisualizerCell.all_cells) > 1:
                # Multiple cells: calculate across all cells
                y_min = float('inf')
                y_max = float('-inf')
                for cell in VisualizerCell.all_cells:
                    for frame in cell.frames:
                        y_min = min(y_min, float(np.min(frame)))
                        y_max = max(y_max, float(np.max(frame)))
            else:
                # Single cell: calculate across this cell's frames
                y_min = float('inf')
                y_max = float('-inf')
                for frame in self.frames:
                    y_min = min(y_min, float(np.min(frame)))
                    y_max = max(y_max, float(np.max(frame)))
            
            if y_max == y_min:
                delta = max(abs(y_min) * 0.1, 0.5)
                y_min -= delta
                y_max += delta
            
            # Store in class variables so all cells can access
            VisualizerCell.y_locked_min = y_min
            VisualizerCell.y_locked_max = y_max
            
            # Update all cells if multiple exist
            if len(VisualizerCell.all_cells) > 1:
                for cell in VisualizerCell.all_cells:
                    cell._update_plot()
            else:
                self._update_plot()
        else:
            # When unlocking, autoscale to current frame
            for cell in VisualizerCell.all_cells:
                cell._update_plot()
        
        self.ax.figure.canvas.draw_idle()
    
    def _start_timer(self):
        """Create and start animation timer for this cell."""
        self.timer = self.ax.figure.canvas.new_timer(interval=30)
        self.timer.add_callback(self._on_timer)
        self.timer.start()
    
    def _on_timer(self):
        """Timer callback for local animation (continuous time)."""
        if not self.playing:
            return
        
        # Advance time by dt * speed (one frame's worth of time per tick)
        self.time_pos += self.dt * self.speed * self.playback_direction
        
        max_time = (len(self.frames) - 1) * self.dt
        
        # Boundary handling for Play Mode (stop at ends)
        if self.playback_direction == 1:
            if self.time_pos >= max_time - 1e-9:
                self.time_pos = max_time
                self._update_plot()
                self.playing = False
                self.playback_direction = 0
                self._update_button_states()
                return
        elif self.playback_direction == -1:
            if self.time_pos <= 1e-9:
                self.time_pos = 0.0
                self._update_plot()
                self.playing = False
                self.playback_direction = 0
                self._update_button_states()
                return
        
        self._update_plot()
    
    def update_from_global_time(self, global_time, global_playing):
        """Update based on global seconds."""
        if global_playing and self.synced_to_global:
            max_time = (len(self.frames) - 1) * self.dt
            self.time_pos = max(0.0, min(global_time, max_time))
            self._update_plot()
    
    def snap_to_global_time(self, global_time):
        """Snap to global seconds."""
        max_time = (len(self.frames) - 1) * self.dt
        self.time_pos = max(0.0, min(global_time, max_time))
        self.synced_to_global = True
        self.playing = False
        self.playback_direction = 0
        self._update_button_states()
        self._update_plot()
    
    def desync_from_global(self):
        """Mark this cell as desynced from global time (user interacted with local controls)."""
        self.synced_to_global = False
        
        # Stop global playback if it is currently active
        if self.main_visualizer and self.main_visualizer.global_playing:
            self.main_visualizer.global_playing = False
            self.main_visualizer.global_playback_direction = 0
            self.main_visualizer._update_global_button_states()

    def reset(self):
        """Reset cell to initial state: time 0, stopped playback."""
        self.time_pos = 0.0
        self.frame_idx = 0
        self.playing = False
        self.playback_direction = 0
        self._update_button_states()
        self._update_plot()

    def stop(self):
        """Stop the timer."""
        if self.timer:
            self.timer.stop()


class MultiGridVisualizer:
    """
    Manages an N×M grid of VisualizerCell objects.
    Uses GridSpec to allocate space for plot and controls per cell.
    Includes optional global play controls when multiple samples are present.
    """
    def __init__(self, rows, cols, folder_paths, initial_speed=1):
        self.rows = rows
        self.cols = cols
        self.folder_paths = folder_paths
        self.initial_speed = initial_speed
        self.cells = []
        
        # Global control state (only used for multiple samples)
        self.global_time = 0.0
        self.global_speed = initial_speed
        self.global_playing = False
        self.global_playback_direction = 0
        self.global_timer = None
        
        # Validate that we have enough folders
        if len(folder_paths) > rows * cols:
            raise ValueError(f"Too many folders ({len(folder_paths)}) for grid size {rows}×{cols}")
        
        # Create figure
        # always reserve a thin bar at the top for global controls
        # reserve a little extra vertical space for the global bar
        self.fig = plt.figure(figsize=(6*cols, 7*rows + 1.5))
        
        # Create GridSpec with an extra row at top for the global bar
        # Each plot row still has a control row beneath it
        gs = plt.GridSpec(rows * 2 + 1, cols, figure=self.fig,
                         height_ratios=[1] + [4, 1] * rows,
                         hspace=0.10, wspace=0.15)
        global_gs = gs[0, :]
        

        
        self.fig.canvas.manager.set_window_title(f"Burgers Visualizer – Grid {rows}×{cols}")
        
        # Reduce outer margins (left, right, top, bottom)
        # increase top margin to give some space above global controls
        self.fig.subplots_adjust(left=0.1, right=0.95, top=0.98, bottom=0.05)
        
        # Load all frames sequentially
        preloaded_data_list = [load_frames(folder_path) for folder_path in folder_paths]
        
        # Create a VisualizerCell for each folder with preloaded data
        for idx, (folder_path, preloaded_data) in enumerate(zip(folder_paths, preloaded_data_list)):
            row = idx // cols
            col = idx % cols
            
            # Plot axis (uses 4 parts of height; always offset by 1 due to reserved bar row)
            gs_row = row * 2 + 1
            ax_plot = self.fig.add_subplot(gs[gs_row, col])
            
            # Control axis (uses 1 part of height)
            ax_controls = self.fig.add_subplot(gs[gs_row + 1, col])
            ax_controls.axis('off')  # Hide the control axis itself
            
            cell = VisualizerCell(ax_plot, ax_controls, folder_path, initial_speed, preloaded_data, main_visualizer=self)
            self.cells.append(cell)
        
        # Hide unused subplots
        for idx in range(len(folder_paths), rows * cols):
            row = idx // cols
            col = idx % cols
            gs_row = row * 2 + 1
            ax_plot = self.fig.add_subplot(gs[gs_row, col])
            ax_controls = self.fig.add_subplot(gs[gs_row + 1, col])
            ax_plot.set_visible(False)
            ax_controls.set_visible(False)
        
        # Create global controls only when we have more than one sample
        if len(folder_paths) > 1:
            self._create_global_controls(global_gs)
        # (single-sample layout already reserves the bar; nothing else needed)
        
        # Initialize locked y-axis limits (since y_lock is True by default)
        self._initialize_locked_limits()
    
    def _create_global_controls(self, gs_global):
        """Create global play controls matching local button dimensions.
        Buttons arranged in bar and slider below like local controls."""
        fig = self.fig
        ax_global = fig.add_subplot(gs_global)
        ax_global.axis('off')

        # derive all dimensions from the bar's actual bbox in figure coordinates
        # so buttons/slider are ALWAYS the same physical size as the local control rows
        bb = ax_global.get_position()
        bar_left_fig = bb.x0
        bar_bottom_fig = bb.y0
        bar_width_fig = bb.width
        bar_height_fig = bb.height          # same ratio as local control rows

        # split bar: top 60% buttons, bottom 35% slider (with 5% padding)
        pad = bar_height_fig * 0.05
        button_height = bar_height_fig * 0.55
        button_bottom = bar_bottom_fig + bar_height_fig * 0.40   # sits in upper part
        slider_height = bar_height_fig * 0.30
        slider_bottom = bar_bottom_fig + pad                     # sits at very bottom

        # square buttons, horizontally centered
        button_width = button_height
        button_spacing = button_width * 0.1
        total_width = button_width * 5 + button_spacing * 4
        centre_x = bar_left_fig + bar_width_fig / 2
        buttons_left = centre_x - total_width / 2

        ax_prev   = fig.add_axes([buttons_left,                                  button_bottom, button_width, button_height])
        ax_rewind = fig.add_axes([buttons_left + (button_width + button_spacing), button_bottom, button_width, button_height])
        ax_play   = fig.add_axes([buttons_left + (button_width + button_spacing)*2, button_bottom, button_width, button_height])
        ax_next   = fig.add_axes([buttons_left + (button_width + button_spacing)*3, button_bottom, button_width, button_height])
        ax_reset  = fig.add_axes([buttons_left + (button_width + button_spacing)*4, button_bottom, button_width, button_height])

        ax_speed  = fig.add_axes([buttons_left, slider_bottom, total_width, slider_height])

        self.global_prev_button = widgets.Button(ax_prev, "<")
        self.global_rewind_button = widgets.Button(ax_rewind, "<<")
        self.global_play_button = widgets.Button(ax_play, ">>")
        self.global_next_button = widgets.Button(ax_next, ">")
        self.global_reset_button = widgets.Button(ax_reset, "R")
        self.global_speed_slider = widgets.Slider(ax_speed, "", 1.0, 100.0, valinit=self.global_speed, valstep=0.25)
        self.global_speed_slider.label.set_x(0.1)

        # scale font size with button dimensions
        fig_width_inches, fig_height_inches = fig.get_size_inches()
        button_height_inches = button_height * fig_height_inches
        fontsize = max(8, button_height_inches * 14)  # scale factor 14, min 8pt
        
        self.global_prev_button.label.set_fontsize(fontsize)
        self.global_rewind_button.label.set_fontsize(fontsize)
        self.global_play_button.label.set_fontsize(fontsize)
        self.global_next_button.label.set_fontsize(fontsize)
        self.global_reset_button.label.set_fontsize(fontsize)
        self.global_speed_slider.label.set_fontsize(fontsize)

        # callbacks
        self.global_prev_button.on_clicked(self._global_prev_frame)
        self.global_rewind_button.on_clicked(self._global_play_backward)
        self.global_play_button.on_clicked(self._global_play_forward)
        self.global_next_button.on_clicked(self._global_next_frame)
        self.global_reset_button.on_clicked(self._global_reset)
        self.global_speed_slider.on_changed(self._global_change_speed)

        # start global timer
        self.global_timer = fig.canvas.new_timer(interval=30)
        self.global_timer.add_callback(self._global_on_timer)
        self.global_timer.start()
    

    
    def _global_prev_frame(self, event):
        """Go to previous seconds step in global time (1 max_dt back)."""
        max_valid_time = max((len(cell.frames) - 1) * cell.dt for cell in self.cells)
        max_dt = max(cell.dt for cell in self.cells)
        
        next_time = self.global_time - max_dt
        # If we would cross 0, loop to end
        self.global_time = max_valid_time if next_time < -1e-9 else max(0.0, next_time)
        
        for cell in self.cells:
            cell.snap_to_global_time(self.global_time)
        self._update_global_button_states()
    
    def _global_next_frame(self, event):
        """Go to next seconds step in global time (1 max_dt forward)."""
        max_valid_time = max((len(cell.frames) - 1) * cell.dt for cell in self.cells)
        max_dt = max(cell.dt for cell in self.cells)
        
        next_time = self.global_time + max_dt
        # If we would go past end, loop to start
        self.global_time = 0.0 if next_time > max_valid_time + 1e-9 else min(max_valid_time, next_time)
        
        for cell in self.cells:
            cell.snap_to_global_time(self.global_time)
        self._update_global_button_states()
    
    def _global_reset(self, event):
        """Reset global time to 0 and stop all playback."""
        self.global_time = 0.0
        self.global_playing = False
        self.global_playback_direction = 0
        
        # Reset all cells
        for cell in self.cells:
            cell.reset()
        
        self._update_global_button_states()
    
    def _global_play_forward(self, event):
        """Play forward (continuous time) or pause."""
        if self.global_playing and self.global_playback_direction == 1:
            self.global_playing = False
            self.global_playback_direction = 0
        else:
            max_valid_time = max((len(cell.frames) - 1) * cell.dt for cell in self.cells)
            if self.global_time >= max_valid_time - 1e-9:
                self.global_time = 0.0
            
            self.global_playing = True
            self.global_playback_direction = 1
            for cell in self.cells:
                cell.snap_to_global_time(self.global_time)
        self._update_global_button_states()
    
    def _global_play_backward(self, event):
        """Play backward (continuous time) or pause."""
        if self.global_playing and self.global_playback_direction == -1:
            self.global_playing = False
            self.global_playback_direction = 0
        else:
            max_valid_time = max((len(cell.frames) - 1) * cell.dt for cell in self.cells)
            if self.global_time <= 1e-9:
                self.global_time = max_valid_time
            
            self.global_playing = True
            self.global_playback_direction = -1
            for cell in self.cells:
                cell.snap_to_global_time(self.global_time)
        self._update_global_button_states()
    
    def _global_change_speed(self, val):
        """Update global speed from slider."""
        self.global_speed = float(val)

    def _update_global_button_states(self):
        """Update global button text to show pause state."""
        if self.global_playing:
            if self.global_playback_direction == 1:
                self.global_play_button.label.set_text("■")
                self.global_rewind_button.label.set_text("<<")
            elif self.global_playback_direction == -1:
                self.global_rewind_button.label.set_text("■")
                self.global_play_button.label.set_text(">>")
        else:
            self.global_play_button.label.set_text(">>")
            self.global_rewind_button.label.set_text("<<")
        self.fig.canvas.draw_idle()
    
    def _global_on_timer(self):
        """Global timer callback for coordinated playback (seconds base)."""
        if not self.global_playing:
            return
        
        max_valid_time = max((len(cell.frames) - 1) * cell.dt for cell in self.cells)
        max_dt = max(cell.dt for cell in self.cells)
        
        # Advance global time in seconds
        delta = self.global_speed * self.global_playback_direction * max_dt
        next_time = self.global_time + delta
        
        # Stop at boundaries
        if self.global_playback_direction == 1 and next_time >= max_valid_time:
            self.global_time = max_valid_time
            for cell in self.cells:
                cell.update_from_global_time(self.global_time, True)
            self.global_playing = False
            self.global_playback_direction = 0
            self._update_global_button_states()
            return
        elif self.global_playback_direction == -1 and next_time <= 0:
            self.global_time = 0.0
            for cell in self.cells:
                cell.update_from_global_time(self.global_time, True)
            self.global_playing = False
            self.global_playback_direction = 0
            self._update_global_button_states()
            return
        
        self.global_time = next_time
        for cell in self.cells:
            cell.update_from_global_time(self.global_time, self.global_playing)
    
    def _initialize_locked_limits(self):
        """Calculate and store the locked y-axis limits based on all cells."""
        y_min = float('inf')
        y_max = float('-inf')
        
        # Calculate across all cells
        for cell in VisualizerCell.all_cells:
            for frame in cell.frames:
                y_min = min(y_min, float(np.min(frame)))
                y_max = max(y_max, float(np.max(frame)))
        
        if y_max == y_min:
            delta = max(abs(y_min) * 0.1, 0.5)
            y_min -= delta
            y_max += delta
        
        VisualizerCell.y_locked_min = y_min
        VisualizerCell.y_locked_max = y_max
        
        # Update all cells to apply the locked limits
        for cell in VisualizerCell.all_cells:
            cell._update_plot()

    def show(self):
        """Display the grid."""
        plt.show()
    
    def cleanup(self):
        """Stop all timers."""
        for cell in self.cells:
            cell.stop()
        if self.global_timer:
            self.global_timer.stop()



def main():
    parser = argparse.ArgumentParser(
        description="Visualize Burgers' equation data in a grid",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single visualizer from integer ID
  python visualize.py 0 --speed 1
  
  # Single visualizer from folder name
  python visualize.py sample_000000 --speed 1
  
  # 2x2 grid with 4 samples using integer IDs
  python visualize.py -r 2 -c 2 0 1 2 3
  
  # 1x2 grid with 2 samples using full paths
  python visualize.py --rows 1 --cols 2 ../training_data/sample_000000 ../training_data/sample_000001
        """
    )
    
    parser.add_argument(
        "folders",
        type=str,
        nargs="+",
        help="Sample specifications: integer ID (e.g., 5→sample_000005), folder name (e.g., sample_000000), or path (e.g., training_data/sample_000000)",
    )
    parser.add_argument(
        "-r", "--rows",
        type=int,
        help="Number of rows in grid (optional; computed if omitted)",
    )
    parser.add_argument(
        "-c", "--cols",
        type=int,
        help="Number of columns in grid (optional; computed if omitted)",
    )
    parser.add_argument(
        "--speed",
        type=float,
        default=1,
        help="Initial playback speed (frames per tick, default: 1)",
    )

    args = parser.parse_args()
    
    # Resolve all folder specifications to actual paths
    resolved_folders = []
    for folder_spec in args.folders:
        try:
            resolved = resolve_folder(folder_spec)
            if resolved is not None:
                resolved_folders.append(resolved)
            else:
                print(f"Skipping ambiguous folder: {folder_spec}")
        except FileNotFoundError as e:
            print(f"Error: {e}")
            return
    
    if not resolved_folders:
        print("Error: No valid folders resolved.")
        return
    
    # derive rows/cols if not provided
    num = len(resolved_folders)
    rows = args.rows
    cols = args.cols
    if rows is None or cols is None:
        # choose grid as square as possible, bias wider: cols >= rows
        import math
        cols = math.ceil(math.sqrt(num))
        rows = math.ceil(num / cols)
    
    # Create and display grid visualizer
    visualizer = MultiGridVisualizer(rows, cols, resolved_folders, args.speed)
    visualizer.show()


if __name__ == "__main__":
    main()
