import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
SCHEME = "Godunov" # match -s argument used with running solver

def load_timestep(folder_name, timestep):
    file_path = os.path.join(DATA_DIR, folder_name, f"timestep_{timestep:05d}.csv")
    if not os.path.exists(file_path):
        return None, None
    df = pd.read_csv(file_path)
    x = df['x'].values - 6.0
    u = df['u'].values
    return x, u

def create_animation():
    viscosities = ["1.0", "0.1", "0.01"]
    colors = {"1.0": "black", "0.1": "red", "0.01": "blue"}
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    lines = {}
    for v in viscosities:
        lines[v], = ax.plot([], [], label=f"$\\nu = {v}$", color=colors[v])
        
    time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, fontsize=18,
                        verticalalignment='top')
    
    ax.set_xlim([-7, 7])
    ax.set_ylim([-1.1, 1.1])
    ax.set_xlabel("x", fontsize=16)
    ax.set_ylabel("u", fontsize=16)
    ax.legend(loc="lower right", fontsize=14)
    ax.set_title("N-wave Initial Condition Evolution", fontsize=16)
    
    # Check max frames
    folder = f"wiki_n_wave_1.0_{SCHEME}"
    max_frames = 0
    # Steps are saved every 10 steps. max time_steps is 5000, so max file is timestep_05000.csv
    # This means up to 501 frames.
    while os.path.exists(os.path.join(DATA_DIR, folder, f"timestep_{max_frames*10:05d}.csv")):
        max_frames += 1

    def init():
        for line in lines.values():
            line.set_data([], [])
        time_text.set_text('')
        return list(lines.values()) + [time_text]

    def animate(i):
        timestep = i * 10
        t = timestep * 0.001
        for v in viscosities:
            folder = f"wiki_n_wave_{v}_{SCHEME}"
            x, u = load_timestep(folder, timestep)
            if x is not None and u is not None:
                lines[v].set_data(x, u)
        
        time_text.set_text(f'$t = {t:.2f}$')
        return list(lines.values()) + [time_text]

    if max_frames == 0:
        print("No data found! Check if the C++ simulation successfully ran.")
        return

    print(f"Generating animation with {max_frames} frames... This might take a minute.")
    ani = animation.FuncAnimation(fig, animate, init_func=init,
                                  frames=max_frames, interval=20, blit=True)
    
    # Save as gif
    output_path = os.path.join(SCRIPT_DIR, "n_wave_animation.gif")
    writer = PillowWriter(fps=30)
    ani.save(output_path, writer=writer)
    print(f"Saved animation to {output_path}")

if __name__ == "__main__":
    create_animation()
