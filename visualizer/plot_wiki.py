import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_ROOT, "data_v2")

def load_first_timestep(folder_name):
    """Loads timestep_00000.csv from the given folder and returns x, u arrays."""
    file_path = os.path.join(DATA_DIR, folder_name, "timestep_00000.csv")
    if not os.path.exists(file_path):
        print(f"Warning: File not found -> {file_path}")
        return None, None
    
    df = pd.read_csv(file_path)
    # Shift x-axis back from [0, 12] to [-6, 6] to match Wikipedia plots
    x = df['x'].values - 6.0
    u = df['u'].values
    return x, u

def plot_test_cases(scheme):
    viscosities = ["1.0", "0.1", "0.01"]
    colors = {"1.0": "black", "0.1": "red", "0.01": "blue"}
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
    fig.canvas.manager.set_window_title(f"Wikipedia Burgers' Equation Initial Conditions - {scheme}")

    # Plot Gaussian
    for v in viscosities:
        folder = f"wiki_gaussian_{v}_{scheme}"
        x, u = load_first_timestep(folder)
        if x is not None and u is not None:
            ax1.plot(x, u, label=f"v = {v}", color=colors[v])
            
    ax1.set_title("Gaussian Initial Condition: $u(x,0) = e^{-x^2/2}$")
    ax1.set_xlabel("x")
    ax1.set_ylabel("u")
    ax1.set_xlim([-7, 7])
    ax1.set_ylim([0, 1.1])
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot N-wave
    for v in viscosities:
        folder = f"wiki_n_wave_{v}_{scheme}"
        x, u = load_first_timestep(folder)
        if x is not None and u is not None:
            ax2.plot(x, u, label=f"v = {v}", color=colors[v])

    ax2.set_title("N-wave Initial Condition: $u(x,0) = e^{-(x-1)^2/2} - e^{-(x+1)^2/2}$")
    ax2.set_xlabel("x")
    ax2.set_ylabel("u")
    ax2.set_xlim([-7, 7])
    ax2.set_ylim([-1.1, 1.1])
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(SCRIPT_DIR, f"wikipedia_comparison_{scheme}.png"))
    print(f"Saved to wikipedia_comparison_{scheme}.png")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot Wikipedia Burgers' equation data")
    parser.add_argument(
        "--scheme",
        type=str,
        default="LaxWendroff",
        choices=["LaxWendroff", "Godunov"],
        help="Numerical scheme used for the simulation",
    )
    args = parser.parse_args()
    plot_test_cases(args.scheme)
