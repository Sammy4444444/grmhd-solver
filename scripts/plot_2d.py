

import sys
import os
import csv
import numpy as np

def plot_2d_file(csv_file="output_2d.csv"):
    if not os.path.exists(csv_file):
        # Fallback to output.csv if output_2d.csv doesn't exist
        if os.path.exists("output.csv"):
            csv_file = "output.csv"
        else:
            print(f"[ERROR] CSV file '{csv_file}' not found.")
            return False

    print(f"[INFO] Reading 2D simulation data from {csv_file}...")

    xs = []
    ys = []
    rhos = []
    presss = []
    alphas = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

        # If 'y' is missing from header, this is a 1D dataset
        if 'y' not in fieldnames:
            print(f"[INFO] File '{csv_file}' does not contain 'y' column (1D dataset).")
            return False

        for row in reader:
            xs.append(float(row.get('x', 0.0)))
            ys.append(float(row.get('y', 0.0)))
            rhos.append(float(row.get('rho', 0.0)))
            presss.append(float(row.get('press', 0.0)))
            alphas.append(float(row.get('alpha', 1.0)))

    if len(xs) == 0:
        print("[ERROR] CSV file is empty.")
        return False

    unique_x = np.unique(xs)
    unique_y = np.unique(ys)
    nx = len(unique_x)
    ny = len(unique_y)

    if nx * ny != len(xs):
        print(f"[WARNING] Data shape mismatch ({len(xs)} points vs {nx}x{ny} grid). Attempting reshape...")

    X = np.reshape(xs, (nx, ny))
    Y = np.reshape(ys, (nx, ny))
    Rho = np.reshape(rhos, (nx, ny))
    Press = np.reshape(presss, (nx, ny))
    Alpha = np.reshape(alphas, (nx, ny))

    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle

        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        fig.suptitle('2D GRMHD Accretion Torus in Kerr-Schild Black Hole Spacetime', fontsize=14, fontweight='bold')

        # 1. Density Rho contour map
        c1 = axes[0].pcolormesh(X, Y, Rho, cmap='inferno', shading='auto')
        fig.colorbar(c1, ax=axes[0], label=r'Rest Mass Density $\rho$')
        axes[0].set_title(r'Density $\rho(x, y)$')
        axes[0].set_xlabel('x [GM/c²]')
        axes[0].set_ylabel('y [GM/c²]')
        axes[0].set_aspect('equal')
        
        # Draw Event Horizon Circle (r ~ 1.9 for M=1, a=0.9375)
        bh_horizon = Circle((0, 0), 1.9, color='cyan', fill=False, linestyle='--', linewidth=1.5, label='Horizon r~1.9')
        axes[0].add_patch(bh_horizon)
        axes[0].legend(loc='upper right')

        # 2. Pressure P contour map
        c2 = axes[1].pcolormesh(X, Y, Press, cmap='magma', shading='auto')
        fig.colorbar(c2, ax=axes[1], label='Fluid Pressure P')
        axes[1].set_title(r'Pressure $P(x, y)$')
        axes[1].set_xlabel('x [GM/c²]')
        axes[1].set_ylabel('y [GM/c²]')
        axes[1].set_aspect('equal')

        # 3. Kerr Lapse Function alpha contour map
        c3 = axes[2].pcolormesh(X, Y, Alpha, cmap='viridis', shading='auto')
        fig.colorbar(c3, ax=axes[2], label=r'Lapse Function $\alpha$')
        axes[2].set_title(r'Kerr Lapse $\alpha(x, y)$')
        axes[2].set_xlabel('x [GM/c²]')
        axes[2].set_ylabel('y [GM/c²]')
        axes[2].set_aspect('equal')

        plt.tight_layout()
        output_img = "grmhd_2d_contour.png"
        plt.savefig(output_img, dpi=300)
        print(f"[SUCCESS] 2D Contour Plot saved to '{output_img}'")
        return True

    except ImportError:
        print("[WARNING] matplotlib is not installed. Data parsed successfully:")
        print(f"Grid size: {nx} x {ny}")
        print(f"x range: [{min(xs)}, {max(xs)}]")
        print(f"y range: [{min(ys)}, {max(ys)}]")
        return True

def main():
    target_file = "output_2d.csv"
    if len(sys.argv) > 1:
        target_file = sys.argv[1]
    plot_2d_file(target_file)

if __name__ == "__main__":
    main()
