


import sys
import os
import csv
import plot_2d

def plot_1d_file(csv_file):
    print(f"[INFO] Rendering 1D plot for {csv_file}...")
    x = []
    rho = []
    press = []
    D = []
    tau = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            x.append(float(row.get('x', 0.0)))
            rho.append(float(row.get('rho', 0.0)))
            press.append(float(row.get('press', 0.0)))
            D.append(float(row.get('D', 0.0)))
            tau.append(float(row.get('tau', 0.0)))

    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        fig.suptitle('GRMHD Solver Simulation Output', fontsize=14, fontweight='bold')

        axes[0, 0].plot(x, rho, 'b-', label=r'$\rho$ (Density)', linewidth=2)
        axes[0, 0].set_ylabel(r'$\rho$')
        axes[0, 0].grid(True, linestyle='--', alpha=0.6)
        axes[0, 0].legend()

        axes[0, 1].plot(x, press, 'r-', label='Pressure (P)', linewidth=2)
        axes[0, 1].set_ylabel('P')
        axes[0, 1].grid(True, linestyle='--', alpha=0.6)
        axes[0, 1].legend()

        axes[1, 0].plot(x, D, 'g-', label='D (Conserved density)', linewidth=2)
        axes[1, 0].set_xlabel('x')
        axes[1, 0].set_ylabel('D')
        axes[1, 0].grid(True, linestyle='--', alpha=0.6)
        axes[1, 0].legend()

        axes[1, 1].plot(x, tau, 'm-', label=r'$\tau$ (Conserved energy)', linewidth=2)
        axes[1, 1].set_xlabel('x')
        axes[1, 1].set_ylabel(r'$\tau$')
        axes[1, 1].grid(True, linestyle='--', alpha=0.6)
        axes[1, 1].legend()

        plt.tight_layout()
        output_img = "grmhd_output.png"
        plt.savefig(output_img, dpi=300)
        print(f"[SUCCESS] 1D Plot saved to '{output_img}'")

    except ImportError:
        print("[WARNING] matplotlib is not installed.")

def main():
    target_file = "output_2d.csv"

    if len(sys.argv) > 1:
        target_file = sys.argv[1]

    if not os.path.exists(target_file):
        if os.path.exists("output_2d.csv"):
            target_file = "output_2d.csv"
        elif os.path.exists("output.csv"):
            target_file = "output.csv"
        else:
            print(f"[ERROR] Neither '{target_file}', 'output_2d.csv' nor 'output.csv' found.")
            sys.exit(1)

    # Read header to determine if data is 1D or 2D
    with open(target_file, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

    if 'y' in fieldnames:
        # 2D dataset detected
        success = plot_2d.plot_2d_file(target_file)
        if not success:
            plot_1d_file(target_file)
    else:
        # 1D dataset detected
        plot_1d_file(target_file)

if __name__ == "__main__":
    main()
