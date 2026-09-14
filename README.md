# GRMHD Solver: General Relativistic Magnetohydrodynamics Engine from Scratch

![C++17](https://img.shields.io/badge/Language-C%2B%2B17-blue.svg)
![Julia](https://img.shields.io/badge/Prototyping-Julia-purple.svg)
![Python](https://img.shields.io/badge/Visualization-Python3-green.svg)
![CMake](https://img.shields.io/badge/Build-CMake%203.14%2B-orange.svg)
![Platform](https://img.shields.io/badge/OS-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)

A high-performance **2D General Relativistic Magnetohydrodynamics (GRMHD)** simulation engine written from scratch. This multi-language software suite combines a high-speed C++ computational core, Julia mathematical tensor prototypes, and Python high-resolution visualization routines to simulate plasma dynamics and accretion disks in curved spacetime around spinning Kerr black holes.

---

##  Physics & Numerical Methods

### 1. Spacetime Geometry & ADM 3+1 Formalism
The simulation operates on a **Kerr-Schild Spacetime** background surrounding a fast-spinning black hole with mass $M = 1.0$ and dimensionless spin parameter $a = 0.9375$. In Cartesian-like Kerr-Schild coordinates $(t, x, y, z)$, the metric decomposes under the ADM 3+1 split:

$$ds^2 = -\alpha^2 dt^2 + \gamma_{ij} (dx^i + \beta^i dt)(dx^j + \beta^j dt)$$

- **Lapse Function**: $\alpha = \frac{1}{\sqrt{1 + f}}$
- **Shift Vector**: $\beta^i = \frac{f}{1 + f} l_i$
- **Spatial Metric**: $\gamma_{ij} = \delta_{ij} + f l_i l_j$
- **4-Metric Determinant**: $\sqrt{-g} = \alpha \sqrt{\gamma} = 1.0$ (everywhere in Kerr-Schild)

where $f = \frac{2 M r^3}{r^4 + a^2 z^2}$ and $l_\mu = (1, \frac{r x + a y}{r^2 + a^2}, \frac{r y - a x}{r^2 + a^2}, \frac{z}{r})$ is the null vector.

### 2. Relativistic Conservation Laws & Conserved Variables
The GRMHD state vector contains primitive variables $P = (\rho, P, v^i, B^i)$ mapped to conserved variables $U = (D, S_i, \tau, B^i)$:

- **Conserved Density**: $D = \alpha \sqrt{\gamma} W \rho = W \rho$
- **Conserved Momentum**: $S_i = (\rho h_{\text{tot}} W^2 v_i - b_0 b_i) \sqrt{-g}$
- **Conserved Energy**: $\tau = (\rho h_{\text{tot}} W^2 - P_{\text{tot}} - b_0^2 - D) \sqrt{-g}$

### 3. 2D Spatial Reconstruction & HLL Riemann Solver
- **2D TVD Spatial Reconstruction**: High-order spatial accuracy is achieved via linear piecewise reconstruction using the **Monotonized Central (MC) Slope Limiter**:
  $$\text{MC}(a, b) = \text{minmod}\left(2a, \frac{a+b}{2}, 2b\right)$$
- **HLL Riemann Solver**: Numerical fluxes across cell interfaces ($X$ and $Y$ directions) are evaluated using local fast-magnetosonic wave speeds transformed into the 3+1 ADM Eulerian frame:
  $$\lambda_{\pm, \text{coord}} = \alpha \lambda_{\pm} - \beta^k$$

---

## Directory Structure

```text
grmhd-solver/
├── src/                      # Core C++17 Computational Engine
│   ├── main.cpp              # 2D Simulation loop & Kerr black hole setup
│   ├── grid.hpp              # 2D Grid mesh layout & Primitive/Conserved conversions
│   ├── metric.hpp            # Kerr-Schild 4D Metric & 3+1 ADM decomposition
│   ├── hll_solver.hpp        # 2D HLL Riemann Solver (X & Y direction fluxes)
│   └── reconstruction.hpp    # 2D TVD Monotonized Central (MC) spatial reconstruction
├── julia/                    # Mathematical & Tensor Prototypes
│   └── test_metric.jl        # Julia prototype script for Kerr metric calculations
├── scripts/                  # Data Processing & Visualization Scripts
│   ├── plot_2d.py            # 2D Heatmap & Contour rendering engine
│   └── plot_results.py       # Auto-detecting multi-dimensional plotter runner
├── CMakeLists.txt            # Modern CMake build system file
├── run.sh                    # Bash automation script for Linux/macOS/Git Bash
├── run.ps1                   # PowerShell automation script for Windows
└── README.md                 # Project documentation
```

---

##  Prerequisites

Ensure you have the following installed on your system:
- **C++ Compiler**: C++17 compatible (GCC 8+, Clang 7+, or MSVC 2019+)
- **Build System**: CMake 3.14 or newer
- **Python**: Python 3.8+ with `matplotlib` and `numpy`
- **Julia** *(Optional)*: Julia 1.6+ for matrix & tensor prototyping

---

##  Quick Start Guide

### 1. One-Click Automated Build & Execution

**Windows (PowerShell):**
```powershell
.\run.ps1
```

**Linux / macOS / Git Bash:**
```bash
chmod +x run.sh
./run.sh
```

### 2. Manual Compilation & Execution (CMake)

```bash
# Create build directory and compile
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release

# Execute 2D Solver
./build/Release/grmhd_solver.exe    # Windows
./build/grmhd_solver                # Linux/macOS

# Render 2D Contour Visualization
python scripts/plot_results.py output_2d.csv
```

---

##  Simulation Output & Visualizations

After completing execution, the solver generates `output_2d.csv` and renders a 3-panel visualization image `grmhd_2d_contour.png`:

1. **Density Map $\rho(x, y)$**: Shows the magnetized plasma accretion torus orbiting around the central black hole, with the Event Horizon marked at $r \approx 1.9 \ GM/c^2$.
2. **Fluid Pressure Map $P(x, y)$**: Highlights pressure gradients and shock structures within the accretion flow.
3. **Kerr Lapse Map $\alpha(x, y)$**: Demonstrates gravitational time dilation slowing down coordinate time near the black hole horizon ($\alpha \to 0$).

---

## License
This project is open-source and available under the MIT License.
