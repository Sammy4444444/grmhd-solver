#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "metric.hpp"
#include "grid.hpp"
#include "hll_solver.hpp"
#include "reconstruction.hpp"

using namespace GRMHD;

int main() {
    std::cout << "=========================================================\n";
    std::cout << " GRMHD Engine: 2D Full Grid Kerr-Schild Black Hole Solver\n";
    std::cout << "=========================================================\n";

    // 1. Black Hole Parameters
    double M = 1.0;
    double a = 0.9375; // Spinning Kerr black hole

    // 2. 2D Computational Grid Setup (X-Y plane: [-10.0, 10.0] x [-10.0, 10.0] GM/c^2)
    int nx = 60;
    int ny = 60;
    double x_min = -10.0, x_max = 10.0;
    double y_min = -10.0, y_max = 10.0;
    Grid grid(nx, ny, x_min, x_max, y_min, y_max);

    std::cout << "[INFO] Black Hole Mass M = " << M << ", Spin a = " << a << "\n";
    std::cout << "[INFO] 2D Grid initialized (" << nx << " x " << ny << " cells). dx = " << grid.dx << ", dy = " << grid.dy << "\n";

    // 3. Precompute 2D Kerr-Schild Metric at Cell Centers & Interfaces
    std::vector<std::vector<Metric>> cell_metrics(nx, std::vector<Metric>(ny));
    std::vector<std::vector<Metric>> interface_metrics_x(nx + 1, std::vector<Metric>(ny));
    std::vector<std::vector<Metric>> interface_metrics_y(nx, std::vector<Metric>(ny + 1));

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            cell_metrics[i][j].set_kerr_schild(grid.x[i], grid.y[j], 0.0, M, a);
        }
    }

    for (int i = 0; i <= nx; ++i) {
        double x_face = x_min + i * grid.dx;
        for (int j = 0; j < ny; ++j) {
            interface_metrics_x[i][j].set_kerr_schild(x_face, grid.y[j], 0.0, M, a);
        }
    }

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            double y_face = y_min + j * grid.dy;
            interface_metrics_y[i][j].set_kerr_schild(grid.x[i], y_face, 0.0, M, a);
        }
    }

    std::cout << "[INFO] 2D Kerr-Schild metrics initialized successfully.\n";

    // 4. Initial Conditions: Magnetized Accretion Ring around Black Hole
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double xc = grid.x[i];
            double yc = grid.y[j];
            double r_cyl = std::sqrt(xc*xc + yc*yc);

            if (r_cyl < 2.0) {
                // Horizon / Central region atmosphere floor
                grid.prim[i][j] = PrimVar(1e-4, 1e-5, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
            } else if (r_cyl >= 3.5 && r_cyl <= 7.5) {
                // Accretion torus ring with orbital velocity
                double v_phi = 0.35;
                double vx = -v_phi * (yc / r_cyl);
                double vy =  v_phi * (xc / r_cyl);
                grid.prim[i][j] = PrimVar(1.0, 0.5, {vx, vy, 0.0}, {0.0, 0.1, 0.0});
            } else {
                // Low density ambient background
                grid.prim[i][j] = PrimVar(0.01, 0.001, {0.0, 0.0, 0.0}, {0.0, 0.01, 0.0});
            }

            grid.cons[i][j] = prim_to_cons(grid.prim[i][j], cell_metrics[i][j]);
        }
    }

    // 5. 2D Time Stepping & TVD Spatial Reconstruction Loop
    double cfl = 0.25;
    double dt = cfl * std::min(grid.dx, grid.dy);
    double t_final = 0.3;
    double t = 0.0;
    int step = 0;

    std::vector<std::vector<PrimVar>> prim_L_x, prim_R_x;
    std::vector<std::vector<PrimVar>> prim_L_y, prim_R_y;

    std::vector<std::vector<ConsVar>> fluxes_x(nx + 1, std::vector<ConsVar>(ny));
    std::vector<std::vector<ConsVar>> fluxes_y(nx, std::vector<ConsVar>(ny + 1));

    std::cout << "[INFO] Running 2D GRMHD time-stepping loop... dt = " << dt << "\n";

    while (t < t_final) {
        if (t + dt > t_final) {
            dt = t_final - t;
        }

        // Step 5a: Outflow Boundary Conditions
        grid.apply_outflow_boundary_conditions();

        // Step 5b: 2D Spatial Reconstruction (MC Limiter)
        reconstruct_2d_tvdslope(grid.prim, prim_L_x, prim_R_x, prim_L_y, prim_R_y, true);

        // Step 5c: Compute 2D X-direction fluxes
        for (int i = 1; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                fluxes_x[i][j] = compute_hll_flux_x(prim_L_x[i][j], prim_R_x[i][j], interface_metrics_x[i][j]);
            }
        }

        // Step 5d: Compute 2D Y-direction fluxes
        for (int i = 0; i < nx; ++i) {
            for (int j = 1; j < ny; ++j) {
                fluxes_y[i][j] = compute_hll_flux_y(prim_L_y[i][j], prim_R_y[i][j], interface_metrics_y[i][j]);
            }
        }

        // Boundary fluxes (outflow)
        for (int j = 0; j < ny; ++j) {
            fluxes_x[0][j] = fluxes_x[1][j];
            fluxes_x[nx][j] = fluxes_x[nx - 1][j];
        }
        for (int i = 0; i < nx; ++i) {
            fluxes_y[i][0] = fluxes_y[i][1];
            fluxes_y[i][ny] = fluxes_y[i][ny - 1];
        }

        // Step 5e: Update Conserved Variables in 2D Grid
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double dtdx = dt / grid.dx;
                double dtdy = dt / grid.dy;

                grid.cons[i][j].D -= dtdx * (fluxes_x[i+1][j].D - fluxes_x[i][j].D)
                                  + dtdy * (fluxes_y[i][j+1].D - fluxes_y[i][j].D);

                for (int k = 0; k < 3; ++k) {
                    grid.cons[i][j].S[k] -= dtdx * (fluxes_x[i+1][j].S[k] - fluxes_x[i][j].S[k])
                                          + dtdy * (fluxes_y[i][j+1].S[k] - fluxes_y[i][j].S[k]);

                    grid.cons[i][j].B[k] -= dtdx * (fluxes_x[i+1][j].B[k] - fluxes_x[i][j].B[k])
                                          + dtdy * (fluxes_y[i][j+1].B[k] - fluxes_y[i][j].B[k]);
                }

                grid.cons[i][j].tau -= dtdx * (fluxes_x[i+1][j].tau - fluxes_x[i][j].tau)
                                     + dtdy * (fluxes_y[i][j+1].tau - fluxes_y[i][j].tau);

                // Update Primitive Variables
                grid.prim[i][j].rho = std::max(1e-4, grid.cons[i][j].D / cell_metrics[i][j].sqrt_minus_g);
                grid.prim[i][j].press = std::max(1e-4, (grid.cons[i][j].tau) * (4.0/3.0 - 1.0));
            }
        }

        t += dt;
        step++;
    }

    std::cout << "[INFO] 2D Simulation completed in " << step << " steps. Final time t = " << t << "\n";

    // 6. Output 2D simulation results to output_2d.csv
    std::ofstream file("output_2d.csv");
    if (file.is_open()) {
        file << "x,y,rho,press,D,tau,alpha\n";
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                file << std::fixed << std::setprecision(6)
                     << grid.x[i] << ","
                     << grid.y[j] << ","
                     << grid.prim[i][j].rho << ","
                     << grid.prim[i][j].press << ","
                     << grid.cons[i][j].D << ","
                     << grid.cons[i][j].tau << ","
                     << cell_metrics[i][j].alpha << "\n";
            }
        }
        file.close();
        std::cout << "[INFO] 2D simulation data successfully exported to output_2d.csv\n";
    } else {
        std::cerr << "[ERROR] Could not open output_2d.csv for writing.\n";
    }

    return 0;
}
