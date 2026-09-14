#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include "metric.hpp"

namespace GRMHD {

/**
 * @brief Primitive Variables structure
 */
struct PrimVar {
    double rho;                     // Rest mass density
    double press;                   // Pressure
    std::array<double, 3> v;        // Spatial velocity (vx, vy, vz)
    std::array<double, 3> B;        // Magnetic field (Bx, By, Bz)

    PrimVar() : rho(1.0), press(1.0), v{0.0, 0.0, 0.0}, B{0.0, 0.0, 0.0} {}
    PrimVar(double r, double p, std::array<double, 3> v_in, std::array<double, 3> b_in)
        : rho(r), press(p), v(v_in), B(b_in) {}
};

/**
 * @brief Conserved Variables structure
 */
struct ConsVar {
    double D;                       // Conserved density: D = gamma_lorentz * rho * sqrt(-g)
    std::array<double, 3> S;        // Conserved momentum vector S_i
    double tau;                     // Conserved energy scalar tau
    std::array<double, 3> B;        // Conserved magnetic field B^i

    ConsVar() : D(0.0), S{0.0, 0.0, 0.0}, tau(0.0), B{0.0, 0.0, 0.0} {}
};

/**
 * @brief Convert Primitive Variables to Conserved Variables under metric m & ideal EOS
 */
inline ConsVar prim_to_cons(const PrimVar& p, const Metric& m, double gamma_eos = 4.0 / 3.0) {
    ConsVar c;
    
    double v2 = p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2];
    v2 = std::min(0.999999, std::max(0.0, v2));
    double W = 1.0 / std::sqrt(1.0 - v2); // Lorentz factor

    double h = 1.0 + gamma_eos / (gamma_eos - 1.0) * (p.press / std::max(1e-10, p.rho)); // Specific enthalpy
    double B2 = p.B[0]*p.B[0] + p.B[1]*p.B[1] + p.B[2]*p.B[2];
    double v_dot_B = p.v[0]*p.B[0] + p.v[1]*p.B[1] + p.v[2]*p.B[2];

    double b0 = W * v_dot_B;
    double b2 = B2 / (W*W) + v_dot_B*v_dot_B;
    double h_tot = h + b2 / std::max(1e-10, p.rho);

    c.D = W * p.rho * m.sqrt_minus_g;

    for (int i = 0; i < 3; ++i) {
        c.S[i] = (p.rho * h_tot * W * W * p.v[i] - b0 * p.B[i]) * m.sqrt_minus_g;
        c.B[i] = p.B[i] * m.sqrt_minus_g;
    }

    double P_tot = p.press + 0.5 * b2;
    c.tau = (p.rho * h_tot * W * W - P_tot - b0 * b0 - c.D) * m.sqrt_minus_g;

    return c;
}

/**
 * @brief Grid class representing 1D / 2D mesh layout
 */
class Grid {
public:
    int nx, ny;
    double x_min, x_max;
    double y_min, y_max;
    double dx, dy;

    std::vector<double> x;
    std::vector<double> y;

    std::vector<std::vector<PrimVar>> prim;
    std::vector<std::vector<ConsVar>> cons;

    Grid(int num_x = 100, int num_y = 100, double xmin = -10.0, double xmax = 10.0, double ymin = -10.0, double ymax = 10.0)
        : nx(num_x), ny(num_y), x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax) {
        
        dx = (x_max - x_min) / nx;
        dy = (ny > 1) ? (y_max - y_min) / ny : 1.0;

        x.resize(nx);
        for (int i = 0; i < nx; ++i) {
            x[i] = x_min + (i + 0.5) * dx;
        }

        y.resize(ny);
        for (int j = 0; j < ny; ++j) {
            y[j] = y_min + (j + 0.5) * dy;
        }

        prim.resize(nx, std::vector<PrimVar>(ny));
        cons.resize(nx, std::vector<ConsVar>(ny));
    }

    /**
     * @brief Apply 2D Outflow (transmissive) boundary conditions along grid edges
     */
    void apply_outflow_boundary_conditions() {
        if (ny <= 1) return;

        for (int j = 0; j < ny; ++j) {
            // X-min and X-max boundaries
            prim[0][j] = prim[1][j];
            prim[nx - 1][j] = prim[nx - 2][j];
        }

        for (int i = 0; i < nx; ++i) {
            // Y-min and Y-max boundaries
            prim[i][0] = prim[i][1];
            prim[i][ny - 1] = prim[i][ny - 2];
        }
    }
};

} // namespace GRMHD

#endif // GRID_HPP
