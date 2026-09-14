#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include "grid.hpp"

namespace GRMHD {

/**
 * @brief Minmod slope limiter function for 2 arguments
 */
inline double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

/**
 * @brief Minmod slope limiter function for 3 arguments
 */
inline double minmod3(double a, double b, double c) {
    if (a > 0.0 && b > 0.0 && c > 0.0) {
        return std::min({a, b, c});
    } else if (a < 0.0 && b < 0.0 && c < 0.0) {
        return std::max({a, b, c});
    }
    return 0.0;
}

/**
 * @brief Monotonized Central (MC) slope limiter
 */
inline double mc_limiter(double delta_left, double delta_right) {
    double delta_central = 0.5 * (delta_left + delta_right);
    return minmod3(2.0 * delta_left, delta_central, 2.0 * delta_right);
}

/**
 * @brief Perform TVD 2nd Order Linear Reconstruction on 1D slice
 */
inline void reconstruct_1d_tvdslope(
    const std::vector<PrimVar>& prim,
    std::vector<PrimVar>& prim_L,
    std::vector<PrimVar>& prim_R,
    bool use_mc = true) 
{
    int n = static_cast<int>(prim.size());
    if (static_cast<int>(prim_L.size()) != n + 1) prim_L.resize(n + 1);
    if (static_cast<int>(prim_R.size()) != n + 1) prim_R.resize(n + 1);

    std::vector<PrimVar> slopes(n);

    for (int i = 0; i < n; ++i) {
        int i_prev = std::max(0, i - 1);
        int i_next = std::min(n - 1, i + 1);

        auto calc_slope = [use_mc](double val_left, double val_center, double val_right) {
            double dL = val_center - val_left;
            double dR = val_right - val_center;
            return use_mc ? mc_limiter(dL, dR) : minmod(dL, dR);
        };

        slopes[i].rho = calc_slope(prim[i_prev].rho, prim[i].rho, prim[i_next].rho);
        slopes[i].press = calc_slope(prim[i_prev].press, prim[i].press, prim[i_next].press);

        for (int k = 0; k < 3; ++k) {
            slopes[i].v[k] = calc_slope(prim[i_prev].v[k], prim[i].v[k], prim[i_next].v[k]);
            slopes[i].B[k] = calc_slope(prim[i_prev].B[k], prim[i].B[k], prim[i_next].B[k]);
        }
    }

    for (int i = 0; i < n - 1; ++i) {
        prim_L[i + 1].rho = std::max(1e-6, prim[i].rho + 0.5 * slopes[i].rho);
        prim_L[i + 1].press = std::max(1e-6, prim[i].press + 0.5 * slopes[i].press);

        prim_R[i + 1].rho = std::max(1e-6, prim[i + 1].rho - 0.5 * slopes[i + 1].rho);
        prim_R[i + 1].press = std::max(1e-6, prim[i + 1].press - 0.5 * slopes[i + 1].press);

        for (int k = 0; k < 3; ++k) {
            prim_L[i + 1].v[k] = prim[i].v[k] + 0.5 * slopes[i].v[k];
            prim_L[i + 1].B[k] = prim[i].B[k] + 0.5 * slopes[i].B[k];

            prim_R[i + 1].v[k] = prim[i + 1].v[k] - 0.5 * slopes[i + 1].v[k];
            prim_R[i + 1].B[k] = prim[i + 1].B[k] - 0.5 * slopes[i + 1].B[k];
        }

        auto clamp_v = [](std::array<double, 3>& v) {
            double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
            if (v2 >= 0.99999) {
                double factor = std::sqrt(0.9999 / v2);
                v[0] *= factor;
                v[1] *= factor;
                v[2] *= factor;
            }
        };
        clamp_v(prim_L[i + 1].v);
        clamp_v(prim_R[i + 1].v);
    }

    prim_L[0] = prim[0];
    prim_R[0] = prim[0];
    prim_L[n] = prim[n - 1];
    prim_R[n] = prim[n - 1];
}

/**
 * @brief Perform 2D Spatial Reconstruction for X-interfaces and Y-interfaces
 */
inline void reconstruct_2d_tvdslope(
    const std::vector<std::vector<PrimVar>>& prim,
    std::vector<std::vector<PrimVar>>& prim_L_x,
    std::vector<std::vector<PrimVar>>& prim_R_x,
    std::vector<std::vector<PrimVar>>& prim_L_y,
    std::vector<std::vector<PrimVar>>& prim_R_y,
    bool use_mc = true)
{
    int nx = static_cast<int>(prim.size());
    int ny = static_cast<int>(prim[0].size());

    // Allocate 2D interface state arrays
    prim_L_x.resize(nx + 1, std::vector<PrimVar>(ny));
    prim_R_x.resize(nx + 1, std::vector<PrimVar>(ny));

    prim_L_y.resize(nx, std::vector<PrimVar>(ny + 1));
    prim_R_y.resize(nx, std::vector<PrimVar>(ny + 1));

    // X-direction reconstruction for each row j
    for (int j = 0; j < ny; ++j) {
        std::vector<PrimVar> row(nx);
        for (int i = 0; i < nx; ++i) row[i] = prim[i][j];

        std::vector<PrimVar> row_L(nx + 1);
        std::vector<PrimVar> row_R(nx + 1);
        reconstruct_1d_tvdslope(row, row_L, row_R, use_mc);

        for (int i = 0; i <= nx; ++i) {
            prim_L_x[i][j] = row_L[i];
            prim_R_x[i][j] = row_R[i];
        }
    }

    // Y-direction reconstruction for each column i
    for (int i = 0; i < nx; ++i) {
        std::vector<PrimVar> col(ny);
        for (int j = 0; j < ny; ++j) col[j] = prim[i][j];

        std::vector<PrimVar> col_L(ny + 1);
        std::vector<PrimVar> col_R(ny + 1);
        reconstruct_1d_tvdslope(col, col_L, col_R, use_mc);

        for (int j = 0; j <= ny; ++j) {
            prim_L_y[i][j] = col_L[j];
            prim_R_y[i][j] = col_R[j];
        }
    }
}

} // namespace GRMHD

#endif // RECONSTRUCTION_HPP
