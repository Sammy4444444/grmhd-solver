#ifndef METRIC_HPP
#define METRIC_HPP

#include <cmath>
#include <array>
#include <algorithm>

namespace GRMHD {

/**
 * @brief Represents 3+1 ADM Metric & 4D Spacetime parameters.
 * Supports Flat Minkowski spacetime and General Relativistic Kerr-Schild spacetime (spinning BH).
 */
class Metric {
public:
    double alpha;                                   // Lapse function
    std::array<double, 3> beta;                     // Shift vector (beta^i)
    std::array<std::array<double, 3>, 3> gamma;     // Spatial metric (gamma_ij)
    std::array<std::array<double, 3>, 3> gamma_inv; // Inverse spatial metric (gamma^ij)
    double det_gamma;                               // Determinant of spatial metric det(gamma_ij)
    double sqrt_minus_g;                            // 4-metric determinant sqrt(-g) = alpha * sqrt(det_gamma)

    // Full 4-metric tensors
    std::array<std::array<double, 4>, 4> g_4d;       // Covariant 4-metric g_mu_nu
    std::array<std::array<double, 4>, 4> g_4d_inv;   // Contravariant 4-metric g^mu_nu
    
    double M_bh;                                    // Black Hole mass
    double a_spin;                                  // Black Hole dimensionless spin parameter

    Metric() {
        set_minkowski();
    }

    /**
     * @brief Set flat Minkowski spacetime metric in Cartesian coordinates.
     */
    void set_minkowski() {
        alpha = 1.0;
        beta = {0.0, 0.0, 0.0};
        det_gamma = 1.0;
        sqrt_minus_g = 1.0;
        M_bh = 0.0;
        a_spin = 0.0;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                gamma[i][j] = (i == j) ? 1.0 : 0.0;
                gamma_inv[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double eta = (mu == nu) ? ((mu == 0) ? -1.0 : 1.0) : 0.0;
                g_4d[mu][nu] = eta;
                g_4d_inv[mu][nu] = eta;
            }
        }
    }

    /**
     * @brief Set Kerr-Schild metric for a spinning black hole at Cartesian coordinates (x, y, z).
     * @param x Cartesian x-coordinate
     * @param y Cartesian y-coordinate
     * @param z Cartesian z-coordinate
     * @param M Mass of the black hole
     * @param a Spin parameter of the black hole (0 <= |a| < M)
     */
    void set_kerr_schild(double x, double y, double z, double M = 1.0, double a = 0.9) {
        M_bh = M;
        a_spin = a;

        double R2 = x*x + y*y + z*z;
        double a2 = a*a;

        // Solve for Kerr radial coordinate r: r^4 - (R^2 - a^2)*r^2 - a^2*z^2 = 0
        double r2 = 0.5 * ((R2 - a2) + std::sqrt((R2 - a2)*(R2 - a2) + 4.0 * a2 * z*z));
        r2 = std::max(1e-6, r2);
        double r = std::sqrt(r2);

        // Kerr-Schild scalar function f = 2 M r^3 / (r^4 + a^2 z^2)
        double f = (2.0 * M * r * r2) / (r2 * r2 + a2 * z * z);

        // Covariant null vector l_mu = (1, l_x, l_y, l_z)
        double lx = (r * x + a * y) / (r2 + a2);
        double ly = (r * y - a * x) / (r2 + a2);
        double lz = z / r;
        std::array<double, 4> l_down = {1.0, lx, ly, lz};

        // Contravariant null vector l^mu = (-1, l_x, l_y, l_z)
        std::array<double, 4> l_up = {-1.0, lx, ly, lz};

        // Full 4-metric g_mu_nu = eta_mu_nu + f * l_mu * l_nu
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double eta = (mu == nu) ? ((mu == 0) ? -1.0 : 1.0) : 0.0;
                g_4d[mu][nu] = eta + f * l_down[mu] * l_down[nu];
            }
        }

        // Full inverse 4-metric g^mu_nu = eta^mu_nu - f * l^mu * l^nu
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double eta = (mu == nu) ? ((mu == 0) ? -1.0 : 1.0) : 0.0;
                g_4d_inv[mu][nu] = eta - f * l_up[mu] * l_up[nu];
            }
        }

        // ADM 3+1 variables
        alpha = 1.0 / std::sqrt(1.0 + f);
        det_gamma = 1.0 + f;
        sqrt_minus_g = 1.0; // In Kerr-Schild coordinates det(g) = -1, so sqrt(-g) = 1.0

        for (int i = 0; i < 3; ++i) {
            beta[i] = (f / (1.0 + f)) * l_down[i + 1];
            for (int j = 0; j < 3; ++j) {
                double delta = (i == j) ? 1.0 : 0.0;
                gamma[i][j] = delta + f * l_down[i + 1] * l_down[j + 1];
                gamma_inv[i][j] = delta - (f / (1.0 + f)) * l_down[i + 1] * l_down[j + 1];
            }
        }
    }

    /**
     * @brief Lower spatial indices: v_i = gamma_{ij} * v^j
     */
    std::array<double, 3> lower_indices(const std::array<double, 3>& v_up) const {
        std::array<double, 3> v_down = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                v_down[i] += gamma[i][j] * v_up[j];
            }
        }
        return v_down;
    }

    /**
     * @brief Raise spatial indices: v^i = gamma^{ij} * v_j
     */
    std::array<double, 3> raise_indices(const std::array<double, 3>& v_down) const {
        std::array<double, 3> v_up = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                v_up[i] += gamma_inv[i][j] * v_down[j];
            }
        }
        return v_up;
    }
};

} // namespace GRMHD

#endif // METRIC_HPP
