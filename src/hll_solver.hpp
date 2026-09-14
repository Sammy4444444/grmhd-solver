#ifndef HLL_SOLVER_HPP
#define HLL_SOLVER_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include "grid.hpp"
#include "metric.hpp"

namespace GRMHD {

/**
 * @brief Calculate physical flux vector F^x(P) in X-direction under 3+1 metric (lapse alpha, shift beta^i)
 */
inline ConsVar compute_flux_x(const PrimVar& p, const Metric& m, double gamma_eos = 4.0 / 3.0) {
    ConsVar flux;

    double v2 = p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2];
    double W = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v2));
    
    double B2 = p.B[0]*p.B[0] + p.B[1]*p.B[1] + p.B[2]*p.B[2];
    double v_dot_B = p.v[0]*p.B[0] + p.v[1]*p.B[1] + p.v[2]*p.B[2];
    
    double b0 = W * v_dot_B;
    double b2 = B2 / (W*W) + v_dot_B*v_dot_B;
    double P_tot = p.press + 0.5 * b2;

    // Transport velocity in x-direction: v_trans = alpha * vx - beta_x
    double v_trans_x = m.alpha * p.v[0] - m.beta[0];

    ConsVar u = prim_to_cons(p, m, gamma_eos);

    // Density flux F^x(D) = (alpha * vx - beta_x) * D
    flux.D = v_trans_x * u.D;

    // Momentum fluxes F^x(S_i)
    for (int i = 0; i < 3; ++i) {
        double delta_1i = (i == 0) ? 1.0 : 0.0;
        flux.S[i] = v_trans_x * u.S[i] + (m.alpha * P_tot * delta_1i - b0 * (p.B[i] / W + p.v[0] * b0)) * m.sqrt_minus_g;
    }

    // Energy flux F^x(tau)
    flux.tau = v_trans_x * u.tau + (m.alpha * P_tot * p.v[0] - b0 * p.B[0] / W) * m.sqrt_minus_g;

    // Induction equation flux F^x(B^i)
    flux.B[0] = 0.0;
    flux.B[1] = (v_trans_x * p.B[1] - (m.alpha * p.v[1] - m.beta[1]) * p.B[0]) * m.sqrt_minus_g;
    flux.B[2] = (v_trans_x * p.B[2] - (m.alpha * p.v[2] - m.beta[2]) * p.B[0]) * m.sqrt_minus_g;

    return flux;
}

/**
 * @brief Calculate physical flux vector F^y(P) in Y-direction under 3+1 metric (lapse alpha, shift beta^i)
 */
inline ConsVar compute_flux_y(const PrimVar& p, const Metric& m, double gamma_eos = 4.0 / 3.0) {
    ConsVar flux;

    double v2 = p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2];
    double W = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v2));
    
    double B2 = p.B[0]*p.B[0] + p.B[1]*p.B[1] + p.B[2]*p.B[2];
    double v_dot_B = p.v[0]*p.B[0] + p.v[1]*p.B[1] + p.v[2]*p.B[2];
    
    double b0 = W * v_dot_B;
    double b2 = B2 / (W*W) + v_dot_B*v_dot_B;
    double P_tot = p.press + 0.5 * b2;

    // Transport velocity in y-direction: v_trans = alpha * vy - beta_y
    double v_trans_y = m.alpha * p.v[1] - m.beta[1];

    ConsVar u = prim_to_cons(p, m, gamma_eos);

    // Density flux F^y(D) = (alpha * vy - beta_y) * D
    flux.D = v_trans_y * u.D;

    // Momentum fluxes F^y(S_i)
    for (int i = 0; i < 3; ++i) {
        double delta_2i = (i == 1) ? 1.0 : 0.0;
        flux.S[i] = v_trans_y * u.S[i] + (m.alpha * P_tot * delta_2i - b0 * (p.B[i] / W + p.v[1] * b0)) * m.sqrt_minus_g;
    }

    // Energy flux F^y(tau)
    flux.tau = v_trans_y * u.tau + (m.alpha * P_tot * p.v[1] - b0 * p.B[1] / W) * m.sqrt_minus_g;

    // Induction equation flux F^y(B^i)
    flux.B[0] = (v_trans_y * p.B[0] - (m.alpha * p.v[0] - m.beta[0]) * p.B[1]) * m.sqrt_minus_g;
    flux.B[1] = 0.0;
    flux.B[2] = (v_trans_y * p.B[2] - (m.alpha * p.v[2] - m.beta[2]) * p.B[1]) * m.sqrt_minus_g;

    return flux;
}

/**
 * @brief Estimate local fast magnetosonic wave speeds transformed into 3+1 coordinate frame along X-axis
 */
inline std::pair<double, double> get_characteristic_speeds_x(const PrimVar& p, const Metric& m, double gamma_eos = 4.0 / 3.0) {
    double v2 = p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2];
    v2 = std::min(0.999999, std::max(0.0, v2));
    double W = 1.0 / std::sqrt(1.0 - v2);

    double h = 1.0 + (gamma_eos / (gamma_eos - 1.0)) * (p.press / std::max(1e-10, p.rho));
    double cs2 = (gamma_eos * p.press) / std::max(1e-10, p.rho * h);
    cs2 = std::min(0.999999, std::max(1e-10, cs2));

    double B2 = p.B[0]*p.B[0] + p.B[1]*p.B[1] + p.B[2]*p.B[2];
    double v_dot_B = p.v[0]*p.B[0] + p.v[1]*p.B[1] + p.v[2]*p.B[2];
    double b2 = B2 / (W * W) + v_dot_B * v_dot_B;

    double vA2 = b2 / (p.rho * h + b2);
    double vf2 = cs2 + vA2 * (1.0 - cs2);
    double vf = std::sqrt(std::min(0.999999, std::max(1e-10, vf2)));

    double vx = p.v[0];
    double lambda_plus = (vx + vf) / (1.0 + vx * vf);
    double lambda_minus = (vx - vf) / (1.0 - vx * vf);

    double lambda_plus_coord = m.alpha * lambda_plus - m.beta[0];
    double lambda_minus_coord = m.alpha * lambda_minus - m.beta[0];

    return {lambda_minus_coord, lambda_plus_coord};
}

/**
 * @brief Estimate local fast magnetosonic wave speeds transformed into 3+1 coordinate frame along Y-axis
 */
inline std::pair<double, double> get_characteristic_speeds_y(const PrimVar& p, const Metric& m, double gamma_eos = 4.0 / 3.0) {
    double v2 = p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2];
    v2 = std::min(0.999999, std::max(0.0, v2));
    double W = 1.0 / std::sqrt(1.0 - v2);

    double h = 1.0 + (gamma_eos / (gamma_eos - 1.0)) * (p.press / std::max(1e-10, p.rho));
    double cs2 = (gamma_eos * p.press) / std::max(1e-10, p.rho * h);
    cs2 = std::min(0.999999, std::max(1e-10, cs2));

    double B2 = p.B[0]*p.B[0] + p.B[1]*p.B[1] + p.B[2]*p.B[2];
    double v_dot_B = p.v[0]*p.B[0] + p.v[1]*p.B[1] + p.v[2]*p.B[2];
    double b2 = B2 / (W * W) + v_dot_B * v_dot_B;

    double vA2 = b2 / (p.rho * h + b2);
    double vf2 = cs2 + vA2 * (1.0 - cs2);
    double vf = std::sqrt(std::min(0.999999, std::max(1e-10, vf2)));

    double vy = p.v[1];
    double lambda_plus = (vy + vf) / (1.0 + vy * vf);
    double lambda_minus = (vy - vf) / (1.0 - vy * vf);

    double lambda_plus_coord = m.alpha * lambda_plus - m.beta[1];
    double lambda_minus_coord = m.alpha * lambda_minus - m.beta[1];

    return {lambda_minus_coord, lambda_plus_coord};
}

/**
 * @brief HLL Riemann Solver for flux computation across X-interface
 */
inline ConsVar compute_hll_flux_x(const PrimVar& pL, const PrimVar& pR, const Metric& m, double gamma_eos = 4.0 / 3.0) {
    ConsVar uL = prim_to_cons(pL, m, gamma_eos);
    ConsVar uR = prim_to_cons(pR, m, gamma_eos);

    ConsVar fL = compute_flux_x(pL, m, gamma_eos);
    ConsVar fR = compute_flux_x(pR, m, gamma_eos);

    auto wave_L = get_characteristic_speeds_x(pL, m, gamma_eos);
    auto wave_R = get_characteristic_speeds_x(pR, m, gamma_eos);

    double S_L = std::min(0.0, std::min(wave_L.first, wave_R.first));
    double S_R = std::max(0.0, std::max(wave_L.second, wave_R.second));

    double diff = S_R - S_L;
    if (diff < 1e-8) {
        S_L -= 0.5 * (1e-8 - diff);
        S_R += 0.5 * (1e-8 - diff);
    }

    ConsVar flux_hll;
    double inv_denom = 1.0 / (S_R - S_L);

    flux_hll.D = (S_R * fL.D - S_L * fR.D + S_L * S_R * (uR.D - uL.D)) * inv_denom;

    for (int i = 0; i < 3; ++i) {
        flux_hll.S[i] = (S_R * fL.S[i] - S_L * fR.S[i] + S_L * S_R * (uR.S[i] - uL.S[i])) * inv_denom;
        flux_hll.B[i] = (S_R * fL.B[i] - S_L * fR.B[i] + S_L * S_R * (uR.B[i] - uL.B[i])) * inv_denom;
    }

    flux_hll.tau = (S_R * fL.tau - S_L * fR.tau + S_L * S_R * (uR.tau - uL.tau)) * inv_denom;

    return flux_hll;
}

/**
 * @brief HLL Riemann Solver for flux computation across Y-interface
 */
inline ConsVar compute_hll_flux_y(const PrimVar& pL, const PrimVar& pR, const Metric& m, double gamma_eos = 4.0 / 3.0) {
    ConsVar uL = prim_to_cons(pL, m, gamma_eos);
    ConsVar uR = prim_to_cons(pR, m, gamma_eos);

    ConsVar fL = compute_flux_y(pL, m, gamma_eos);
    ConsVar fR = compute_flux_y(pR, m, gamma_eos);

    auto wave_L = get_characteristic_speeds_y(pL, m, gamma_eos);
    auto wave_R = get_characteristic_speeds_y(pR, m, gamma_eos);

    double S_L = std::min(0.0, std::min(wave_L.first, wave_R.first));
    double S_R = std::max(0.0, std::max(wave_L.second, wave_R.second));

    double diff = S_R - S_L;
    if (diff < 1e-8) {
        S_L -= 0.5 * (1e-8 - diff);
        S_R += 0.5 * (1e-8 - diff);
    }

    ConsVar flux_hll;
    double inv_denom = 1.0 / (S_R - S_L);

    flux_hll.D = (S_R * fL.D - S_L * fR.D + S_L * S_R * (uR.D - uL.D)) * inv_denom;

    for (int i = 0; i < 3; ++i) {
        flux_hll.S[i] = (S_R * fL.S[i] - S_L * fR.S[i] + S_L * S_R * (uR.S[i] - uL.S[i])) * inv_denom;
        flux_hll.B[i] = (S_R * fL.B[i] - S_L * fR.B[i] + S_L * S_R * (uR.B[i] - uL.B[i])) * inv_denom;
    }

    flux_hll.tau = (S_R * fL.tau - S_L * fR.tau + S_L * S_R * (uR.tau - uL.tau)) * inv_denom;

    return flux_hll;
}

} // namespace GRMHD

#endif // HLL_SOLVER_HPP
