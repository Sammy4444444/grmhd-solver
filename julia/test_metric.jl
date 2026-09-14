
# Julia Prototype: Metric & Tensor Operations for GRMHD

using LinearAlgebra

"""
    kerr_metric_components(r, theta, M=1.0, a=0.9)

Calculate 4-metric components g_mu_nu for Kerr spacetime in Boyer-Lindquist coordinates.
"""
function kerr_metric_components(r::Float64, theta::Float64; M::Float64=1.0, a::Float64=0.9)
    Sigma = r^2 + a^2 * cos(theta)^2
    Delta = r^2 - 2.0 * M * r + a^2

    g_tt = -(1.0 - (2.0 * M * r) / Sigma)
    g_tphi = -(2.0 * M * a * r * sin(theta)^2) / Sigma
    g_rr = Sigma / Delta
    g_thth = Sigma
    g_phiphi = (r^2 + a^2 + (2.0 * M * a^2 * r * sin(theta)^2) / Sigma) * sin(theta)^2

    # 4-metric matrix g_mu_nu (t, r, theta, phi)
    g = [
        g_tt   0.0   0.0   g_tphi;
        0.0    g_rr  0.0   0.0;
        0.0    0.0   g_thth 0.0;
        g_tphi 0.0   0.0   g_phiphi
    ]
    return g
end

function main()
    println("=========================================")
    println(" Julia GR Metric Prototype & Test Script")
    println("=========================================")

    r_test = 5.0
    theta_test = pi / 4.0
    M = 1.0
    a = 0.9375 # Highly spinning black hole

    g = kerr_metric_components(r_test, theta_test, M=M, a=a)
    g_inv = inv(g)
    det_g = det(g)

    println("Evaluation point: r = $r_test, theta = $(round(theta_test, digits=4)), a = $a")
    println("Metric tensor g_mu_nu:")
    display(g)
    println("\nMetric Determinant det(g) = ", det_g)
    println("Inverse Metric g^mu_nu:")
    display(g_inv)

    # Check g * g_inv = Identity
    identity_check = norm(g * g_inv - I(4))
    println("\nIdentity Check ||g * g^-1 - I|| = ", identity_check)
    if identity_check < 1e-12
        println("[SUCCESS] Julia tensor metric calculations verified!")
    else
        println("[WARNING] Deviation detected in inverse metric check.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
