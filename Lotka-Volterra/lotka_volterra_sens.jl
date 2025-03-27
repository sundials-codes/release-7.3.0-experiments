# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------

using OrdinaryDiffEq, SciMLSensitivity, LinearAlgebra
using QuadGK, ForwardDiff, Zygote
using Plots

# Lotka Volterra with 4 parameters
function f(du, u, p, t)
    du[1] = dx = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = dy = -p[3] * u[2] + p[4] * u[1] * u[2]
end

function g(u, p)
    # tmp = (sum(u) .^ 2) ./ 2
    # (u[1].^2 + u[2].^2) ./ 2
    sum(((1 .- u) .^ 2) ./ 2)
end

function dgdu!(out, u, p, t)
    # # g = (u[1].^2 + u[2].^2) ./ 2
    # out[1] = u[1] + u[2]
    # out[2] = u[1] + u[2]

    # g = sum(((1 .- u) .^ 2) ./ 2)
    out[1] = u[1] - 1.0
    out[2] = u[2] - 1.0
end

function dgdp!(out, u, p, t)
    out[1] = 0
    out[2] = 0
    out[3] = 0
    out[4] = 0
end

function dgdu_discrete!(out, u, p, t, i)
    dgdu!(out, u, p, t)
end

# --------- Setup Forward Problem

# Problem specification
p = [1.5, 1.0, 3.0, 1.0]
u0 = [1.0; 1.0]
tspan = (0.0, 1.0)
rtol = 1e-14
atol = 1e-14
method = Vern9()

# Solve forward problem
prob = ODEProblem(f, u0, tspan, p)
sol = solve(prob, method, abstol=atol, reltol=rtol)
println("OrdinaryDiffEq computed u(t_f): ", sol.u[end])

plot(sol)
savefig("lotka_volterra_plot.png")

# --------- Setup Adjoint Problem

# # Solve adjoint problem with continuous cost functional
# res1 = adjoint_sensitivities(sol, method, dgdu_continuous=dgdu!, dgdp_continuous=dgdp!, g=g, abstol=atol, reltol=rtol)
# println("Continuous SciMLSensitivity computed sensitivities: ", res1)

# Solve adjoint problem with discrete cost functional
ts = tspan
res1 = adjoint_sensitivities(sol, method, t=ts, dgdu_discrete=dgdu_discrete!, abstol=atol, reltol=rtol)
println("Discrete SciMLSensitivity computed sensitivities: ", res1)

# --------- Validate against ForwardDiff and Zygote

# Discrete functional
function G(up)
    tmp_prob = remake(prob, u0=up[1:2], p=up[3:end])
    sol = solve(tmp_prob, method, abstol=atol, reltol=rtol, saveat=ts,
                sensealg=SensitivityADPassThrough())
    A = convert(Array, sol)
    g(A, up[3:end])
end
res2 = ForwardDiff.gradient(G, [u0; p])
println("Discrete ForwardDiff computed sensitivities: ", res2)

function G(u, p)
    tmp_prob = remake(prob, u0=u, p=p)
    sol = solve(tmp_prob, method, abstol=atol, reltol=rtol, saveat=ts)
    A = convert(Array, sol)
    g(A, p)
end
res3 = Zygote.gradient(G, u0, p)
println("Discrete Zygote (reverse mode) computed sensitivities: ", res3)
