# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# This example solves the Lotka-Volterra ODE with four parameters,
#     u = [dx/dt] = [ p_0*x - p_1*x*y  ]
#         [dy/dt]   [ -p_2*y + p_3*x*y ].
#
# The initial condition is u(t_0) = 1.0 and we use the parameters
# p  = [1.5, 1.0, 3.0, 1.0].
# ---------------------------------------------------------------

using Pkg

Pkg.add([
    PackageSpec(name="OrdinaryDiffEq", version="6.58.2"),
    PackageSpec(name="SciMLSensitivity", version="7.46.0"),
    PackageSpec(name="ForwardDiff", version="0.10.38"),
    PackageSpec(name="Zygote", version="0.6.77"),
    PackageSpec(name="Plots", version="1.40.13")
])

using LinearAlgebra, OrdinaryDiffEq, SciMLSensitivity, ForwardDiff, Zygote, Plots

# Lotka Volterra with 4 parameters
function f(du, u, p, t)
    du[1] = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = -p[3] * u[2] + p[4] * u[1] * u[2]
end

function g(u, p)
    sum(((1 .- u) .^ 2) ./ 2)
end

function dgdu_continuous!(out, u, p, t)
    # g = sum(((1 .- u) .^ 2) ./ 2)
    out[1] = u[1] - 1.0
    out[2] = u[2] - 1.0
end

function dgdu_discrete!(out, u, p, t, i)
    out[1] = u[1] - 1.0
    out[2] = u[2] - 1.0
end

function dgdp!(out, u, p, t)
    out[1] = 0
    out[2] = 0
    out[3] = 0
    out[4] = 0
end

# --------- Setup Forward Problem

# Problem specification
p = [1.5, 1.0, 3.0, 1.0]
u0 = [1.0; 1.0]
tspan = (0.0, 10.0)
rtol = 1e-14
atol = 1e-14
method = Tsit5()

# Integrate with OrdinaryDiffEq
prob = ODEProblem(f, u0, tspan, p)
sol = solve(prob, method, abstol=atol, reltol=rtol)
println("OrdinaryDiffEq computed ||u(t_f)||: ", norm(sol.u[end]))

# Plot forward solution
plot(sol)
savefig("lotka_volterra_plot.png")

# --------- Setup Adjoint Problem

ts = tspan

# Solve adjoint problem with SciMLSensitivity
res1 = adjoint_sensitivities(sol, method, t=ts, dgdu_discrete=dgdu_discrete!, abstol=atol, reltol=rtol)
println("Discrete SciMLSensitivity computed sensitivities L2 norm: ", norm(res1))

# Solve adjoint problem with forward-mode automatic differentiation
function G(up)
    tmp_prob = remake(prob, u0=up[1:2], p=up[3:end])
    sol = solve(tmp_prob, method, abstol=atol, reltol=rtol, saveat=ts,
        sensealg=SensitivityADPassThrough())
    A = convert(Array, sol)
    g(A, up[3:end])
end
res2 = ForwardDiff.gradient(G, [u0; p])
println("Discrete ForwardDiff computed sensitivities L2 norm: ", norm(res2[1:2]), " ", norm(res2[3:end]))

# Solve adjoint problem with reverse-mode automatic differentiation
function G(u, p)
    tmp_prob = remake(prob, u0=u, p=p)
    sol = solve(tmp_prob, method, abstol=atol, reltol=rtol, saveat=ts)
    A = convert(Array, sol)
    g(A, p)
end
res3 = Zygote.gradient(G, u0, p)
println("Discrete Zygote (reverse mode) computed sensitivities L2 norm: ", norm(res3[1]), " ", norm(res3[2]))

