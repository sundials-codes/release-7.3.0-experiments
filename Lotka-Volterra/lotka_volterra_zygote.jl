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
    PackageSpec(name="Printf"),
    PackageSpec(name="ArgParse"),
    PackageSpec(name="OrdinaryDiffEq"),
    PackageSpec(name="SciMLSensitivity"),
    PackageSpec(name="ForwardDiff"),
    PackageSpec(name="Zygote"),
    PackageSpec(name="Plots")
])

using ArgParse, Printf, LinearAlgebra, OrdinaryDiffEq, SciMLSensitivity, ForwardDiff, Zygote, Plots

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

# --------- Command line argument handling

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--order"
        help = "Order of the ODE solver method (e.g., 3 for BS3, 4 for RK4, 5 for Tsit5)"
        arg_type = Int
        default = 3
        "--dt"
        help = "Step size (default=0.0 - which means adaptive)"
        arg_type = Float64
        default = 0.0
    end
    return parse_args(s)
end

args = parse_commandline()
method_order = args["order"]
dt = args["dt"]
adaptive = false

if method_order == 3
    method = BS3()
elseif method_order == 4
    method = RK4()
elseif method_order == 5
    method = Tsit5()
else
    error("Unsupported method order: $method_order. Supported orders are 3 (BS3), 4 (RK4), 5 (Tsit5).")
end

if dt == 0.0
    adaptive = true
else
    adaptive = false
end

# --------- Setup Forward Problem

# Problem specification
p = [1.5, 1.0, 3.0, 1.0]
u0 = [1.0; 1.0]
tspan = (0.0, 10.0)
rtol = 1e-14
atol = 1e-14

# Integrate with OrdinaryDiffEq
prob = ODEProblem(f, u0, tspan, p)
if adaptive
    sol = solve(prob, method, adaptive=adaptive, abstol=atol, reltol=rtol)
else
    sol = solve(prob, method, dt=dt, adaptive=adaptive)
end
println("OrdinaryDiffEq computed solution: ", sol.u[end])
sol_nrm = norm(sol.u[end])

# Plot forward solution
plot(sol)
savefig("lotka_volterra_plot.png")

# --------- Setup Adjoint Problem

ts = tspan

# Solve adjoint problem with reverse-mode automatic differentiation
function G(u, p)
    tmp_prob = remake(prob, u0=u, p=p)
    if adaptive
        sol = solve(tmp_prob, method, adaptive=adaptive, abstol=atol, reltol=rtol, saveat=ts)
    else
        sol = solve(tmp_prob, method, dt=dt, adaptive=adaptive, saveat=ts)
    end
    A = convert(Array, sol)
    g(A, p)
end
res3 = Zygote.gradient(G, u0, p)
@printf("L2 Norm of u, dg/du0, dg/dp: %.16e, %.16e, %.16e\n", sol_nrm, norm(res3[1]), norm(res3[2]))

