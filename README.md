# SUNDIALS v7.3.0 Release Experiments

This repository contains experiments used in the SUNDIALS v7.3.0 release paper.

## Prerequisites

A C and C++ compiler as well as [CMake](https://cmake.org/) are needed to build the SUNDIALS-based experiments.
For the Julia experiment, the [Julia language](https://julialang.org/) is required.

The specific compilers and library versions used for the release paper are:

GNU gcc/g++:
```shell
$ cc --version
cc (GCC) 8.5.0 20210514 (Red Hat 8.5.0-26)
Copyright (C) 2018 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

Julia:
```shell
$ julia --version
julia version 1.8.3
```

## Building

To build the SUNDIALS-based C/C++ experiments:

```shell
git clone https://github.com/sundials-codes/release-7.3.0-experiments.git && cd release-7.3.0-experiments
cmake -S . -B builddir
cmake --build builddir -j8 
```

The Julia experiment is self-contained and will install the Julia packages it requires.

## Running

### Gray-Scott

To run the Gray-Scott problem as done for the paper:

```shell
$ 
```

### Lotka Volterra

To run the ARKODE version as done for the paper:

```shell
$ cd builddir
$ ./lotka_volterra_arkode
Initial condition:
 1.000000000000000e+00
 1.000000000000000e+00
Forward Solution:
 1.026344767571519e+00
 9.096910781383613e-01
ARKODE Stats for Forward Solution:
Current time                  = 10.0009999999999
Steps                         = 10001
Step attempts                 = 10001
Stability limited steps       = 0
Accuracy limited steps        = 0
Error test fails              = 0
NLS step fails                = 0
Inequality constraint fails   = 0
Initial step size             = 0.001
Last step size                = 0.001
Current step size             = 0.001
Explicit RHS fn evals         = 40005
Implicit RHS fn evals         = 0
NLS iters                     = 0
NLS fails                     = 0
NLS iters per step            = 0
LS setups                     = 0

Adjoint terminal condition:
 2.634476757151916e-02
-9.030892186163875e-02
 0.000000000000000e+00
 0.000000000000000e+00
 0.000000000000000e+00
 0.000000000000000e+00
Adjoint Solution:
 2.990884811622611e-01
-1.393393475855864e-02
 6.209407790826785e-01
 6.811859220563070e-02
 1.690058792635961e-01
 2.721500531537252e-01

SUNAdjointStepper Stats:
Num backwards steps           = 10001
Num recompute passes          = 5000
```

To run the CVODES version as done for the paper:
```shell
$ ./lotka_volterra_cvodes
Forward Solution at t = 10:
 1.026344801589461e+00
 9.096910639240370e-01
CVODES Stats for Forward Solution:
Current time                  = 10.0005022384865
Steps                         = 1542
Error test fails              = 12
NLS step fails                = 0
Initial step size             = 3.74913085977546e-06
Last step size                = 0.0065673662835467
Current step size             = 0.0065673662835467
Last method order             = 5
Current method order          = 5
Stab. lim. order reductions   = 0
RHS fn evals                  = 2640
NLS iters                     = 2637
NLS fails                     = 0
NLS iters per step            = 1.71011673151751
LS setups                     = 0
Jac fn evals                  = 0
LS RHS fn evals               = 1555
Prec setup evals              = 0
Prec solves                   = 0
LS iters                      = 1555
LS fails                      = 0
Jac-times setups              = 0
Jac-times evals               = 1555
LS iters per NLS iter         = 0.58968524838832
Jac evals per NLS iter        = 0
Prec evals per NLS iter       = 0
Root fn evals                 = 0
Adjoint terminal condition:
 2.634480158946095e-02
-9.030893607596302e-02
 0.000000000000000e+00
 0.000000000000000e+00
 0.000000000000000e+00
 0.000000000000000e+00
Adjoint Solution at t = 0:
 2.995996061359468e-01
-1.415181444744058e-02
 6.219617379471348e-01
 6.800142024299936e-02
 1.691765073199833e-01
 2.725607544066027e-01
CVODES Stats for Adjoint Integration:
Current time                  = 0
Steps                         = 2681
Error test fails              = 32
NLS step fails                = 0
Initial step size             = -1.62278558779604e-13
Last step size                = -3.81920641808941e-05
Current step size             = -3.81920641808941e-05
Last method order             = 4
Current method order          = 5
Stab. lim. order reductions   = 0
RHS fn evals                  = 3494
NLS iters                     = 3493
NLS fails                     = 0
NLS iters per step            = 1.30287206266319
LS setups                     = 0
Jac fn evals                  = 0
LS RHS fn evals               = 2294
Prec setup evals              = 0
Prec solves                   = 0
LS iters                      = 2294
LS fails                      = 0
Jac-times setups              = 0
Jac-times evals               = 2294
LS iters per NLS iter         = 0.656742055539651
Jac evals per NLS iter        = 0
Prec evals per NLS iter       = 0
Root fn evals                 = 0
Quad fn evals                 = 2712
Quad error test fails         = 30
```

To run the Julia version as done for the paper:

```shell
$ cd Lotka-Volterra
$ julia lotka_volterra_zygote.jl
OrdinaryDiffEq computed ||u(t_f)||: 1.3714668933550804
Discrete SciMLSensitivity computed sensitivities L2 norm: 0.7644151699219316
Discrete ForwardDiff computed sensitivities L2 norm: 0.2999335581117122 0.7031147934185299
Discrete Zygote (reverse mode) computed sensitivities L2 norm: 0.2999335581120665 0.703114793419319
```
