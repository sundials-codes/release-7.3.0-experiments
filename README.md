# SUNDIALS v7.4.0 Release Experiments

This repository contains experiments used in the SUNDIALS v7.4.0 release paper.

## Prerequisites

A C and C++ compiler, [CMake](https://cmake.org/), and [OpenMP](https://www.openmp.org/) are needed to build the SUNDIALS-based experiments.
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
git clone https://github.com/sundials-codes/release-7.4.0-experiments.git && cd release-7.4.0-experiments
cmake -S . -B builddir
cmake --build builddir -j8 
```

The Julia experiment is self-contained and will install the Julia packages it
requires.

## Running

### Gray–Scott

Slurm batch scripts are provided to run the Gray–Scott experiments. To generate
a reference solution and the operator splitting results, set the CMake option
`SUNDIALS_PRECISION` to `EXTENDED` and execute the following commands:

```shell
cd Gray-Scott
mkdir data
sbatch ref.sh
sbatch splitting.sh
```

Finally, return the `SUNDIALS_PRECISION` to `DOUBLE` for the LSRK experiments
and run

```shell
sbatch lsrk.sh
```

The solutions and runtimes are writtin into the `data` directory. A Python
script is included to convert these into CSV tables:

```shell
python3 postprocess.py --method Splitting --grid_pts_1d 128 --plot
python3 postprocess.py --method LSRK --grid_pts_1d 1024 --plot
```

### Lotka Volterra

To run the ARKODE version as done for the paper:

```shell
$ cd builddir
$ ./run_experiment.sh
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
