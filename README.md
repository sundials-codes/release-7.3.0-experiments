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
$
```

## Building

To build the SUNDIALS-based C/C++ experiments:

```bash
git clone https://github.com/sundials-codes/release-7.3.0-experiments.git && cd release-7.3.0-experiments
cmake -S . -B builddir
cmake --build builddir -j8 
```

The Julia experiment is self-contained and will install the Julia packages it requires.

## Running

