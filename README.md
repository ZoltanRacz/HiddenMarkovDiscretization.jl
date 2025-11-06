# HiddenMarkovDiscretization.jl

[![Build Status](https://github.com/ZoltanRacz/HiddenMarkovDiscretization.jl/actions/workflows/CI.yml/badge.svg?event=push)](https://github.com/ZoltanRacz/HiddenMarkovDiscretization.jl/actions/workflows/CI.yml?event=push)
[![codecov](https://codecov.io/gh/ZoltanRacz/HiddenMarkovDiscretization.jl/graph/badge.svg?token=SV1YNHSXY1)](https://codecov.io/gh/ZoltanRacz/HiddenMarkovDiscretization.jl)

Efficient finite-state Hidden Markov Model approximation for Markov processes with continuous support. This repository implements the method developed in [Janssens and McCrary (2024)](https://github.com/SeanMcCrary/HMM_Discretization/blob/main/JanssensMcCrary_FiniteState.pdf) into Julia. For a corresponding package written in MATLAB, see the code by the authors at [https://github.com/SeanMcCrary/HMM_Discretization](https://github.com/SeanMcCrary/HMM_Discretization).

## Installation

Install `HiddenMarkovDiscretization.jl` using the Julia package manager:
```julia
import Pkg
Pkg.add("HiddenMarkovDiscretization")
```

## Example