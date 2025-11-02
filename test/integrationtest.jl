using .HMMDiscretization
using Test

include("examples/VAR.jl")

cp = VAR(B = [1.0 0.0; 0.0 1.0], Σ = [1.0 0.0; 0.0 1.0])

np = HMMNumericalParameters(T = 10^7, m = 11, T0 = 100)

