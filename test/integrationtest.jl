using HiddenMarkovDiscretization
using Test

include("examples/VAR.jl")

#cp = VAR(B = ones(1,1)*0.9)

cp = VAR(B = [0.7 0.2; 0.2 0.7], Σ = [1.0 0.0; 0.0 1.0])

np = HMDNumericalParameters(T = 4, N = 10000, m = 5, T0 = 100)

sim = simulate_continuous(cp, np)
@test sim isa Array{Float64, 3}

d0 = HiddenMarkovDiscretization.default_dp(np,sim)

@test d0 isa HiddenMarkovDiscretization.HMDDiscretizedParameters{Float64}
@test maximum(abs.(sum(d0.Π, dims = 2) .- ones(np.m,1)))<10^-5

d = discretization(np,sim)
@test d isa HiddenMarkovDiscretization.HMDDiscretizedParameters{Float64}

simd = HiddenMarkovDiscretization.simulate_discrete(d, np)

mc = HiddenMarkovDiscretization.moment_comparison(sim,simd)
@test abs(mc[1][2,2] - mc[1][2,3]) < 0.05

#using Plots
#scatter(d0.μ[:,1], d0.μ[:,2])
#scatter!(d.μ[:,1], d.μ[:,2])
#scatter(d.μ)
#plot(sim')
#scatter(d0.μ[:, 1])
#scatter!(d.μ[:,1])