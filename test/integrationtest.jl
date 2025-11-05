# useful lines for testing manually, while developing. Install TestEnv in your main environment. When running the first time, activate and instantiate the test environment before restarting Julia and using TestEnv. For more info check: https://github.com/JuliaTesting/TestEnv.jl/blob/main/README.md
#using TestEnv
#TestEnv.activate()
using HMMDiscretization
using Test

include("examples/VAR.jl")

#cp = VAR(B = ones(1,1)*0.9)

cp = VAR(B = [0.7 0.2; 0.2 0.7], Σ = [1.0 0.0; 0.0 1.0])

np = HMMNumericalParameters(T = 4, N = 10000, m = 15, T0 = 100)

sim = simulate_continuous(cp, np)

d0 = HMMDiscretization.default_dp(np,sim)

@test d0 isa HMMDiscretization.HMMDiscretizedParameters{Float64}
@test maximum(abs.(sum(d0.Π, dims = 2) .- ones(np.m,1)))<10^-5

d = discretization(np,sim)

simd = HMMDiscretization.simulate_discrete(d, np)

simd0 = HMMDiscretization.simulate_discrete(d0, np)

mc = HMMDiscretization.moment_comparison(sim,simd)

mc0 = HMMDiscretization.moment_comparison(sim,simd0)


@test sim isa Array{Float64, 3}
@test d isa HMMDiscretization.HMMDiscretizedParameters{Float64}

#using Plots
#scatter(d0.μ[:,1], d0.μ[:,2])
#scatter!(d.μ[:,1], d.μ[:,2])
#scatter(d.μ)
#plot(sim')
#scatter(d0.μ[:, 1])
#scatter!(d.μ[:,1])