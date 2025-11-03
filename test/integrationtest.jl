# useful lines for testing manually, while developing. Install TestEnv in your main environment. When running the first time, activate and instantiate the test environment before restarting Julia and using TestEnv. For more info check: https://github.com/JuliaTesting/TestEnv.jl/blob/main/README.md
#using TestEnv
#TestEnv.activate()
using HMMDiscretization
using Test

include("examples/VAR.jl")

cp = VAR(B = [0.9 0.0; 0.0 0.9], Σ = [1.0 0.0; 0.0 1.0])

#cp = VAR(B = [0.4 0.4; 0.4 0.4], Σ = [1.0 0.0; 0.0 1.0])

np = HMMNumericalParameters(T = 10^3, m = 11, T0 = 100)

sim = simulate_continuous(cp, np)

d = discretization(np,sim)


@test sim isa Array{Float64}

#using Plots
#scatter(d.μ)
#plot(sim')
