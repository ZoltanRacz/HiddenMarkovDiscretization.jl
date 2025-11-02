using HMMDiscretization
using Test

include("examples/VAR.jl")

@testset "Unit test of global estimation" begin
    include("integrationtest.jl")
end