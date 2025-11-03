module HMMDiscretization

using LinearAlgebra
using Statistics: mean, var, cov, quantile
using Parameters: @with_kw, @unpack # For keywords in types
using DocStringExtensions: FIELDS, TYPEDSIGNATURES, TYPEDEF # For easier documentation
using Distributions
using Clustering # for starting grid and transition matrix


export 
    HMMContinuousSpaceModel,
    HMMPreallocatedContainers,
    HMMNumericalParameters,

    dimnum,
    simulate_continuous,
    simulate_continuous!,
    discretization

include("time_invariant.jl")

end