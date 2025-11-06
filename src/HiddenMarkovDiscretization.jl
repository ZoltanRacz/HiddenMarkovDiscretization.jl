module HiddenMarkovDiscretization

using LinearAlgebra
using Statistics: mean, var, cov, quantile
using Parameters: @with_kw, @unpack # For keywords in types
using DocStringExtensions: FIELDS, TYPEDSIGNATURES, TYPEDEF # For easier documentation
using Distributions
using Clustering # for starting grid and transition matrix
using DataFrames # output format in diagnostics


export 
    HMDContinuousSpaceModel,
    HMDPreallocatedContainers,
    HMDNumericalParameters,

    dimnum,
    simulate_continuous,
    simulate_continuous!,
    discretization

include("time_invariant.jl")
include("diagnostics.jl")

end