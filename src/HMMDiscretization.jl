module HMMDiscretization

using LinearAlgebra: diagm, tr, norm, cond, diag
using Statistics: mean, var, cov, quantile
using Parameters: @with_kw, @unpack # For keywords in types
using DocStringExtensions: FIELDS, TYPEDSIGNATURES, TYPEDEF # For easier documentation
using Distributions

export 
    HMMContinuousSpaceModel,
    HMMPreallocatedContainers,
    HMMNumericalParameters

include("time_invariant.jl")

end