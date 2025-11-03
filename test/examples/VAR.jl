using Parameters: @with_kw, @unpack
using Distributions
using Random
using LinearAlgebra

@with_kw struct VAR{T<:AbstractFloat} <: HMMContinuousSpaceModel
    "coefficient matrix"
    B::Array{T,2}
    "error terms"
    Σ::Array{T,2} = Matrix(Diagonal(ones(size(B, 1))))
end

struct VARPrealCont{T<:AbstractFloat} <: HMMPreallocatedContainers
    "y_{t-1}"
    y0::Vector{T}
    "without shock"
    ybar::Vector{T}
    "shocks"
    ε::Vector{T}
end

function HMMDiscretization.dimnum(model::VAR)
    return size(model.B, 1)
end

function HMMDiscretization.HMMPreallocatedContainers(model::VAR)
    return VARPrealCont(zeros(dimnum(model)), zeros(dimnum(model)), zeros(dimnum(model)))
end

function HMMDiscretization.simulate_continuous!(sim::AbstractArray, prealcont::VARPrealCont, model::VAR, t::Integer)
    @unpack B, Σ = model
    @unpack y0, ybar, ε = prealcont
    k = dimnum(model)
    errordistr = MvNormal(Σ)

    for l in 1:k
        y0[l] = sim[l, t-1]
    end

    mul!(ybar, B, y0)
    
    rand!(errordistr, ε)
    
    for l in 1:k
        sim[l, t] = ybar[l] + ε[l]
    end
end

