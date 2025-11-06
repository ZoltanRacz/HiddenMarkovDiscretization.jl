using Parameters: @with_kw, @unpack
using Distributions
using Random
using LinearAlgebra

@with_kw struct VAR{T<:AbstractFloat} <: HMDContinuousSpaceModel
    "coefficient matrix"
    B::Array{T,2}
    "error terms"
    Σ::Array{T,2} = Matrix(Diagonal(ones(size(B, 1))))
end

struct VARPrealCont{T<:AbstractFloat} <: HMDPreallocatedContainers
    "y_{t-1}"
    y0::Vector{T}
    "without shock"
    ybar::Vector{T}
    "shocks"
    ε::Vector{T}
end

function HiddenMarkovDiscretization.dimnum(model::VAR)
    return size(model.B, 1)
end

function HiddenMarkovDiscretization.HMDPreallocatedContainers(model::VAR)
    return VARPrealCont(zeros(dimnum(model)), zeros(dimnum(model)), zeros(dimnum(model)))
end

function HiddenMarkovDiscretization.simulate_continuous!(sim::AbstractArray, prealcont::VARPrealCont, model::VAR, n::Integer, t::Integer)
    @unpack B, Σ = model
    @unpack y0, ybar, ε = prealcont
    k = dimnum(model)
    errordistr = MvNormal(Σ)

    for l in 1:k
        y0[l] = sim[l, n, t-1]
    end

    mul!(ybar, B, y0)
    
    rand!(errordistr, ε)
    
    for l in 1:k
        sim[l, n, t] = ybar[l] + ε[l]
    end
end

