using Parameters: @with_kw, @unpack

@with_kw struct VAR{T<:AbstractFloat} <: HMMContinuousSpaceModel
    "coefficient matrix"
    B::Array{T,2}
    "error terms"
    Σ::Array{T,2} = Identity(size(B, 1))
end

struct VARPrealCont{T<:AbstractFloat} <: HMMPreallocatedContainers
    "y_{t-1}"
    y0::Vector{T}
    "y_{t}"
    y1::Vector{T}
end

function dimnum(model::VAR)
    return size(model.B, 1)
end

function HMMDiscretization.HMMPreallocatedContainers(model::VAR)
    return VARPrealCont(y0=zeros(dimnum(model)), y1=zeros(dimnum(model)))
end

function simulate_continuous!(sim::AbstractArray, prealcont::VARPrealCont, model::VAR, t::Integer)
    @unpack y0, y1 = prealcont
    k = dimnum(model)

    for l in 1:k
        y0[l] = sim[l, t-1]
    end



    for l in 1:k
        sim[l, t] = y1[l] 
    end
end