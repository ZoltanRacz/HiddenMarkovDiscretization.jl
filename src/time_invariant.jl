"""
$(TYPEDEF)
# Description
Supplies parameters to the continuous space model to be discretized."""
abstract type HMMContinuousSpaceModel end

"""
$(TYPEDEF)
# Description
Supplies containers to simulate continuous space model."""
abstract type HMMPreallocatedContainers end

"""
$(TYPEDEF)
# Description
Stores numerical parameters for the discretization procedure.
# Fields
$(FIELDS)
"""
@with_kw struct HMMNumericalParameters{S<:Integer,U<:AbstractFloat}
    "number of time periods to simulate, excluding burn-in"
    T::S
    "number of time periods for burn-in"
    T0::S
    "number of discrete states"
    m::S
    "tolerance"
    ϵ::U = 10^-7
    "maximal number of iterations"
    maxiter::S = 10^5
end

"""
$(TYPEDEF)
# Description
Stores parameters of a discretized Hidden Markov Model.

 * `k` is the number of dimensions of the stochastic process
 * `m` is the number of desired discrete states
# Fields
$(FIELDS)
"""
@with_kw struct HMMDiscretizedParameters{T<:AbstractFloat}
    "m*m transition matrix between discrete states. Π[k,j] denotes the P(x_{t+1} = j | x_{t} = k)"
    Π::Array{T,2}
    "m*k matrix containing the k-long grid for each discrete states"
    μ::Array{T,2}
    "k-long vector containing the standard errors around the hidden state, for each dimension"
    σ::Vector{T}
end

function simulate_continuous(model::HMMContinuousSpaceModel, numpar::HMMNumericalParameters; prealcont::HMMPreallocatedContainers=HMMPreallocatedContainers(model), y0::AbstractVector=zeros(dimnum(model)))
    @unpack T, T0 = numpar
    k = dimnum(model)
    sim = Array{Float64}(undef,k, T + T0)
    for l in 1:k
        sim[l, 1] = y0[l]
    end
    for t in 2:(T+T0)
        simulate_continuous!(sim, prealcont, model, t)
    end

    return sim[:, (T0+1):end]
end

function dimnum(model::HMMContinuousSpaceModel)
    throw(error("a separate method has to be written for $(typeof(model))"))
end

function simulate_continuous!(sim::AbstractArray, prealcont::HMMPreallocatedContainers, model::HMMContinuousSpaceModel, t::Integer)
    throw(error("a separate method has to be written for $(typeof(model))"))
end

function default_dp(numpar::HMMNumericalParameters, ys::AbstractArray)
    @unpack m = numpar
    clust = kmeans(ys,m)
    d_state = assignments(clust)
    state_count = counts(clust)

    μ = Matrix(clust.centers')

    Π = zeros(m,m)
    for t in 2:size(ys,2)
        Π[d_state[t-1],d_state[t]] += 1
    end
    state_count[d_state[end]] -= 1
    for mi in 1:m
        Π[mi,:] = Π[mi,:]./ state_count[mi]
    end    

    σ = vec(sqrt.(mean((ys .- clust.centers[:,d_state]).^2, dims = 2)))

    return HMMDiscretizedParameters(Π, μ, σ)
end

function discretization(numpar::HMMNumericalParameters, ys::AbstractArray; dp_prev::HMMDiscretizedParameters = default_dp(numpar,ys))
    return dp_prev
    k = size(ys,1)
    xs = similar(ys)

    d
end
