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
    "number of individuals to simulate"
    N::S
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
    "Kullback-Leibler divergence"
    KL::Vector{T}
end

KL(dp::HMMDiscretizedParameters) = dp.KL[1]

function simulate_continuous(model::HMMContinuousSpaceModel, numpar::HMMNumericalParameters; prealcont::HMMPreallocatedContainers=HMMPreallocatedContainers(model), y0::AbstractArray=zeros(dimnum(model), numpar.N))
    @unpack T, T0, N = numpar
    k = dimnum(model)
    sim = Array{Float64}(undef, k, N, T + T0)
    for n in 1:N
        for l in 1:k
            sim[l, n, 1] = y0[l, n]
        end
    end
    for t in 2:(T+T0)
        for n in 1:N
            simulate_continuous!(sim, prealcont, model, n, t)
        end
    end

    return sim[:, :, (T0+1):end]
end

function dimnum(model::HMMContinuousSpaceModel)
    throw(error("a separate method has to be written for $(typeof(model))"))
end

function simulate_continuous!(sim::AbstractArray, prealcont::HMMPreallocatedContainers, model::HMMContinuousSpaceModel, n::Integer, t::Integer)
    throw(error("a separate method has to be written for $(typeof(model))"))
end

function default_dp(numpar::HMMNumericalParameters, ys::AbstractArray)
    @unpack m = numpar

    ys_flat = hcat([ys[:, n, :] for n in axes(ys, 2)]...)
    clust = kmeans(ys_flat, m)
    d_state = assignments(clust)
    state_count = counts(clust)

    μ = Matrix(clust.centers')

    Π = zeros(m, m)
    inds = LinearIndices((size(ys, 3), size(ys, 2)))

    for n in axes(ys, 2)
        state_count[d_state[inds[1, n]]] -= 1

        for t in 2:size(ys, 3)
            i = inds[t, n]
            Π[d_state[i-1], d_state[i]] += 1
        end
    end

    for mi in 1:m
        Π[mi, :] = Π[mi, :] ./ state_count[mi]
    end

    σ = vec(sqrt.(mean((ys_flat .- clust.centers[:, d_state]) .^ 2, dims=2)))

    return HMMDiscretizedParameters(Π, μ, σ, [999.0])
end

function discretization(numpar::HMMNumericalParameters, ys::AbstractArray; dp_prev::HMMDiscretizedParameters=default_dp(numpar, ys))
    @unpack maxiter, m, ϵ, N = numpar

    k = size(ys, 1)
    T = size(ys, 3)
    dp_next = deepcopy(dp_prev)
    αs = Array{Float64,3}(undef, m, N, T)
    βs = Array{Float64,3}(undef, m, N, T)
    γs = Array{Float64,3}(undef, m, N, T)
    φs = fill(Normal(), (m, k))
    δs = Matrix{Float64}(undef, (1,m))

    iter = 1
    dif = 100.0

    while iter < maxiter && dif > ϵ
        println("iter is $iter,dif is $dif")
        E_step!(αs, βs, γs, φs, δs, dp_prev, ys)
        M_step!(dp_next, αs, βs, γs, φs, δs, dp_prev, ys)
        dif = abs(KL(dp_next) - KL(dp_prev))
        dp_prev = deepcopy(dp_next)
        println(dp_next)
        iter += 1
    end

    return dp_next
end

function φ(φs, ys, mi, n, t)
    a = 1.0
    for ki in axes(ys, 1)
        a *= pdf(φs[mi, ki], ys[ki, n, t])
    end
    return a
end

function E_step!(αs, βs, γs, φs, δs, dp_prev, ys)
    T = size(ys, 3)
    N = size(ys, 2)
    k = size(ys, 1)
    m = size(αs, 1)
    @unpack μ, Π, σ = dp_prev
    for ki in 1:k
        for mi in 1:m
            φs[mi, ki] = Normal(μ[mi, ki], σ[ki])
        end
    end

    mul!(δs, ones(1, m), inv(UniformScaling(1) - Π + ones(m, m)))
    #println((φs))
    for n in 1:N
        for mi in 1:m
            #    println(φ(φs, ys, mi, 1))
            αs[mi, n, 1] = φ(φs, ys, mi, n, 1) * δs[mi]
            βs[mi, n, end] = 1.0
        end
    end
    #println(αs[:,1])
    for t in 1:(T-1)
        for n in 1:N
            for j in 1:m
                a = 0.0
                for k in 1:m
                    a += αs[k, n, t] * Π[k, j]
                end
                αs[j, n, t+1] = a * φ(φs, ys, j, n, t + 1)
            end
        end
    end
    #println(αs[:,40])
    for t in (T-1):-1:1
        for n in 1:N
            for k in 1:m
                b = 0.0
                for j in 1:m
                    b += βs[j, n, t+1] * Π[k, j] * φ(φs, ys, j, n, t + 1)
                end
                βs[k, n, t] = b
            end
        end
    end
    for t in 1:T
        for n in 1:N
            γsums = 0.0
            for k in 1:m
                γs[k, n, t] = αs[k, n, t] * βs[k, n, t]
                γsums += γs[k, n, t]
            end
            for k in 1:m
                γs[k, n, t] /= γsums
            end
        end
    end
end

function M_step!(dp_next, αs, βs, γs, φs, δs, dp_prev, ys)
    T = size(ys, 3)
    N = size(ys, 2)
    k = size(ys, 1)
    m = size(αs, 1)

    for ki in 1:k
        for mi in 1:m
            μ_num = 0.0
            μ_denom = 0.0
            for t in 1:T
                for n in 1:N
                    μ_num += ys[ki, n, t] * γs[mi, n, t]
                    μ_denom += γs[mi, n, t]
                end
            end
            dp_next.μ[mi, ki] = μ_num / μ_denom
        end
    end

    for ki in 1:k
        a = 0.0
        for mi in 1:m
            for t in 1:T
                for n in 1:N
                    a += (ys[ki, n, t] - dp_next.μ[mi, ki])^2 * γs[mi, n, t]
                end
            end
        end
        dp_next.σ[ki] = sqrt(a / (T * N))
    end

    for k in 1:m
        rowsum = 0.0
        for j in 1:m
            dp_next.Π[k, j] = 0.0
            for t in 1:(T-1)
                for n in 1:N
                    dp_next.Π[k, j] += βs[j, n, t+1] * αs[k, n, t] * dp_prev.Π[k, j] * φ(φs, ys, j, n, t + 1)
                end
            end
            rowsum += dp_next.Π[k, j]
        end
        for j in 1:m
            dp_next.Π[k, j] /= rowsum
        end
    end

    ll = 0.0
    for t in 1:T
        for n in 1:N
            llt = 0.0
            for mi in 1:m
                llt += δs[mi] * φ(φs, ys, mi, n, t)
            end
            ll += log(max(llt, 10^-9))
        end        
    end
    dp_next.KL[1] = -ll / (T * N)
end