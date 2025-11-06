function moment_comparison(sim_c::AbstractArray, sim_d::AbstractArray; unimoments::Bool=true, covariances::Bool=true, max_lag::Integer=1, varnames::Vector{String}=["y" * string(i) for i in axes(sim_c, 1)])
    ds = DataFrame[]

    if unimoments
        sim_c_flat = [vec(sim_c[ki, :, :]) for ki in axes(sim_c, 1)]
        sim_d_flat = [vec(sim_d[ki, :, :]) for ki in axes(sim_d, 1)]
        momentnames = ["Mean", "Standard Deviation", "Skewness", "Kurtosis"]
        funs = [mean, std, skewness, kurtosis]
        rownames = String[]
        momvals_c = Float64[]
        momvals_d = Float64[]
        for ki in eachindex(sim_c_flat)
            for fi in eachindex(funs)
                push!(rownames, "$(momentnames[fi]) of " * varnames[ki])
                push!(momvals_c, funs[fi](sim_c_flat[ki]))
                push!(momvals_d, funs[fi](sim_d_flat[ki]))
            end
        end
        push!(ds, DataFrame(:Moments => rownames, Symbol("Continuous Simulation") => round.(momvals_c, digits=3), Symbol("Discretized Simulation") => round.(momvals_d, digits=3)))
    end

    if covariances
        rownames = String[]
        momvals_c = Any[]
        momvals_d = Any[]
        for lag in 0:max_lag
            sim_c_flat_0 = hcat([sim_c[:, n, (1+lag):end] for n in axes(sim_c, 2)]...)
            sim_c_flat_lag = hcat([sim_c[:, n, 1:(end-lag)] for n in axes(sim_c, 2)]...)
            push!(momvals_c, round.(cov(sim_c_flat_0, sim_c_flat_lag, dims=2), digits=3))
            sim_d_flat_0 = hcat([sim_d[:, n, (1+lag):end] for n in axes(sim_d, 2)]...)
            sim_d_flat_lag = hcat([sim_d[:, n, 1:(end-lag)] for n in axes(sim_d, 2)]...)
            push!(momvals_d, round.(cov(sim_d_flat_0, sim_d_flat_lag, dims=2), digits=3))
            push!(rownames, "Covariance of y_t with y_t-" * string(lag))
        end

        push!(ds, DataFrame(:Moments => rownames, Symbol("Continuous Simulation") => momvals_c, Symbol("Discretized Simulation") => momvals_d))
    end

    return ds
end
