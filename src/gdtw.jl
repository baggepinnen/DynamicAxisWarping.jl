# A cache to preallocate everything for the GDTW distance
struct GDTWWorkspace{T1, T2, T3}
    τ::T2
    l::T1
    u::T1
    l_prev::T1
    u_prev::T1
    l₀::T1
    u₀::T1
    min_costs::T2
    costs::T3
end

"""
    GDTWWorkspace{T}(M, N)

Creates a cache of numeric type `T` for use in [`gdtw`](@ref).
"""
function GDTWWorkspace(::Type{T}, M, N) where {T}
    GDTWWorkspace{Vector{T}, Matrix{T}, Array{T, 3}}(
        zeros(T, M, N), zeros(T, N), zeros(T, N),
        zeros(T, N), zeros(T, N), zeros(T, N),
        zeros(T, N), zeros(T, M, N), zeros(T, M, M, N)
    )
end

GDTWWorkspace(M, N) = GDTWWorkspace(Float64, M, N)

# refine the bounds, as described in Section 4.2 of DB19
function refine!(l_current, u_current, l_prev, u_prev, l₀, u₀, warp; η)
    @avx for i in eachindex(l_current, u_current, l_prev, u_prev, l₀, u₀, warp)
        δ = η * (u_prev[i] - l_prev[i]) / 2
        l_current[i] = max(warp[i] - δ, l₀[i])
        u_current[i] = min(warp[i] + δ, u₀[i])
    end
    return nothing
end

# inital choices of `l` and `u`, modified from Eq (7) of DB19
function inital_bounds!(l, u, t, smin, smax, symmetric)
    # We need to loosen the bounds to account for floating point error
    # Otherwise these bounds can be too tight and disallow valid moves
    # which can lead to very wrong results.
    smin = .99*smin
    smax = 1.01*smax

    @inbounds for i in eachindex(t, l, u)
        # You must be able to get to `warp[i]` in time `t[i]`, so
        #   `smax >= warp[i] / t[i] >= smin`
        # This gives a lower and upper bound on `warp[i]`, i.e., on `τ[:, i]`.
        # You must be able to get to `1` from `warp[i]` in time `1-t[i]`, so
        #   `smin <= (1-warp[i])/(1-t[i]) <= smax`
        # this gives another lower and upper bound.
        # The resulting bounds:
        lower = max(smin * t[i], 1 - smax * (1 - t[i]))
        upper = min(smax * t[i], 1 - smin * (1 - t[i]))

        if symmetric
            # We need to apply the above bounds to `ψ(s) = 2s - ϕ(s) = 2t[i] - warp[i]`
            # as well Which leads to the following.
            l[i] = max(lower, 2*t[i] - upper)
            u[i] = min(upper, 2*t[i] - lower)
        else
            l[i] = lower
            u[i] = upper
        end
    end
    return nothing
end

# an inplace update to `τ` with the bounds `l` and `u`.
# produces `τ` as described in Section 4.1 of DB19.
function update_τ!(τ, t, M, l, u)
    N = length(t)
    @assert size(τ) == (M, N)
    @inbounds for t = 1:N, j = 1:M
        τ[j, t] = l[t] + ((j - 1) / (M - 1)) * (u[t] - l[t])
    end
    return nothing
end

"""
    gdtw(
        x,
        y,
        ::Type{T}  = Float64;
        symmetric::Bool = true,
        M::Int     = 100,
        N::Int     = 100,
        t          = range(T(0), stop = T(1), length = N),
        cache::GDTWWorkspace = GDTWWorkspace(T, M, length(t)),
        λcum       = T(0.01),
        λinst      = T(0.01),
        η          = T(1 / 8),
        max_iters  = 3,
        metric     = (x, y) -> norm(x - y),
        Rcum       = abs2,
        smin::Real = T(0.001),
        smax::Real = T(5.0),
        Rinst      = symmetric  ?
                        ϕ′ -> ( (smin <= ϕ′ <= smax)
                            && (smin <= 2 - ϕ′ <= smax) ) ? (ϕ′-1)^2 : typemax(T)
                                :
                        ϕ′ -> (smin <= ϕ′ <= smax) ? (ϕ′-1)^2 : typemax(T),
        verbose    = false,
        warp       = zeros(T, length(t)),
        callback   = nothing,
    ) where T -> cost, ϕ, ψ

Computes a general DTW distance following [DB19](https://arxiv.org/abs/1905.12893).

Aims to find `ϕ(s)` to minimize

    ∫ metric(x(ϕ(s)), y(ψ(s))) + λinst*Rinst(ϕ'(s) - 1) + λcum*Rcum(ϕ(s) - s) ds

over the interval `s ∈ [0,1]`, where `ψ(s) = 2s - ϕ(s)` (if `symmetric=true`) or `ψ(s) = s`
(if `symmetric = false`). The integral is discretized in time into `N` points (or according
to the times `t`, if `t` is specified). Additionally, the possible values obtained by `ϕ`
(and hence `ψ`) at each possible time `s` are discretized into `M` points.

If `max_iters > 1`, then after solving the doubly-discretized problem to obtain the optimal `ϕ`,
the problem is solved again by choosing a new discretization of `M` possible values
of `ϕ(s)` in an interval (whose width is governed by the parameter `η`) around the
previous optimal value. This is repeated until the problem has been solved `max_iters`
times in total. Setting `verbose=true` prints the cost at each iteration; a "high enough"
value of `max_iters` can be chosen by inspecting when the cost stabilizes sufficiently.

The parameters are:

* `x`: the continuous time signal to warp (see [`LinearInterpolation`](@ref) for generating such a signal from discrete data)
* `y`: the continuous-time signal to warp to
* `T`: the numeric type to be used in the problem
* `symmetric`: if true, `ψ(s) = 2s - ϕ(s)`, otherwise `ψ(s) = s`.
* `t`: the discretization of time on `[0,1]`; either `t` or `N` should be specified
* `M`: the discretization of the values obtained by the warping path
* `metric`:  a function `metric(u,v) -> ℝ` to compute differences between the signals at a time point (such as a Distances.jl distance)
* `Rcum`: penalty function on the cumulative warp
* `Rinst`: penalty function on the instantaenous warping. Should be infinite outside of `[smin, smax]`.
* `smin`, `smax`: minimum and maximum allowed instantaenous warping. Should have `smin > 0` and `smin < smax`.
* `λcum`, `λinst`: the regularization constants for `Rcum` and `Rinst`, respectively

The following may be pre-allocated and reused between distance computations with the same `M` and `N` (or `length(t)`).

* `cache`: a cache of matrices and vectors, generated by `GDTW.GDTWWorkspace{T}(N,M)`

"""
function gdtw(args...; kwargs...)
    data = prepare_gdtw(args...; kwargs...)
    cost = iterative_gdtw!(data)
    return cost, gdtw_warpings(data)...
end

"""
    prepare_gdtw(x, y; kwargs...)

Creates a NamedTuple of parameters, using the same keyword argments as `dist`.
A preprocessing step before calling `iterative_gdtw!`.
"""
function prepare_gdtw(
    x,
    y,
    ::Type{T}  = Float64;
    symmetric::Bool = true,
    M::Int     = 100,
    N::Int     = 100,
    t::AbstractVector{T} = range(T(0), stop = T(1), length = N),
    cache::GDTWWorkspace = GDTWWorkspace(T, M, length(t)),
    λcum::T    = T(0.01),
    λinst::T   = T(0.01),
    η::T       = T(1 / 8),
    max_iters  = 3,
    metric     = (x, y) -> norm(x - y),
    Rcum       = abs2,
    smin::T    = T(0.001),
    smax::T    = T(5.0),
    Rinst      = symmetric  ?
                    ϕ′ -> ( (smin <= ϕ′ <= smax)
                        && (smin <= 2 - ϕ′ <= smax) ) ? (ϕ′-1)^2 : typemax(T)
                            :
                    ϕ′ -> (smin <= ϕ′ <= smax) ? (ϕ′-1)^2 : typemax(T),
    verbose::Bool = false,
    warp       = zeros(T, length(t)),
    callback   = nothing,
) where T
    N = length(t)

    (M > N / smax) || @warn "`M <= N / smax`; problem may be infeasible" M N smax


    @unpack l₀, l_prev, l,  u₀, u_prev, u, τ = cache
    inital_bounds!(l₀, u₀, t, smin, smax, symmetric)
    l_prev .= l₀
    u_prev .= u₀
    u .= u₀
    l .= l₀
    update_τ!(τ, t, M, l, u)

    function node_weight(j, s)
        s == length(t) && return zero(T)
        Rval = Rcum(τ[j, s] - t[s])
        yval = symmetric ? 2*t[s] - τ[j, s] : t[s]
        (t[s+1] - t[s])*(metric(x(τ[j, s]), y(yval)) + λcum * Rval)
    end

    @inline function edge_weight((j, s), (k, s2))
        s + 1 ≠ s2 && return typemax(T)
        ϕ′ = (τ[k, s+1] - τ[j, s]) / (t[s+1] - t[s])
        (t[s+1] - t[s]) * (λinst * Rinst(ϕ′))
    end

    return (
        iter        = Ref(1),
        N           = N,
        M           = M,
        τ           = τ,
        node_weight = node_weight,
        edge_weight = edge_weight,
        η           = η,
        max_iters   = max_iters,
        t           = t,
        smin        = smin,
        smax        = smax,
        callback    = callback,
        verbose     = verbose,
        metric      = metric,
        cache       = cache,
        warp        = warp,
        symmetric   = symmetric,
    )
end

"""
    iterative_gdtw!(data; max_iters = data.max_iters, verbose = data.verbose) -> cost

Runs the GDTW algorithm with iterative refinement for `max_iters` iterations,
returning the resulting cost. Uses the [`GDTWWorkspace`](@ref) in `data.cache` and
updates the `data.warp` vector. Here, `data` is usually obtained by [`prepare_gdtw`](@ref).

This can be called multiple times on the same `data` with a higher value of `max_iters`
to refine a calculation without starting over, which can be useful for checking convergence.

## Example

```julia
data = prepare_gdtw(x,y; max_iters = 3)
cost3 = iterative_gdtw!(data) # performs 3 iterations, returns the cost
cost10 = iterative_gdtw!(data; max_iters = 10) # performs 7 more iterations, returns the cost
ϕ, ψ = gdtw_warpings(data)

# Equivalently,
cost10, ϕ, ψ = gdtw(x,y; max_iters = 10)
```
"""
function iterative_gdtw!(data; max_iters = data.max_iters, verbose = data.verbose)
    @unpack N, M, τ, η, iter, t, callback = data
    @unpack cache, warp, symmetric = data
    if iter[] > max_iters
        @warn "`iter[] > max_iters`; no iterations performed." iter[] max_iters
        return zero(eltype(warp))
    end
    local cost

    # First iteration is special cased because
    # we can't send the cost and warp to the callback yet,
    # and we can quit early if we don't need to do refinement.
    if iter[] == 1
        if callback !== nothing
            callback((iter=1, t=t, τ=τ))
        end

        cost = single_gdtw!(data)
        verbose && @info "Iteration" iter[] cost
        iter[] += 1
        max_iters == 1 && return cost
    end

    @unpack l_prev, u_prev, l, u, l₀, u₀ = cache

    while iter[] <= max_iters
        l_prev .= l
        u_prev .= u
        refine!(l, u, l_prev, u_prev, l₀, u₀, warp; η=η)
        update_τ!(τ, t, M, l, u)
        cost = single_gdtw!(data)
        if callback !== nothing
            callback((iter = iter, t = t, τ = τ, warp = warp, cost = cost))
        end
        verbose && @info "Iteration" iter[] cost

        iter[] += 1
    end

    return cost
end

"""
    gdtw_warpings(data) -> ϕ, ψ

Computes the interpolations from a `data` `NamedTuple`
with entries for the time points `t`, warping points `warp`,
and a boolean `symmetric`.
"""
function gdtw_warpings(data)
    @unpack t, warp, symmetric = data
    if symmetric
        ψ = LinearInterpolation(2*t - warp, t)
    else
        ψ = LinearInterpolation(t, t)
    end
    ϕ = LinearInterpolation(warp, t)
    return ϕ, ψ
end

## Dynamic programming to compute the distance

function single_gdtw!(data::T) where {T}
    @unpack N, M, node_weight, cache, edge_weight, τ, warp = data
    @unpack min_costs, costs = cache
    calc_costs!(min_costs, costs, N, M, node_weight, edge_weight)
    cost = min_costs[end, end]
    trackback!(warp, costs, τ)
    return cost
end

function calc_costs!(min_costs, costs, N, M, node_weight::F1, edge_weight::F2) where {F1,F2}
    @boundscheck checkbounds(min_costs, 1:M, 1:N)
    @boundscheck checkbounds(costs, 1:M, 1:M, 1:N)

    @inbounds begin
        min_costs .= node_weight.(1:M, permutedims(1:N))
        # t = 2 case
        for j = 1:M
            costs[1, j, 2] = min_costs[1, 1] + edge_weight((1, 1), (j, 2))
            min_costs[j, 2] += costs[1, j, 2]
        end
        for t = 3:N
            for j = 1:M
                mi = typemax(eltype(costs))
                for k = 1:M
                    c = min_costs[k, t-1] + edge_weight((k, t - 1), (j, t))
                    costs[k, j, t] = c
                    mi = ifelse(c < mi, c, mi)
                end
                min_costs[j, t] += mi
            end
        end
    end
    return nothing
end

function trackback!(warp, costs, τ)
    (M, _, N) = size(costs)
    @boundscheck checkbounds(costs, 1:M, 1:M, 1:N)
    c = M
    @inbounds for t = N:-1:3
        warp[t] = τ[c, t]
        c = argmin(@views costs[:, c, t])
    end
    warp[2] = τ[c, 2]
    warp[1] = τ[1, 1]

    return nothing
end

## Interpolations


"""
LinearInterpolation(x::AbstractVector) -> Function

Provides a linear interpolation of `x` on the interval `[0,1]`.
"""
struct LinearInterpolation{Tx,Tt} <: Function
    x::Tx
    t::Tt
    function LinearInterpolation(x::Tx, ts::Ts) where {Tx,Ts}
        issorted(ts) || throw(ArgumentError("Time parameter `ts` must be sorted in increasing order."))
        T = eltype(ts)
        t = (ts .- T(first(ts))) ./ T( last(ts) - first(ts))
        Tt = typeof(t)
        new{Tx,Tt}(x, t)
    end
end

LinearInterpolation(x) = LinearInterpolation(x, axes(x, ndims(x)))

function (xt::LinearInterpolation)(s)
    x = xt.x
    0 <= s <= 1 || return zero(x[!, 1])
    t = xt.t
    i = searchsortedlast(t, s)
    (i == 0) && return x[!, 1]
    (i == lastlength(x)) && return x[!, lastlength(x)]
    (s == t[i]) && return x[!, i]
    weight = (s - t[i]) / (t[i+1] - t[i])
    omw = 1 - weight
    x[!, i] .* omw .+ x[!, i+1] .* weight
end
