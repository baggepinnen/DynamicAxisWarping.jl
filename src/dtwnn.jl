struct DTWWorkspace{T,AT<:AbstractArray,D,N}
    q::AT
    dist::D
    r::Int
    l::Vector{T}
    u::Vector{T}
    l_buff::Vector{T}
    u_buff::Vector{T}
    cb::Vector{T}
    c1::Vector{T}
    c2::Vector{T}
    normalizer::N
    function DTWWorkspace(q::AbstractArray{QT}, dist, r::Int, normalizer=Val(Nothing)) where QT
        T      = floattype(QT)
        m      = lastlength(q)
        n      = 2r + 1
        l      = zeros(T, m)
        u      = zeros(T, m)
        l_buff = zeros(T, m)
        u_buff = zeros(T, m)
        cb     = zeros(T, m)
        c1     = zeros(T, n)
        c2     = zeros(T, n)
        new{T, typeof(q), typeof(dist), typeof(normalizer)}(q, dist, r, l, u, l_buff, u_buff, cb, c1, c2, normalizer)
    end
end


struct DTWSearchResult{QT,C,D} <: AbstractSearchResult{QT}
    q::QT
    cost::C
    loc::Int
    prunestats
    dists::D
end

SlidingDistancesBase.value(r::DTWSearchResult) = r.cost
SlidingDistancesBase.location(r::DTWSearchResult) = r.loc
SlidingDistancesBase.payload(r::DTWSearchResult) = r.dists
SlidingDistancesBase.target(r::DTWSearchResult) = r.q

Base.findmin(results::Vector{<:DTWSearchResult}) = (i=argmin(results); (results[i].cost,i))
Base.findmax(results::Vector{<:DTWSearchResult}) = _findres(results, >)
Base.minimum(results::Vector{<:DTWSearchResult}) = findmin(results)[1]
Base.maximum(results::Vector{<:DTWSearchResult}) = findmax(results)[1]

function _findres(results::Vector{<:DTWSearchResult}, comp)
    mapreduce((a,b)->comp(a[1], b[1]) ? a : b, enumerate(results)) do (i,r)
        minimum(r.dists), i
    end
end

function lower_upper_envs!(w::DTWWorkspace{T}, q, bsf, query = false) where {T}
    du, dl = Deque{Int}(), Deque{Int}()
    push!(du, 0)
    push!(dl, 0)
    r = w.r
    if query
        u,l = w.u, w.l
    else
        u,l = w.u_buff, w.l_buff
    end
    m = lastlength(w.q)
    for i = 1:m-1
        if i > r
            u[i-r] = q[first(du)+1]
            l[i-r] = q[first(dl)+1]
        end
        if q[i+1] > q[i]
            pop!(du)
            while (!isempty(du) && q[i+1] > q[last(du)+1])
                pop!(du)
            end
        else
            pop!(dl)
            while (!isempty(dl) && q[i+1] < q[last(dl)+1])
                pop!(dl)
            end
        end
        push!(du, i)
        push!(dl, i)
        if i == 2r + 1 + first(du)
            popfirst!(du)
        elseif (i == 2r + 1 + first(dl))
            popfirst!(dl)
        end
    end
    for i = m:m+r
        u[i-r] = q[first(du)+1]
        l[i-r] = q[first(dl)+1]
        if i - first(du) >= 2r + 1
            popfirst!(du)
        end
        if i - first(dl) >= 2r + 1
            popfirst!(dl)
        end
    end
end

function lb_endpoints(dist, q, buffer, best_so_far; kwargs...)
    m = lastlength(q)

    x1 = buffer[!,1]
    y1 = buffer[!,m]
    lb = dist(q[!,1], x1; kwargs...) + dist(q[!,m], y1; kwargs...)
    lb >= best_so_far && return lb

    x2 = buffer[!,2]
    d = min(dist(x2, q[!,1]; kwargs...), dist(x1,q[!,2]; kwargs...), dist(x2,q[!,2]); kwargs...)
    lb += d
    lb >= best_so_far && return lb

    y2 = buffer[!,m-1]
    d = min(dist(y2, q[!,m]; kwargs...), dist(y1,q[!,m-1]; kwargs...), dist(y2,q[!,m-1]); kwargs...)
    lb += d
    lb >= best_so_far && return lb

    return lb
    # TODO: can add more comparisons here
end

function lb_env!(w::DTWWorkspace{T}, buffer, best_so_far; kwargs...) where T
    lb = zero(T)
    q, dist, u, l = w.q, w.dist, w.u, w.l
    for i in 1:lastlength(q)
        x = buffer[!,i] # This function only supports data with natural ordering
        d = zero(T)
        if x > u[i]
            d = dist(x, u[i]; kwargs...)
        elseif x < l[i]
            d = dist(x, l[i]; kwargs...)
        end
        lb += d
        w.cb[i] = d
        lb > best_so_far && return lb
    end
    return lb
end

function rev_cumsum!(cb)
    @inbounds for k = length(cb)-1:-1:1
        cb[k] = cb[k+1] + cb[k]
    end
end

"""
    search_result = dtwnn(q, y, dist, rad; kwargs...)

Compute the nearest neighbor to `q` in `y`.

# Arguments:
- `q`: query (the short time series)
- `y`: data ( the long time series)
- `dist`: distance
- `rad`: radius
- `showprogress`: Defaults to true
- `prune_endpoints = true`: use endpoint heuristic
- `prune_envelope  = true`: use envelope heuristic
- `bsf_multiplier  = 1`: If > 1, require lower bound to exceed `bsf_multiplier*best_so_far`.
- `saveall = false`: compute a dense result (takes longer, no early stopping methods used). If false, then a vector of lower bounds on the distance is stored in `search_result.dists`, if true, all distances are computed and stored.
- `avoid`: If an integer index (or set of indices) is provided, this index will be avoided in the search. This is useful in case `q` is a part of `y`.
"""
function dtwnn(q, y, dist, rad; normalizer=Val(Nothing), kwargs...)
    n = normalizer isa Val ? normalizer : Val(normalizer)
    n, q, y = setup_normalizer(n, q, y)
    w = DTWWorkspace(q, dist, rad, n)
    dtwnn(w, y; kwargs...)
end

function dtwnn(w::DTWWorkspace{T}, y::AbstractArray;
    prune_endpoints = true,
    prune_envelope  = true,
    saveall         = false,
    bsf_multiplier  = 1,
    transportcost   = 1,
    showprogress    = true,
    avoid           = nothing,
    kwargs...) where T


    w.normalizer !== Val(Nothing) && !isa(y, AbstractNormalizer) && @warn("Normalizer in use but `y` is not wrapped in a normalizer object. This will result in highly suboptimal performance", maxlog=10)
    bsf_multiplier >= 1 || throw(DomainError("It does not make sense to have the bsf_multiplier < 1"))
    best_so_far = typemax(T)
    best_loc    = 1
    q           = w.q
    m           = lastlength(q)
    my          = actuallastlength(y)
    my >= m || throw(ArgumentError("q must be shorter than y, swap inputs."))
    onedim      = ndims(q) == 1 && eltype(q) <: Real
    onedim && prune_envelope && lower_upper_envs!(w, q, best_so_far, true) # Result stored in w

    # Counters to keep track of how many times lb helps
    prune_end   = 0
    prune_env   = 0
    dists       = fill(typemax(T), my-m+1)

    prog = Progress((my-m)÷max(my÷100,1), dt=1, desc="DTW NN")
    # @inbounds @showprogress 1.5 "DTW NN" for it = 1:my-m
    @inbounds for it = 1:my-m+1
        showprogress && it % max(my÷100,1) == 0 && next!(prog)
        advance!(y)
        avoid !== nothing && it ∈ avoid && continue
        bsf = bsf_multiplier*best_so_far
        ym = getwindow(y, m, it)
        if prune_endpoints && !saveall
            lb_end = lb_endpoints(w.dist, w.q, ym, bsf; kwargs...)
            if lb_end > bsf
                prune_end += 1
                continue
            end
        end
        if onedim && prune_envelope && !saveall # This bound only works when there is a natural ordering
            # lower_upper_envs!(w, ym, bsf) # This step is only required for reverse bound
            lb_env = lb_env!(w, ym, bsf; kwargs...) # updates w.cb
            rev_cumsum!(w.cb)
            if lb_env > bsf
                prune_env += 1
                continue
            end
        end
        # If we get here, we must normalize the entire y
        buffern = normalize(w.normalizer, ym) # This only normalizes what's not already normalized

        newdist = dtw_cost(q, buffern, w.dist, w.r;
            cumulative_bound = w.cb,
            best_so_far      = saveall ? typemax(T) : bsf,
            s1               = w.c1,
            s2               = w.c2,
            transportcost    = transportcost,
            kwargs...
        )
        dists[it] = newdist
        if newdist < best_so_far
            best_so_far = newdist
            best_loc = it
        end
    end
    prunestats = (prune_end=prune_end, prune_env=prune_env)
    DTWSearchResult(q, best_so_far, best_loc, prunestats, dists)
end

struct Neighbor{T}
    i::Int
    d::T
end

Base.isless(n1::Neighbor,n2::Neighbor) = isless(n1.d, n2.d)
Base.isless(n1, n2::Neighbor) = isless(n1, n2.d)
Base.isless(n1::Neighbor,n2) = isless(n1.d, n2)

"""
    dists, inds = sparse_distmat(y::Vector{<:AbstractVector{T}}, k, dist, radius; kwargs...) where T

Compute the `k` nearest neighbors between signals in `y`, corresponding to the `k` smallest entries in each row of the pairwise distance matrix. The return values are vectors of length-k vectors with the calculated distances and neighbor indices.

#Arguments:
- `y`: Vector of vectors containing the signals
- `k`: number of neighbors
- `dist`: the inner metric, e.g., `SqEuclidean()`
- `kwargs`: these are sent to `dtw_cost`.
"""
function sparse_distmat(y::Vector{<:AbstractVector{S}}, k, dist, rad; kwargs...) where S
    T = floattype(S)
    N = length(y)
    INDS = [zeros(Int, k) for _ in 1:N]
    DISTS = [zeros(T, k) for _ in 1:N]
    for i = 1:N
        bsf = typemax(T)
        dists = BinaryMaxHeap{Neighbor{T}}()
        for j = 1:N
            j == i && continue
            d = lb_endpoints(dist, y[i], y[j], bsf; kwargs...)
            if d < bsf
                d = dtw_cost(y[i], y[j], dist, rad; best_so_far = bsf, kwargs...)
            end
            push!(dists, Neighbor(j,T(d)))
            if length(dists) > k
                bsf = pop!(dists)
            end
        end

        for j = k:-1:1
            n = pop!(dists)
            INDS[i][j] = n.i
            DISTS[i][j] = n.d
        end
    end
    DISTS, INDS
end
