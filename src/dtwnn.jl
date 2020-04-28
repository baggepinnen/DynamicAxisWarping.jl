struct DTWWorkspace{T,AT<:AbstractArray,D}
    q::AT
    dist::D
    r::Int
    buffer::AT
    l::Vector{T}
    u::Vector{T}
    l_buff::Vector{T}
    u_buff::Vector{T}
    cb::Vector{T}
    c1::Vector{T}
    c2::Vector{T}
end

function DTWWorkspace(q::AbstractArray{QT}, dist, r::Int) where QT
    T      = floattype(QT)
    m      = lastlength(q)
    n      = 2r + 1
    buffer = similar(q)
    l      = zeros(T, m)
    u      = zeros(T, m)
    l_buff = zeros(T, m)
    u_buff = zeros(T, m)
    cb     = zeros(T, m)
    c1     = zeros(T, n)
    c2     = zeros(T, n)
    DTWWorkspace(q, dist, r, buffer, l, u, l_buff, u_buff, cb, c1, c2)
end

struct DTWSearchResult
    cost
    loc
    prunestats
    dists
end

function lower_upper_envs!#(w,buffer, bsf, query=false)
end

function lb_endpoints(w, buffer, best_so_far; kwargs...)
    dist = w.dist
    q = w.q
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

# function lb_env!(w::DTWWorkspace{T}, buffer, best_so_far; kwargs...) where T
#     lb = zero(T)
#     q, dist, u, l = w.q, w.dist, w.u, w.l
#     for i in eachindex(q)
#         x = buffer[i] # This function only supports data with natural ordering
#         if x > u[i]
#             d = dist(x, u[i]; kwargs...)
#         elseif x < l[i]
#             d = dist(x, l[i]; kwargs...)
#         end
#         lb += d
#         w.cb[i] = d
#         lb > best_so_far && return lb
#     end
#     return lb
# end
#
# function rev_cumsum!(cb)
#     for k = length(cb)-1:-1:1
#         cb[k] = cb[k+1] + cb[k]
#     end
# end

"""
    search_result = dtwnn(q, y, dist, rad; kwargs...)

Compute the nearest neighbor to `q` in `y`.

# Arguments:
- `q`: query (the short time series)
- `y`: data ( the long time series)
- `dist`: distance
- `rad`: radius
- `prune_endpoints = true`: use endopint heuristic
- `prune_envelope = false`: use envelope heuristic
- `bsf_multiplier = 1`: If > 1, require lower bound to exceed `bsf_multiplier*best_so_far`. Will only have effect if `saveall = false`.
- `saveall = false`: compute a dense result (takes longer, no early stopping methods used). If false, then a vector of lower bounds on the distance is stored in `search_result.dists`, if true, all distances are computed and stored.
"""
function dtwnn(q, y, dist, rad; kwargs...)
    w = DTWWorkspace(q, dist, rad)
    dtwnn(w, y; kwargs...)
end

function dtwnn(w::DTWWorkspace{T}, data::AbstractArray;
    prune_endpoints = true,
    prune_envelope = false,
    saveall = false,
    bsf_multiplier = 1,
    kwargs...) where T


    bsf_multiplier >= 1 || throw(ArgumentError("It does not make sense to have the bsf_multiplier < 1"))
    best_so_far = typemax(T)
    best_loc    = 1
    q           = w.q
    # TODO: normalize q
    m           = lastlength(q)
    md          = lastlength(data)
    md >= m || throw(ArgumentError("q must be shorter than y, swap inputs."))
    onedim      = ndims(q) == 1 && eltype(q) <: Real
    onedim && prune_envelope && lower_upper_envs!(w, q, best_so_far, true)

    # Counters to keep track of how many times lb helps
    prune_end   = 0
    prune_env   = 0
    dists = fill(typemax(T), md-m)

    @showprogress 1.5 "DTW NN" for it = 1:md-m
        bsf = bsf_multiplier*best_so_far
        buffer = data[!, (1:m) .+ (it-1)]
        # TODO: normalize
        # they copy d into circular array t at i%m and (i%m)+m, they then use t in place of buffer in many places to save time on renormalization
        if prune_endpoints && !saveall
            lb_end = lb_endpoints(w, buffer, bsf; kwargs...)
            if lb_end > bsf
                prune_end += 1
                continue
            end
        end
        if onedim && prune_envelope && !saveall # This bound only works when there is a natural ordering
            lower_upper_envs!(w, buffer, bsf)
            lb_env = lb_env!(w, buffer, bsf; kwargs...) # updates w.cb
            rev_cumsum!(w.cb)
            if lb_env > bsf
                prune_env += 1
                continue
            end
        end

        newdist = dtw_cost( buffer, q, w.dist, w.r;
            cumulative_bound = w.cb,
            best_so_far      = saveall ? typemax(T) : bsf,
            s1               = w.c1,
            s2               = w.c2,
            kwargs...
        )
        dists[it] = newdist
        if newdist < best_so_far
            best_so_far = newdist
            best_loc = it
        end
    end
    prunestats = (prune_end=prune_end, prune_env=prune_env)
    DTWSearchResult(best_so_far, best_loc, prunestats, dists)
end
