struct DTWWorkspace{T,AT<:AbstractArray,D,N}
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
    normalizer::N
end

function DTWWorkspace(q::AbstractArray{QT}, dist, r::Int, normalizer=Nothing) where QT
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
    DTWWorkspace(q, dist, r, buffer, l, u, l_buff, u_buff, cb, c1, c2, normalizer)
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
- `prune_endpoints = true`: use endpoint heuristic
- `prune_envelope = false`: use envelope heuristic
- `bsf_multiplier = 1`: If > 1, require lower bound to exceed `bsf_multiplier*best_so_far`. Will only have effect if `saveall = false`.
- `saveall = false`: compute a dense result (takes longer, no early stopping methods used). If false, then a vector of lower bounds on the distance is stored in `search_result.dists`, if true, all distances are computed and stored.
"""
function dtwnn(q, y, dist, rad; normalizer=Val(Nothing), kwargs...)
    n = normalizer isa Val ? normalizer : Val(normalizer)
    q, y = setup_normalizer(n, q, y)
    w = DTWWorkspace(q, dist, rad, n)
    dtwnn(w, y; kwargs...)
end

function dtwnn(w::DTWWorkspace{T}, y::AbstractArray;
    prune_endpoints = true,
    prune_envelope  = false,
    saveall         = false,
    bsf_multiplier  = 1,
    kwargs...) where T


    bsf_multiplier >= 1 || throw(DomainError("It does not make sense to have the bsf_multiplier < 1"))
    best_so_far = typemax(T)
    best_loc    = 1
    q           = w.q
    m           = lastlength(q)
    my          = lastlength(y)
    my >= m || throw(ArgumentError("q must be shorter than y, swap inputs."))
    onedim      = ndims(q) == 1 && eltype(q) <: Real
    onedim && prune_envelope && lower_upper_envs!(w, q, best_so_far, true) # Result stored in w

    # Counters to keep track of how many times lb helps
    prune_end   = 0
    prune_env   = 0
    dists = fill(typemax(T), my-m)

    prog = Progress((my-m)÷10, 1)
    # @inbounds @showprogress 1.5 "DTW NN" for it = 1:my-m
    @inbounds for it = 1:my-m
        advance!(y)
        bsf = bsf_multiplier*best_so_far
        ym = y[!, (1:m) .+ (it-1)] # if y isa Normalizer, this is a noop
        if prune_endpoints && !saveall
            lb_end = lb_endpoints(w, ym, bsf; kwargs...)
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

        newdist = dtw_cost( buffern, q, w.dist, w.r;
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
        it % 10 == 0 && next!(prog)
    end
    prunestats = (prune_end=prune_end, prune_env=prune_env)
    DTWSearchResult(best_so_far, best_loc, prunestats, dists)
end


# q is normalized first   ( we do this in wrapper that creates workspace)
# then create envelope of q
# sort q
# store sorted q and u/l in qo,ql,qu
# read in buffer
# create envelop of buffer ( no normalization done before or inside here, but this step is only required if one does the reverse env bound as well)
# reset norm accumulators
# for i
#  advance accumulator with buffer[i] ( in effect, populate it fully before doing anything)
#  update circular array
# now they have an if statement that is only enterd when i >= m, inside that
#    mean,std = ...
#    calc lbs(mean, std)
#    if no pruning succeeded, then fully z-normalize t ( how does t relate to our z.x? length(t) == m
#    what they send into dtw is the fully znormalized t ( collect!(buffer, z) perhaps?)
# advance accumulators


# I wonder how much better it is to abandon the Z-normalization rather than using SIMD to normalize the entire thing in one go. If one has calculated the entire envelop, then the entire t has already been normalized, why redo it?
