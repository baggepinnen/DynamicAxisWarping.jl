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
    T = QT <: Number ? QT : Float64
    m           = length(q)
    n           = 2r + 1
    buffer      = zeros(QT, m)
    l           = zeros(T, m)
    u           = zeros(T, m)
    l_buff      = zeros(T, m)
    u_buff      = zeros(T, m)
    cb          = zeros(T, m)
    c1          = zeros(T, n)
    c2          = zeros(T, n)
    DTWWorkspace(q, dist, r, buffer, l, u, l_buff, u_buff, cb, c1, c2)
end

struct DTWSearchResult
    cost
    loc
    prunestats
end

function lower_upper_envs!(w,buffer, bsf, query=false)

end

function lb_endpoints(w, buffer, best_so_far)
    dist = w.dist
    q = w.q
    m = lastlength(q)

    x1 = buffer[!,1]
    y1 = buffer[!,m]
    lb = dist(q[!,1], x1) + dist(q[!,m], y1)
    lb >= best_so_far && return lb

    x2 = buffer[!,2]
    d = min(dist(x2, q[!,1]), dist(x1,q[!,2]), dist(x2,q[!,2]))
    lb += d
    lb >= best_so_far && return lb

    y2 = buffer[!,m-1]
    d = min(dist(y2, q[!,m]), dist(y1,q[!,m-1]), dist(y2,q[!,m-1]))
    lb += d
    lb >= best_so_far && return lb

    return lb
    # TODO: can add more comparisons here
end

function lb_env!(w::DTWWorkspace{T}, buffer, best_so_far) where T
    lb = zero(T)
    q, dist, u, l = w.q, w.dist, w.u, w.l
    for i in eachindex(q)
        x = buffer[i] # This function only supports data with natural ordering
        if x > u[i]
            d = dist(x, u[i])
        elseif x < l[i]
            d = dist(x, l[i])
        end
        lb += d
        w.cb[i] = d
        lb > best_so_far && return lb
    end
    return lb
end

function rev_cumsum!(cb)
    for k = length(cb)-1:-1:1
        cb[k] = cb[k+1] + cb[k]
    end
end


function dtwnn(w::DTWWorkspace{T}, data::AbstractArray) where T

    best_so_far = typemax(T)
    best_loc    = 1
    q           = w.q
    # TODO: normalize q
    m           = lastlength(q)
    onedim      = ndims(q) == 1 && eltype(q) <: Real
    onedim && lower_upper_envs!(w, q, best_so_far, true)

    # Counters to keep track of how many times lb helps
    prune_end   = 0
    prune_env   = 0


    for it = 1:length(data)-m+1
        buffer = data[!, it:it+m-1]
        # TODO: normalize
        # they copy d into circular array t at i%m and (i%m)+m, they then use t in place of buffer in many places to save time on renormalization
        lb_end = lb_endpoints(w, buffer, best_so_far)
        if lb_end > best_so_far
            prune_end += 1
            continue
        end
        if onedim && false # This bound only works when there is a natural ordering
            lower_upper_envs!(w, buffer, best_so_far)
            lb_env = lb_env!(w, buffer, best_so_far) # updates w.cb
            rev_cumsum!(w.cb)
            if lb_env > best_so_far
                prune_env += 1
                continue
            end
        end

        newdist = dtw_cost( buffer, q, w.dist, w.r,
            cumulative_bound = w.cb,
            best_so_far      = best_so_far,
            cost             = w.c1,
            cost_prev        = w.c2,
        )
        if newdist < best_so_far
            best_so_far = newdist
            best_loc = it

        end
    end
    prunestats = (prune_end=prune_end, prune_env=prune_env)
    DTWSearchResult(best_so_far, best_loc, prunestats)
end
