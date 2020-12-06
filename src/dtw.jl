#####################################
#     Basic interface functions     #
#####################################

"""
    cost,i1,i2 = dtw(seq1, seq2, [dist=SqEuclidean])
    cost,i1,i2 = dtw(seq1, seq2, dist, i2min, i2max)

Perform dynamic-time warping to measure the distance between two sequences.

Find a set of indices (`i1`,`i2`) that align two series (`seq1`,`seq2`) by
dynamic axis warping. Also returns the distance (after warping) according to
the SemiMetric `dist`, which defaults to squared Euclidean distance (see
Distances.jl). If `seq1` and `seq2` are matrices, each column is considered
an observation.

If `i2min/max` are provided, do DTW to align `seq1` and `seq2` confined to a window. Vectors `i2min` and
`i2max` specify (inclusive) lower and upper bounds for `seq2` for each index in
`seq1`. Thus, `i2min` and `i2max` are required to be the same length as `seq1`.

See also [`dtw_cost`](@ref) and [`dtwnn`](@ref).
"""
function dtw(args...; kwargs...)
    D = dtw_cost_matrix(args...; kwargs...)
    return trackback(D)
end

##############################
#  Cost matrix computations  #
##############################

Distances.pairwise(d::PreMetric, s1::AbstractVector, s2::AbstractVector; dims=2) = evaluate.(Ref(d), s1, transpose(s2))

function Distances.pairwise(d::PreMetric, s1::AbstractArray, s2::AbstractArray; dims=2)
    [evaluate(d, s1[!,i], s2[!,j]) for i in 1:lastlength(s1), j in lastlength(s2)]
end

@inbounds function dtw_cost_matrix(seq1::AbstractArray{T}, seq2::AbstractArray{T}, dist::SemiMetric = SqEuclidean();
    transportcost=1) where T
    # Build the cost matrix
    m = lastlength(seq2)
    n = lastlength(seq1)

    # Initialize first column and first row
    D = pairwise(dist, seq2, seq1, dims=2)
    @assert size(D) == (m,n)

    for r=2:m
        D[r,1] += D[r-1,1]
    end
    for c=2:n
        D[1,c] += D[1,c-1]
    end

    # Complete the cost matrix
    for c = 2:n
        for r = 2:m
            best_neighbor_cost = min(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1])
            D[r, c] += best_neighbor_cost
        end
    end

    return D
end

Base.@propagate_inbounds function dtw_cost_matrix(
    seq1::AbstractArray{T},
    seq2::AbstractArray{T},
    dist::SemiMetric,
    i2min::AbstractVector{U},
    i2max::AbstractVector{U};
    transportcost = 1
) where {T,U<:Integer}
    n = lastlength(seq1) # of columns in cost matrix
    m = lastlength(seq2) # of rows in cost matrix
    Base.@boundscheck begin
        n == length(i2min) || throw(ArgumentError("i2min does not match length of seq1."))
        n == length(i2max) || throw(ArgumentError("i2max does not match length of seq1."))
        1 == i2min[1]      || throw(ArgumentError("i2min must start at 1."))
        m == i2max[end]    || throw(ArgumentError("i2max must end at length(seq2), was $(i2max[end]) ≂̸ $(m)"))
    end

    # Build the (n x m) cost matrix into a WindowedMatrix, because it's ragged.
    # That type gives efficient storage with convenient [r,c] indexing and returns
    # Inf when accessed outside the window.
    D = WindowedMatrix(i2min, i2max, Inf)

    # First column first
    D[1, 1] = evaluate(dist, seq1[!, 1], seq2[!, 1])
    for r = 2:i2max[1]
        D[r, 1] = D[r-1, 1] + evaluate(dist, seq1[!, 1], seq2[!, r])
    end

    # Complete the cost matrix from columns 2 to m.
    for c = 2:n
        for r = i2min[c]:i2max[c]
            best_neighbor_cost = min(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1])
            D[r, c] = best_neighbor_cost + evaluate(dist, seq1[!, c], seq2[!, r])
        end
    end

    return D
end

########################################
#  Find Best Path through Cost Matrix  #
########################################

"""
    cost,cols,rows = trackback(D::Matrix)

Given the cost matrix `D`, computes the optimal track from end to beginning.
Returns `cols` and `rows` which are vectors respectively holding the track.
"""
function trackback(D::AbstractMatrix{T}) where {T<:Number}

    # initialize trackback throught rows/columns
    r, c       = size(D)
    rows, cols = Int[r], Int[c]

    # estimate that we'll need N⋅logN elements
    N  = max(r, c)
    sz = 2 * N
    sizehint!(rows, sz)
    sizehint!(cols, sz)

    # do trackback
    @inbounds while r > 1 && c > 1
        tb, r, c = indmin3(D[r-1, c-1], D[r-1, c], D[r, c-1], r, c)
        push!(rows, r)
        push!(cols, c)
    end
    # Possibly either r>1 or c>1 at this point (but not both).
    # Add the unfinished part of the track to reach [1,1]
    for r = r-1:-1:1
        push!(rows, r)
        push!(cols, 1)
    end
    for c = c-1:-1:1
        push!(rows, 1)
        push!(cols, c)
    end
    return D[end, end], reverse(cols), reverse(rows)
end





"""
    dtw_cost(a::AbstractArray, b::AbstractArray, dist::Distances.SemiMetric, r::Int; best_so_far = Inf, cumulative_bound = Zeros(length(a)))

Perform dynamic time warping to measure the distance between two sequences.

Calculate the DTW cost between `a` and `b` with maximum warping radius `r`. You may provide values of `best_so_far` and `cumulative_bound` in order to enable early stopping.

# Keyword arguments:
- `best_so_far`: The best cost value obtained so far (optional)
- `cumulative_bound`: A vector the same length as a and b (optional)
- `s1`: Optional storage vector of length 2r+1, can be used to save allocations.
- `s2`: Optional storage vector of length 2r+1, can be used to save allocations.

Providing the two vectors `s1, s2` does not save very much time, but it makes the function completely allocation free. Can be useful in a threaded context.

See also [`dtw`](@ref) and [`dtwnn`](@ref).
This code was inspired by https://www.cs.ucr.edu/~eamonn/UCRsuite.html
"""
function dtw_cost(
    a::AbstractArray{QT},
    b::AbstractArray,
    dist::F,
    r::Int;
    transportcost = 1,
    best_so_far = typemax(floattype(QT)),
    cumulative_bound = Zeros(lastlength(a)),
    s1 = fill(typemax(floattype(QT)), 2r + 1),
    s2 = fill(typemax(floattype(QT)), 2r + 1),
    kwargs...
) where {QT, F <: Union{Distances.SemiMetric, Function}}

    T = floattype(QT)
    cost      = s1 # just change the name of these variables
    cost_prev = s2

    # Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    m, mb = lastlength(a), lastlength(b)
    mb == m || throw(ArgumentError("a and b must have the same length, got $m and $mb. To compare two series of different lengths, use function dtw"))
    length(cumulative_bound) == m || throw(ArgumentError("cumulative_bound and a must have the same length."))
    length(s1) == 2r+1 || throw(ArgumentError("s1 must be length 2r+1."))
    length(s2) == 2r+1 || throw(ArgumentError("s2 must be length 2r+1."))

    local k

    for i = 0:m-1
        k = max(0, r - i)
        min_cost = typemax(T)

        for j = max(0, i - r):min(m - 1, i + r)
            if j == 0 && i == 0
                cost[k+1] = dist(a[!,1], b[!,1]; kwargs...)
                min_cost = cost[k+1]
                k += 1
                continue
            end
            y = (j - 1 < 0) || (k - 1 < 0)     ? typemax(T) : cost[k]
            x = (i - 1 < 0) || (k + 1 > 2 * r) ? typemax(T) : transportcost*cost_prev[k+2]
            z = (i - 1 < 0) || (j - 1 < 0)     ? typemax(T) : transportcost*cost_prev[k+1]

            cost[k+1] = min(x, y, z) + dist(a[!,i+1], b[!,j+1]; kwargs...)

            # Find minimum cost in row for early stopping
            if cost[k+1] < min_cost
                min_cost = cost[k+1]
            end
            k += 1
        end

        # We can abandon early if the current cumulative distace with lower bound together are larger than best_so_far
        if ((i + r) < (m - 1)) && (min_cost + cumulative_bound[i+r+1] >= best_so_far)
            return min_cost + cumulative_bound[i+r+1]
        end

        cost_prev, cost = cost, cost_prev
    end

    # the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    final_dtw = cost_prev[k]
    return T(final_dtw)
end




@inbounds function soft_dtw_cost_matrix(seq1::AbstractArray, seq2::AbstractArray, dist::SemiMetric = SqEuclidean(); γ = 1,
    transportcost=1, radius=nothing)
    # Build the cost matrix
    m = lastlength(seq2)
    n = lastlength(seq1)

    # Initialize first column and first row
    D = pairwise(dist, seq2, seq1, dims=2)
    @assert size(D) == (m,n)

    for r=2:m
        D[r,1] += D[r-1,1]
    end
    for c=2:n
        D[1,c] += D[1,c-1]
    end

    # Complete the cost matrix
    if radius === nothing
        for c = 2:n
            for r = 2:m
                D[r, c] += softmin(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1], γ)
            end
        end
    else
        for c = 2:n
            for r = 2:m
                if abs(c-r) > radius
                    D[r, c] += 1/γ
                    # continue
                end    
                D[r, c] += softmin(transportcost*D[r-1, c], D[r-1, c-1], transportcost*D[r, c-1], γ)
            end
        end
    end

    return D
end



"""
    soft_dtw_cost(args...; γ = 1, kwargs...)

Perform Soft DTW. This is a differentiable version of DTW. The "distance" returned by this function is quite far from a true distance and can be negative. A smaller value of `γ` makes the distance closer to the standard DTW distance.


To differentiate w.r.t. the first argument, try
```julia
using ReverseDiff
da = ReverseDiff.gradient(a->soft_dtw_cost(a,b), a)
```

Ref: "Soft-DTW: a Differentiable Loss Function for Time-Series" https://arxiv.org/pdf/1703.01541.pdf

#Arguments:
- `args`: same as for [`dtw`](@ref)
- `γ`: The smoothing factor. A small value means less smoothing and a result closer to [`dtw_cost`](@ref)
- `kwargs`: same as for [`dtw`](@ref)
"""
function soft_dtw_cost(args...; γ = 1, kwargs...)
    D = soft_dtw_cost_matrix(args...; γ = γ, kwargs...)
    D[end,end]
end


@fastmath @inline function softmin(a, b, c, γ)
    γ = -γ
    a,b,c = a/γ, b/γ, c/γ
    maxv = max(a,b,c)
    ae,be,ce = exp(a - maxv), exp(b - maxv), exp(c - maxv)
    γ*(log(ae+be+ce) + maxv)
end

const LVB = LoopVectorization.VectorizationBase
@inline function softmin(a::T, b::T, c::T, γ) where T <: Union{Float64, Float32}
    γ = -γ
    ninvγ = one(T) / γ
    v = LVB.Vec{4,T}(a, b, c, typemax(T))
    v = v * ninvγ
    @fastmath maxv = min(a,b,c) * ninvγ
    ve = exp(v - maxv) * LVB.Vec{4,T}(one(T), one(T), one(T), zero(T))
    γ*(log(LVB.vsum(ve)) + maxv)
end
