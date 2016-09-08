#####################################
#     Basic interface functions     #
#####################################

"""
    cost,i1,i2 = dtw(seq1, seq2, [dist=SqEuclidean])

Find a set of indices (`i1`,`i2`) that align two time series (`seq1`,`seq2`) by
dynamic time warping. Also returns the distance (after warping) according to
the SemiMetric `dist`, which defaults to squared Euclidean distance (see
Distances.jl). If `seq1` and `seq2` are matrices, each column is considered
an observation.
"""
function dtw(
        seq1::Sequence,
        seq2::Sequence,
        dist::SemiMetric=SqEuclidean()
    )

    D = dtw_cost_matrix(seq1, seq2, dist)
    return trackback(D)
end

"""
    cost,i1,i2 = dtw(seq1,seq2,i2min,i2max,[dist=SqEuclidean])

Do DTW to align `seq1` and `seq2` confined to a window. Vectors `i2min` and
`i2max` specify (inclusive) lower and upper bounds for `seq2` for each index in
`seq1`. Thus, `i2min` and `i2max` are required to be the same length as `seq1`.
"""
function dtw(
        seq1::Sequence,
        seq2::Sequence,
        i2min::AbstractVector,
        i2max::AbstractVector,
        dist::SemiMetric = SqEuclidean()
    )

    D = dtw_cost_matrix(seq1, seq2, i2min, i2max, dist)
    return trackback(D)
end

##############################
#  Cost matrix computations  #
##############################
function dtw_cost_matrix{N,T}(
        seq1::Sequence{N,T},
        seq2::Sequence{N,T},
        dist::SemiMetric=SqEuclidean()
    )
    # Build the cost matrix
    m = length(seq2)
    n = length(seq1)
    D = zeros(T, m, n)

    # Initialize first column and first row
    D[1,1] = evaluate(dist, seq1[1], seq2[1])
    for r=2:m
        D[r,1] = D[r-1,1] + evaluate(dist, seq1[1], seq2[r])
    end
    for c=2:n
        D[1,c] = D[1,c-1] + evaluate(dist, seq1[c], seq2[1])
    end

    # Complete the cost matrix
    for c=2:n
        for r=2:m
            best_neighbor_cost = min(D[r-1,c], D[r-1,c-1], D[r,c-1])
            D[r,c] = best_neighbor_cost + evaluate(dist, seq1[c], seq2[r])
        end
    end

    return D
end

function dtw_cost_matrix{N,T,U<:Integer}(
        seq1::Sequence{N,T},
        seq2::Sequence{N,T},
        i2min::AbstractVector{U},
        i2max::AbstractVector{U},
        dist::SemiMetric = SqEuclidean()
    )
    m = length(seq2) # of rows in cost matrix
    n = length(seq1) # of columns in cost matrix
    n == length(i2min) || throw(ArgumentError("i2min does not match length of seq1."))
    n == length(i2max) || throw(ArgumentError("i2max does not match length of seq1."))
    1 == i2min[1] || throw(ArgumentError("i2min must start at 1."))
    m == i2max[end] || throw(ArgumentError("i2max must end at length(seq2)."))

    # Build the (n x m) cost matrix into a WindowedMatrix, because it's ragged.
    # That type gives efficient storage with convenient [r,c] indexing and returns
    # Inf when accessed outside the window.
    D = WindowedMatrix(i2min, i2max, Inf)

    # First column first
    D[1,1] = evaluate(dist, seq1[1], seq2[1])
    for r=2:i2max[1]
        D[r,1] = D[r-1,1]  + evaluate(dist, seq1[1], seq2[r])
    end

    # Complete the cost matrix from columns 2 to m.
    for c=2:n
        for r=i2min[c]:i2max[c]
            best_neighbor_cost = min(D[r-1,c], D[r-1,c-1], D[r,c-1])
            D[r,c] = best_neighbor_cost + evaluate(dist, seq1[c], seq2[r])
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
function trackback{T<:Number}(D::AbstractMatrix{T})

    # initialize trackback throught rows/columns
    r,c = size(D)
    rows,cols = Int[r],Int[c]

    # estimate that we'll need N⋅logN elements
    N = max(r,c)
    sz = 2*N
    sizehint!(rows,sz)
    sizehint!(cols,sz)

    # do trackback
    while r > 1 && c > 1
        tb,r,c = indmin3(D[r-1,c-1], D[r-1,c], D[r,c-1], r, c)
        push!(rows,r)
        push!(cols,c)
    end
    # Possibly either r>1 or c>1 at this point (but not both). 
    # Add the unfinished part of the track to reach [1,1]
    for r=r-1:-1:1
        push!(rows,r)
        push!(cols,1)
    end
    for c=c-1:-1:1
        push!(rows,1)
        push!(cols,c)
    end
    return D[end,end], reverse(cols), reverse(rows)
end


###############################################
# Wrapper functions for abstract array inputs #
###############################################

function dtw(
        seq1::AbstractArray,
        seq2::AbstractArray,
        dist::SemiMetric=SqEuclidean()
    )
    dtw(Sequence(seq1),Sequence(seq2),dist)
end

function dtw(
        seq1::AbstractArray,
        seq2::AbstractArray,
        i2min::AbstractVector,
        i2max::AbstractVector,
        dist::SemiMetric = SqEuclidean()
    )
    dtw(Sequence(seq1),Sequence(seq2),i2min,i2max,dist)
end
