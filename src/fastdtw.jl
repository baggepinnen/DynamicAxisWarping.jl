
"""
    cost,i1,i2 = fastdtw(seq1,seq2,radius,[dist=SqEuclidean])

Computes FastDTW approximation to the DTW, described in Salvador & Chan,
Intelligent Data Analysis (2007).
"""
function fastdtw(
        seq1::AbstractVector,
        seq2::AbstractVector,
        radius::Int,
        dist::SemiMetric=SqEuclidean()
    )

    MinSize = max(radius + 2, 10)
    N1 = length(seq1)
    N2 = length(seq2)
    if N1 <= MinSize || N2 <= MinSize
        return (dtw(seq1, seq2, dist))
    end

    # Call recursively on a pair of sequences half this length
    compressed1 = compress2(seq1)
    compressed2 = compress2(seq2)
    _cost, lowrescol, lowresrow = fastdtw(compressed1, compressed2, radius, dist)

    # Now resample that path to the finer resolution, find the correct
    # window around it, and get the DTW given that window.
    hirescol, hiresrow = expandpath(lowrescol, lowresrow, N1, N2)
    idx2min, idx2max = computewindow(hirescol, hiresrow, radius)
    cost1, newcol, newrow = dtw(seq1, seq2, idx2min, idx2max, dist)
end


# Given a path through low-res space, generate an approximate path
# through high-res space. It should have dimension Ncol x Nrow

@inbounds function expandpath(lowrescol, lowresrow, Ncol, Nrow)
    @assert div(Ncol+1,2) == lowrescol[end]
    @assert div(Nrow+1,2) == lowresrow[end]
    Np = length(lowrescol)
    @assert Np == length(lowresrow)

    hirescol = zeros(eltype(lowrescol), 2*Np)
    hiresrow = zeros(eltype(lowresrow), 2*Np)
    hirescol[1] = hiresrow[1] = c = r = 1
    for i=1:Np-1
        # Select plan according to the next move in lowres path.
        if lowrescol[i+1] == lowrescol[i]  # Next move is up
            r += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            r += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r

        elseif lowresrow[i+1] == lowresrow[i] # Next move is sideways
            c += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            c += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r

        else  # Next move is diagonal.
            c += 1; r += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            c += 1; r += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r
        end
    end
    hirescol[end] = Ncol
    hiresrow[end] = Nrow
    # When expanding to an odd numbered size, it's possible to repeat
    # the last step.  Fix that:
    if hirescol[end]==hirescol[end-1] && hiresrow[end]==hiresrow[end-1]
        hirescol = hirescol[1:end-1]
        hiresrow = hiresrow[1:end-1]
    end
    hirescol, hiresrow
end

# yshort = compress2(y)
#   Returns a shortened time series that is half the length of the input sequence.
#   The length of the compressed sequence is always even.
function compress2(seq::AbstractVector)
    # Navg = div(length(seq), 2)
    evenseq = 0.5*(seq[1:2:end-1]+seq[2:2:end])
    if length(seq)%2 == 1
        return vcat(evenseq, [seq[end]])
    end
    evenseq
end



# Given the lists of (col,row) indices for the optimal path, compute a "window"
# around that path of the given radius.
# Returns (rowmin, rowmax), each a vector of length pathcols[end], representing
# for each column, the minimum and maximum row numbers used in that column.

@inbounds function computewindow(pathcols, pathrows, radius)
    Np = length(pathcols)
    @assert Np == length(pathrows)
    Ncol = pathcols[end]
    Nrow = pathrows[end]

    # Find the min/max row at each column in the path.
    pathmin = zeros(Int, Ncol)
    pathmax = zeros(Int, Ncol)
    for i=1:Np
        c,r = pathcols[i], pathrows[i]
        pathmax[c] = r
        if pathmin[c] == 0
            pathmin[c] = r
        end
    end

    # The window in each column for "radius" r starts at the pathmin
    # of the rth-previous column and ends at the pathmax of the
    # rth-next column, plus (in each case) the radius.
    if radius < Ncol-1 && radius < Nrow-1
        rowmin = vcat(fill(1,radius), pathmin[1:end-radius] .- radius)
        rowmax = vcat(pathmax[radius+1:end] .+ radius, fill(Nrow,radius))

        # Window values must be in the range [1:Nrow].
        for c=1:Ncol
            if rowmin[c]<1; rowmin[c]=1; end
            if rowmax[c]>Nrow; rowmax[c]=Nrow; end
        end
    else
        rowmin = fill(1,Ncol)
        rowmax = fill(Nrow,Ncol)
    end
    rowmin, rowmax
end
