#
# Define a "windowed matrix" type, where a maximum and minimum row is specified for
# each column.  Trying to read the matrix outside that range yields a default value
# (this defaults to Inf, though it can be specified when constructing the object).
# Trying to write outside that range throws a BoundsError().
#
# This is provided for efficient implementation of restricted Dynamic Time Warping
# algorithms.
#
# Joe Fowler
# NIST Boulder Laboratories
# December 2014
#

struct WindowedMatrix{T<:Real} <: AbstractArray{T,2}
    nrow::Int
    ncol::Int
    ncells::Int

    cost::Vector{T}
    rowmin::Vector{Int}
    rowmax::Vector{Int}
    rowspercol::Vector{Int}
    idxcol::Vector{Int}  # Index (in cost) of 1st element in each column
    defaultval::T

    function WindowedMatrix{T}(
        rmin::Vector{Int},
        rmax::Vector{Int},
        default::T,
    ) where {T<:Real}
        rowmin = copy(rmin)
        rowmax = copy(rmax)
        rowspercol = rowmax .- rowmin .+ 1
        ncells = sum(rowspercol)
        nrow = maximum(rowmax)
        ncol = length(rowmax)
        cost = zeros(T, ncells)
        idxcol = 1 .+ vcat([0], cumsum(rowspercol[1:end-1])) # index of 1st element per column

        new(nrow, ncol, ncells, cost, rowmin, rowmax, rowspercol, idxcol, default)
    end
end


# If the default matrix value is given, then the matrix takes on its type
function WindowedMatrix(rmin::Vector{Int}, rmax::Vector{Int}, default::T) where {T<:Real}
    WindowedMatrix{T}(rmin, rmax, default)
end

# If no default matrix value is given, it will be Inf, and matrix will hold
# Float64 values.
WindowedMatrix(rmin::Vector{Int}, rmax::Vector{Int}) =
    WindowedMatrix{Float64}(rmin, rmax, Inf)

Base.size(W::WindowedMatrix) = W.nrow, W.ncol
Base.size(W::WindowedMatrix, i) = size(W)[i]

@inline function Base.getindex(W::WindowedMatrix, r::Integer, c::Integer)
    Base.@boundscheck if c < 1 || c > W.ncol || r < W.rowmin[c] || r > W.rowmax[c]
        return W.defaultval
    end

    offset = r - W.rowmin[c]
    W.cost[W.idxcol[c]+offset]
end


@inline function Base.setindex!(W::WindowedMatrix, val, r::Integer, c::Integer)
    Base.@boundscheck if c < 1 || c > W.ncol || r < W.rowmin[c] || r > W.rowmax[c]
        throw(BoundsError())
    end

    offset = r - W.rowmin[c]
    W.cost[W.idxcol[c]+offset] = val
    return W
end
