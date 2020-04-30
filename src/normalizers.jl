import Base: @propagate_inbounds, @boundscheck, getindex
import LinearAlgebra: normalize

abstract type AbstractNormalizer{T,N} <: AbstractArray{T,N} end

advance!(x) = 0 # Catches everything that's not a normalizer

setup_normalizer(::Val{Nothing}, q, y) = q, y

normalize(::Val{Nothing}, q) = q


normalize(T::Type, q) = normalize(Val(T), q)

"""
    ZNormalizer{T} <: AbstractNormalizer{T, 1}

Utility object that normalizes an input array on the fly.
Works like a vector, index into window with `z[i]`, index into normalized window with `z[!,i]`.

# Arguments:
- `x::Vector{T}`: The vecetor to operate on
- `n::Int`: The lenght of a window
- `μ::T`: mean over windo
- `σ::T`: std over window
- `s::T`: sum over windoe
- `ss::T`: sum^2 over window
- `i::Int`: the current index
- `buffer::Vector{T}`: stores calculated normalized values, do not access this, it might not be complete
- `bufi::Int`: index of last normalized value
"""
mutable struct ZNormalizer{T} <: AbstractNormalizer{T,1}
    x::Vector{T}
    n::Int
    μ::T
    σ::T
    s::T
    ss::T
    i::Int
    buffer::Vector{T}
    bufi::Int
end

function ZNormalizer(x::AbstractArray{T}, n) where T
    @assert length(x) >= n
    s = ss = zero(T)
    @inbounds @simd for i in 1:n
        s += x[i]
        ss += x[i]^2
    end
    μ = s/n
    σ = sqrt(ss/n - μ^2)
    buffer = similar(x, n)
    ZNormalizer(x, n, μ, σ, s, ss, 0, buffer, 0)
end

function normalize(::Val{ZNormalizer}, q)
    q = q .- mean(q)
    q ./= std(q, corrected=false, mean=0)
end

@inline @propagate_inbounds function normalize(::Val{ZNormalizer}, z::ZNormalizer)
    if z.bufi == z.n
        return z.buffer
    end
    for i = z.bufi+1:z.n
        z[!, i, true] # This populates the buffer
    end
    return z.buffer
end

setup_normalizer(z::Val{ZNormalizer}, q, y) = normalize(z, q), ZNormalizer(y, length(q))

@propagate_inbounds function advance!(z::ZNormalizer{T}) where T

    if z.i == 0
        return z.i = 1
    end
    @boundscheck if z.i + z.n > length(z.x)
        return z.i += 1
    end
    z.bufi = 0

    # Remove old point
    x = z.x[z.i]
    z.s -= x
    z.ss -= x^2

    # Add new point
    x = z.x[z.i+z.n]
    z.s += x
    z.ss += x^2
    z.μ = z.s/z.n
    z.σ = sqrt(z.ss/z.n - z.μ^2)
    z.i += 1
end


@inline @propagate_inbounds function getindex(z::ZNormalizer, i)
    @boundscheck 1 <= i <= z.n || throw(BoundsError(z,i))
    xi = i+z.i-1
    @boundscheck xi <= length(z.x) || throw(BoundsError(z,i))
    z.x[xi]
end

@inline @propagate_inbounds function getindex(z::ZNormalizer, ::typeof(!), i, inorder = i == z.bufi + 1)
    y = (z[i]-z.μ) / z.σ
    if inorder
        z.bufi = i
        z.buffer[i] = y
    end
    y
end

@inline @propagate_inbounds function getindex(z::ZNormalizer, ::typeof(!), i::AbstractRange)
    @boundscheck (i[1] == z.i && length(i) == z.n) || throw(ArgumentError("ZNormalizers can only be indexed by ranges corresponding to their current state. Got range $i but state was $(z.i)"))
    z
end

Statistics.mean(z::ZNormalizer) = z.μ
Statistics.std(z::ZNormalizer) = z.σ

lastlength(z::ZNormalizer) = length(z.x)
Base.length(z::ZNormalizer) = z.n
Base.size(z::ZNormalizer) = (z.n,)
