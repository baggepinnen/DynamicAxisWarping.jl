import Base: @propagate_inbounds, @boundscheck, getindex

abstract type AbstractNormalizer{T,N} <: AbstractArray{T,N} end

advance!(x) = 0 # Catches everything that's not a normalizer

setup_normalizer(::Nothing, q, y) = q, y

mutable struct ZNormalizer{T} <: AbstractNormalizer{T,1}
    x::Vector{T}
    n::Int
    μ::T
    σ::T
    s::T
    ss::T
    i::Int
end

function normalize(::ZNormalizer, q)
    q = q .- mean(q)
    q ./= std(q, corrected=false, mean=0)
end

setup_normalizer(z::Type{ZNormalizer}, q, y) = normalize(z, q), z(y, length(q))

function ZNormalizer(x::AbstractArray{T}, n) where T
    @assert length(x) >= n
    s = ss = zero(T)
    @inbounds @simd for i in 1:n
        s += x[i]
        ss += x[i]^2
    end
    μ = s/n
    σ = sqrt(ss/n - μ^2)
    ZNormalizer(x, n, μ, σ, s, ss, 1)
end

@propagate_inbounds function advance!(z::ZNormalizer{T}) where T

    @boundscheck if z.i + z.n > length(z.x)
        return z.i += 1
    end

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

@inline @propagate_inbounds function getindex(z::ZNormalizer, ::typeof(!), i)
    (z[i]-z.μ) / z.σ
end

Statistics.mean(z::ZNormalizer) = z.μ
Statistics.std(z::ZNormalizer) = z.σ

Base.length(z::ZNormalizer) = z.n
Base.size(z::ZNormalizer) = (z.n,)
