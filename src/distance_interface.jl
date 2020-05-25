# methods for estimating dtw #

abstract type DTWDistance{D <: Union{Function, Distances.PreMetric}} end


"""
    struct DTW{D} <: DTWDistance{D}

# Arguments:
- `radius`: The maximum allowed deviation of the matching path from the diagonal
- `dist`: Inner distance
- `transportcost` If >1, an additional penalty factor for non-diagonal moves is added.
"""
Base.@kwdef struct DTW{D} <: DTWDistance{D}
    "The maximum allowed deviation of the matching path from the diagonal"
    radius::Int
    dist::D = SqEuclidean()
    "If >1, an additional penalty factor for non-diagonal moves is added."
    transportcost::Float64 = 1.0
    DTW(r,dist=SqEuclidean(),transportcost=1) = new{typeof(dist)}(r,dist,transportcost)
end

"""
    struct SoftDTW{D, T} <: DTWDistance{D}

# Arguments:
- `γ`: smoothing parameter
- `dist`
- `transportcost`
"""
Base.@kwdef struct SoftDTW{D,T} <: DTWDistance{D}
    γ::T
    "The maximum allowed deviation of the matching path from the diagonal"
    dist::D = SqEuclidean()
    "If >1, an additional penalty factor for non-diagonal moves is added."
    transportcost::Float64 = 1.0
    SoftDTW(γ, dist=SqEuclidean(),transportcost=1) = new{typeof(dist), typeof(γ)}(γ,dist,transportcost)
end


"""
    struct FastDTW{D} <: DTWDistance{D}

- `radius`
- `dist` inner distance
"""
Base.@kwdef struct FastDTW{D} <: DTWDistance{D}
    radius::Int
    dist::D = SqEuclidean()
    FastDTW(r, dist=SqEuclidean()) = new{typeof(dist)}(r, dist)
end



Distances.evaluate(d::DTW, x, y) = dtw_cost(x, y, d.dist, d.radius)
Distances.evaluate(d::SoftDTW, x, y) = soft_dtw_cost(x, y, d.dist, γ=d.γ)
Distances.evaluate(d::FastDTW, x, y) =
    fastdtw(x, y, d.dist, d.radius)[1]

distpath(d::DTW, x, y) = dtw(x, y, d.dist)
distpath(d::DTW, x, y, i2min::AbstractVector, i2max::AbstractVector) =
    dtw(x, y, i2min, i2max, d.dist)
distpath(d::FastDTW, x, y) = fastdtw(x, y, d.dist, d.radius)

(d::DTWDistance)(x,y) = Distances.evaluate(d,x,y)

"""
    distance_profile(d::DTWDistance, Q::AbstractArray{S}, T::AbstractArray{S}; kwargs...) where S

Optimized method for computing the distance profile using DTW distances. kwargs are sent to [`dtwnn`](@ref).
"""
function SlidingDistancesBase.distance_profile(d::DTWDistance, Q::AbstractArray{S}, T::AbstractArray{S}; kwargs...) where S
    m = lastlength(Q)
    n = lastlength(T)
    n >= m || throw(ArgumentError("Q cannot be longer than T"))
    l = n-m+1
    res = dtwnn(Q, Y, d.dist, d.radius; saveall=true, kwargs...)
    res.dists
end
