# methods for estimating dtw #
abstract type DTWMethod end

struct DTW <: DTWMethod
    "The maximum allowed deviation of the matching path from the diagonal"
    radius::Int
    "If >1, an additional penalty factor for non-diagonal moves is added."
    transportcost::Float64
    DTW(r,transportcost=1) = new(r,transportcost)
end


struct FastDTW <: DTWMethod
    radius::Int
end

# distance interface #
Base.@kwdef struct DTWDistance{M<:DTWMethod,D<:SemiMetric} <: SemiMetric
    method::M
    dist::D = SqEuclidean()
end

DTWDistance(m::DTWMethod) = DTWDistance(m, SqEuclidean())

Distances.evaluate(d::DTWDistance{DTW}, x, y) = dtw_cost(x, y, d.dist, d.method.radius)
Distances.evaluate(d::DTWDistance{FastDTW}, x, y) =
    fastdtw(x, y, d.dist, d.method.radius)[1]

distpath(d::DTWDistance{DTW}, x, y) = dtw(x, y, d.dist)
distpath(d::DTWDistance{DTW}, x, y, i2min::AbstractVector, i2max::AbstractVector) =
    dtw(x, y, i2min, i2max, d.dist)
distpath(d::DTWDistance{FastDTW}, x, y) = fastdtw(x, y, d.dist, d.method.radius)

(d::DTWDistance)(x,y) = Distances.evaluate(d,x,y)
