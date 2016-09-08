# methods for estimating dtw #
abstract DTWMethod

type ClassicDTW <: DTWMethod end

immutable FastDTW <: DTWMethod
    radius::Int
end

# distance interface #
immutable DTWDistance{M<:DTWMethod, D<:SemiMetric} <: SemiMetric
    method::M
    dist::D
    DTWDistance(m::M,d::D) = new(m,d)
end

function DTWDistance{M<:DTWMethod,D<:SemiMetric}(
        d::D=SqEuclidean(),
        m::M=ClassicDTW()
    )
    DTWDistance{M,D}(m,d)
end

function DTWDistance{M<:DTWMethod,D<:SemiMetric}(
        m::M,
        d::D=SqEuclidean()
    )
    DTWDistance{M,D}(m,d)
end

Distances.evaluate(d::DTWDistance{ClassicDTW}, x, y) = dtw(x,y,d.dist)[1]
Distances.evaluate(d::DTWDistance{FastDTW}, x, y) = fastdtw(x,y,d.method.radius,d.dist)[1]

distpath(d::DTWDistance{ClassicDTW}, x, y) = dtw(x,y,d.dist)
distpath(d::DTWDistance{FastDTW}, x, y) = fastdtw(x,y,d.method.radius,d.dist)
