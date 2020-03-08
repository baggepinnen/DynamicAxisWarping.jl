using Base.Cartesian

"""
    Sequence{N,T}(::AbstractArray)

A `Sequence` is a thin wrapper around a multi-dimensional array that
allows you to index it like a vector. For example, the following
code shows a 10-dimensional time series, sampled at 100 time points:

```
data = randn(10,100)
seq = Sequence(data)
seq[1] == data[:,1] # true
```

This generalizes to higher-order arrays. Consider a time series
where a 10x10 matrix of data is collected at 100 time points:

```
data = randn(10,10,100)
seq = Sequence(data)
seq[1] == data[:,:,1] # true
```
"""
struct Sequence{N,T} <: AbstractArray{T,1}
    val::AbstractArray{T,N}
end

Base.size(x::Sequence{N}) where N = (size(x.val,N),)
Base.eltype(x::Sequence{N,T}) where {N, T} = T

@generated function Base.getindex(x::Sequence{N}, i) where N
    :( x.val[@ntuple($N, (n-> n==$N ? i : Colon()))...] )
end

@generated function Base.setindex!(x::Sequence{N}, val, i) where N
    :( x.val[@ntuple($N, (n-> n==$N ? i : Colon()))...] = val )
end

function seq_to_array(seq::Sequence{N}) where N
    return cat(seq..., dims=N)
end
