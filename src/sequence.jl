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

Base.size{N}(x::Sequence{N}) = (size(x.val,N),)
Base.eltype{N,T}(x::Sequence{N,T}) = T

@generated function Base.getindex{N}(x::Sequence{N}, i)
    :( x.val[@ntuple($N, (n-> n==$N ? i : Colon()))...] )
end

@generated function Base.setindex!{N}(x::Sequence{N}, val, i)
    :( x.val[@ntuple($N, (n-> n==$N ? i : Colon()))...] = val )
end

# convert a sequence back into an array
# in case sequences have different lenght, throw error
function seq_to_array{T<:Sequence}(seq::AbstractVector{T})
  len_seq = length(seq[1])
  arr = zeros(typeof(seq[1][1]),len_seq,length(seq))
  for i=1:length(seq)
    length(seq[i]) != len_seq ? error("Sequences do not have the same length, cannot construct array") : nothing
    arr[:,i] = seq[i][:]
  end
  return arr
end
