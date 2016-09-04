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
type Sequence{N,T} <: AbstractArray{T,1}
    val::AbstractArray{T,N}
end

@generated Base.size{N}(x::Sequence{N}) = :((size(x.val,N),))

@generated function Base.getindex{N}(x::Sequence{N},i::Int)
    :( x.val[@ntuple($N, (n-> n==N ? i : Colon()))...] )
end

"""
    SequenceArray(::AbstractArray)

A `SequenceArray` holds a collection of time series that all have the
same number of time samples. It is a thin wrapper around an array 
that allows you to index it like a matrix, with the first dimension
corresponding to time point and the second dimension corresponding to
the time series.

The simplest example is a collection of 1-dimensional time series. The
following example shows a `SequenceArray` holding 20 time series at
100 time points:

```
data = randn(100,20)
seqs = SequenceArray(data)
data[t,s] == seqs[t,s] # gives sequence s at time t
```

In this case a `SequenceArray` behaves like a typical matrix object. The
real utility of the `SequenceArray` comes when you have multi-dimensional
time series. The next example shows a collection of 10-dimensional time
series:

```
data = randn(10,100,20)
seqs = SequenceArray(data)
data[:,t,s] == seqs[t,s] # gives sequence s at time t
```

This generalizes to higher-order arrays (e.g. collections of time-series
where each observation is a matrix).
"""
type SequenceArray{N,T} <: AbstractArray{T,2}
    val::AbstractArray{T,N}

    function SequenceArray{N,T}(X::AbstractArray{T,N})
        N<2 && throw(ArgumentError("SequenceArray requires a matrix or higher-order array as input."))
        new(X)
    end
end
SequenceArray{N,T}(X::AbstractArray{T,N}) = SequenceArray{N,T}(X)

@generated Base.size{N}(x::SequenceArray{N}) = :((size(x.val,N-1),size(x.val,N)))

@inline Base.getindex(x::SequenceArray,i) = getindex(x,ind2sub(x,i)...)

@generated function Base.getindex{N}(x::SequenceArray{N},t,s)
  quote
    idx = @ntuple($N,
        n-> begin
            if n==$(N-1)
                t
            elseif n==$N
                s
            else
                Colon()
            end
        end
    )
    x.val[idx...]
  end
end

"""
A `SequenceList` holds a collection of time series that may have a 
different number of time samples. Underneath, it stores a vector
of array views.
"""
type SequenceList end

"""
SequenceCollection

A sequences collection holds a dataset containing multiple time-series.
"""
typealias SequenceCollection Union{SequenceArray,SequenceList}
