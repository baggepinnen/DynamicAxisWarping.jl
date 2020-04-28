# DynamicAxisWarping.jl

[![CI](https://github.com/baggepinnen/DynamicAxisWarping.jl/workflows/CI/badge.svg)](https://github.com/baggepinnen/DynamicAxisWarping.jl/actions)
[![codecov](https://codecov.io/gh/baggepinnen/DynamicAxisWarping.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/baggepinnen/DynamicAxisWarping.jl)


Dynamic Time Warping (DTW) and related algorithms in Julia.

- **Warning:** This package is under active development and the API is likely to break.
- This package supports arbitrary metrics and arbitrary "spaces", i.e., as long as you are passing a vector or higher dimensional array of something that your distance can operate on, you're good to go. Time is always considered to be the last dimension.

This package is not registered. Install using:
```julia
using Pkg
pkg"add https://github.com/baggepinnen/DynamicAxisWarping.jl"
```

## Simple usage

Inputs of dimension larger than 1 will be treated as sequences where time is in the last dimension. When using higher-dimensional series, make sure the provided distance accepts them.

Any distance implementing the [Distances.jl](https://github.com/JuliaStats/Distances.jl) interface works, as well as functions on the form `dist(x,y) -> ℝ`.

```julia
using DynamicAxisWarping, Distances, Plots
cost, i1, i2 = dtw(a,b, [dist=SqEuclidean()]; transportcost = 1)
cost, i1, i2 = fastdtw(a,b, dist, radius)
cost = dtw_cost(a, b, dist, radius) # Optimized method that only returns cost. Supports early stopping, see docstring. Can be made completely allocation free.

# dtw supports arbitrary upper and lower bound vectors constraining the warping path.
imin,imax = radiuslimits(5,20,20), plot([imin imax])
dtw(a, b, dist, imin, imax) # Cost eqivalent to dtw_cost(a, b, dist, 5)
```

## Plotting
```julia
dtwplot(a,b, [dist=SqEuclidean()]; transportcost = 1)
matchplot(a,b, [dist=SqEuclidean()])
```

## Find a short pattern in a long time series
The function `dtwnn` searches for a pattern in a long time series. It's currently not super well optimized, and *does not normalize* the data over each window (see roadmap issue).

```julia
using DynamicAxisWarping, Distances
radius = 5
a      = sin.(0.1 .* (1:100))     .+ 0.1 .* randn.()
b      = sin.(0.1 .* (1:100_000)) .+ 0.1 .* randn.()
res    = dtwnn(a, b, SqEuclidean(), radius, saveall=false, bsf_multiplier=1) # takes about 0.1s # DynamicAxisWarping.DTWSearchResult(0.4625287975222824, 73452, (prune_end = 79108, prune_env = 0))
plot([a b[eachindex(a) .+ (res.loc-1)]])
```

- `saveall` allows you to return the entire distance profile. This will take longer time to compute.
- `bsf_multiplier = 1`: If > 1, require lower bound to exceed `bsf_multiplier*best_so_far`. Will only have effect if `saveall = false`. This allows you to find several nearby points without having to compute the entire distance profile.

## Clustering and barycenter averaging
```julia
barycenter = dba(vector_of_arrays)
result     = dbaclust(data, nclust, ClassicDTW())
```
Note that `dba` is known to not always produce the best barycenters. See, e.g., ["Soft-DTW: a Differentiable Loss Function for Time-Series"](https://arxiv.org/pdf/1703.01541.pdf) or ["Spatio-Temporal Alignments: Optimal transport through space and time"](https://arxiv.org/pdf/1910.03860.pdf) for a method that produces better barycenters at the expense of a much higher computational cost.


### Optimizations
The following optimizations are implemented.
- [x] Endpoint lower bound pruning
- [ ] Envelope lower bound pruning
- [x] DTW early termination
- [ ] Online normalization
- [ ] Sorting of query series



### `transportcost`
`transportcost` adds an additional penalty multiplier for "transporting", i.e., deviations from the Euclidean matching. The standard DTW distance does not consider this added cost and the default is 1. A value greater than 1 multiplies the cost of moving horizontally or vertically in the coupling matrix, promoting a diagnoal move, corresponding to the standard Euclidean matching. The influence of the transport cost can be visualized with
```julia
a = sin.(1:100); b = sin.(1:100) .+ randn.();
dtwplot(a,b, transportcost=1)    # Default
dtwplot(a,b, transportcost=1.01) # Should be "more diagnoal"
dtwplot(a,b, transportcost=1.1)  # Should be almost completely diagnoal
```
You can try a `transportcost < 1` as well, but then it is preferable to make weird alignments and I'm not sure how much sense that would make.

## Combine with optimal transport
The distance between two datapoints can be any distance supporting the [Distances.jl](https://github.com/JuliaStats/Distances.jl/) interface.

See the file [`frequency_warping.jl`](https://github.com/baggepinnen/DynamicAxisWarping.jl/blob/master/examples/frequency_warping.jl) ([notebook](https://nbviewer.jupyter.org/github/baggepinnen/julia_examples/blob/master/frequency_warping.ipynb)) for an example combining dynamic time warping with optimal transport along the frequency axis for spectrograms. This example makes use of [SpectralDistances.jl](https://github.com/baggepinnen/SpectralDistances.jl).

#### Acknowledgements

This package is a fork of https://github.com/ahwillia/TimeWarp.jl which is no longer maintained.

Special thanks to Joseph Fowler ([@joefowler](https://github.com/joefowler)) who contributed a substantial portion of this code.

[build-img]: https://travis-ci.org/baggepinnen/DynamicAxisWarping.jl.svg?branch=master
[build-url]: https://travis-ci.org/baggepinnen/DynamicAxisWarping.jl
