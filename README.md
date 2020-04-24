# DynamicAxisWarp.jl

[![CI](https://github.com/baggepinnen/DynamicAxisWarp.jl/workflows/CI/badge.svg)](https://github.com/baggepinnen/DynamicAxisWarp.jl/actions)
[![codecov](https://codecov.io/gh/baggepinnen/DynamicAxisWarp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/baggepinnen/DynamicAxisWarp.jl)


Dynamic Time Warping (DTW) and related algorithms in Julia.

This package is a fork of https://github.com/ahwillia/DynamicAxisWarp.jl which is no longer maintained.

This package isn't officially registered. Install using:

```julia
using Pkg
pkg"add https://github.com/baggepinnen/DynamicAxisWarp.jl"
```

## Usage

```julia
using DynamicAxisWarp, Distances, Plots
cost, i1, i2 = dtw(a,b, [dist=SqEuclidean()]; transportcost = 1)
cost, i1, i2 = fastdtw(a,b, [dist=SqEuclidean()])
dtwplot(a,b, [dist=SqEuclidean()]; transportcost = 1)
centers, clustids, result = dbaclust(data, nclust, FastDTW(10))
imin,imax = radiuslimits(5,20,20), plot([imin imax])
```
`transportcost` adds an additional penalty multiplier for "transporting", i.e., deviations from the Euclidean matching. The standard DTW distance does not consider this added cost and the default is 1. A value greater than 1 multiplies the cost of moving horizontally or vertically in the coupling matrix, promoting a diagnoal move, corresponding to the standard Euclidean matching. The influence of the transport cost can be visualized with
```julia
a = sin.(1:100); b = sin.(1:100) .+ randn.();
dtwplot(a,b, transportcost=1)    # Default
dtwplot(a,b, transportcost=1.01) # Should be "more diagnoal"
dtwplot(a,b, transportcost=1.1)  # Should be almost completely diagnoal
```
You can try a `transportcost < 1` as well, but then it is preferable to make weird alignments and I'm not sure how much sense that would make.

See also function `dba` for barycenter averaging, but note that `dba` is known to not always produce the best barycenters. See, e.g., ["Soft-DTW: a Differentiable Loss Function for Time-Series"](https://arxiv.org/pdf/1703.01541.pdf) or ["Spatio-Temporal Alignments: Optimal transport through space and time"](https://arxiv.org/pdf/1910.03860.pdf) for a method that produces better barycenters at the expense of a much higher computational cost.

#### Combine with optimal transport
See the file [`frequency_warping.jl`](https://github.com/baggepinnen/DynamicAxisWarp.jl/blob/master/examples/frequency_warping.jl) for an example combining dynamic time warping with optimal transport along the frequency axis for spectrograms. This example makes use of [SpectralDistances.jl](https://github.com/baggepinnen/SpectralDistances.jl).

#### Acknowledgements

Special thanks to Joseph Fowler ([@joefowler](https://github.com/joefowler)) who contributed a substantial portion of this code.

[build-img]: https://travis-ci.org/baggepinnen/DynamicAxisWarp.jl.svg?branch=master
[build-url]: https://travis-ci.org/baggepinnen/DynamicAxisWarp.jl
