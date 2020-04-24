# DynamicAxisWarp.jl

[![Build Status][build-img]][build-url]
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)


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
```
`transportcost` adds an additional penalty multiplier for "transporting", i.e., deviations from the Euclidean matching. The standard DTW distance does not consider this added cost and the default is 1. A value greater than 1 multiplies the cost of moving horizontally or vertically in the coupling matrix, promoting a diagnoal move, corresponding to the standard Euclidean matching. The influence of the transport cost can be visualized with
```julia
a = sin.(1:100); b = sin.(1:100) .+ randn.();
dtwplot(a,b, transportcost=1) # Default
dtwplot(a,b, transportcost=1.01) # Should be "more diagnoal"
dtwplot(a,b, transportcost=1.1) # Should be almost completely diagnoal
```
You can try a `transportcost < 1` as well, but then it is preferable to make weird alignments and I'm not sure how much sense that would make.

See also function `dba` for barycenter averaging, but note that `dba` is known to not always produce the best barycenters. See, e.g., ["Spatio-Temporal Alignments: Optimal transport through space and time"](https://arxiv.org/pdf/1910.03860.pdf) for a method that produces better barycenters at the expense of a much higher computational cost.

#### Acknowledgements

Special thanks to Joseph Fowler ([@joefowler](https://github.com/joefowler)) who contributed a substantial portion of this code.

[build-img]: https://travis-ci.org/baggepinnen/DynamicAxisWarp.jl.svg?branch=master
[build-url]: https://travis-ci.org/baggepinnen/DynamicAxisWarp.jl
