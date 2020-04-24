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
using DynamicAxisWarp, Distances
cost, i1, i2 = dtw(a,b, [dist=SqEuclidean()]; transportcost = 0)
cost, i1, i2 = fastdtw(a,b, [dist=SqEuclidean()])
dtwplot(a,b, [dist=SqEuclidean()]; transportcost = 0)
```
`transportcost` adds an additional penalty for "transporting", i.e., deviations from the Euclidean matching. The standard DTW distance does not consider this added cost and the default is 0.

See also function `dba` for barycenter averaging, but note that `dba` is known to not always produce the best barycenters. See, e.g., ["Spatio-Temporal Alignments: Optimal transport through space and time"](https://arxiv.org/pdf/1910.03860.pdf) for a method that produces better barycenters at the expense of a much higher computational cost.

#### Acknowledgements

Special thanks to Joseph Fowler ([@joefowler](https://github.com/joefowler)) who contributed a substantial portion of this code.

[build-img]: https://travis-ci.org/baggepinnen/DynamicAxisWarp.jl.svg?branch=master
[build-url]: https://travis-ci.org/baggepinnen/DynamicAxisWarp.jl
