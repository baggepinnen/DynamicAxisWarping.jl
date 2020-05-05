module DynamicAxisWarping

import ProgressMeter
using Statistics
using StatsBase
using Distances
using FillArrays
using ProgressMeter
using UnsafeArrays
using DataStructures
using SlidingDistancesBase
import SlidingDistancesBase: floattype, lastlength
using Plots, Plots.PlotMeasures # Plots required both for @layout and for the margins

export dtw,
       dtw_cost,
       DTWDistance,
       DTWMethod,
       DTW,
       FastDTW,
       distpath,
       dba,
       dbaclust,
       dtw_cost_matrix,
       DBAResult,
       fastdtw,
       radiuslimits,
       dtwnn,
       DTWWorkspace

export ZNormalizer,
       normalize

export stomp

include("utils.jl")

include("distance_interface.jl")
include("normalizers.jl")
include("dtw.jl")
include("dtwnn.jl")
include("dba.jl")
include("dbaclust.jl")
include("windowed_matrix.jl")
include("fastdtw.jl")
include("datasets/datasets.jl")
include("plots.jl")

end # module
