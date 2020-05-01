module DynamicAxisWarping

import ProgressMeter
using Statistics
using StatsBase
using Parameters
using Distances
using FillArrays
using ProgressMeter
using UnsafeArrays
using DataStructures
using LoopVectorization
using FFTW
using Plots, Plots.PlotMeasures

export dtw,
       dtw_cost,
       DTWDistance,
       DTWMethod,
       ClassicDTW,
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
include("mp.jl")

end # module
