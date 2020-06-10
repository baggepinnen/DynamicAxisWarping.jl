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
import SlidingDistancesBase: floattype, lastlength, setup_normalizer
using LoopVectorization
using Plots, Plots.PlotMeasures # Plots required both for @layout and for the margins
using Requires
using UnPack

export dtw,
       dtw_cost,
       soft_dtw_cost,
       DTWDistance,
       DTW,
       SoftDTW,
       FastDTW,
       distpath,
       dba,
       dbaclust,
       dtw_cost_matrix,
       soft_dtw_cost_matrix,
       DBAResult,
       fastdtw,
       radiuslimits,
       dtwnn,
       DTWWorkspace,
       DTWSearchResult,
       GDTWCache,
       sparse_distmat,
       gdtw,
       LinearInterpolation

export ZNormalizer,
       IsoZNormalizer,
       normalize

export stomp

include("utils.jl")

include("distance_interface.jl")
include("dtw.jl")
include("gdtw.jl")
include("dtwnn.jl")
include("dba.jl")
include("dbaclust.jl")
include("windowed_matrix.jl")
include("fastdtw.jl")
include("datasets/datasets.jl")
include("plots.jl")

function __init__()
    @require MatrixProfile = "24e37439-14ec-4097-bda3-6a65822e2305" include("matrix_profile.jl")
end


end # module
