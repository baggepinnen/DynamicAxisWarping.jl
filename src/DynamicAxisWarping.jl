module DynamicAxisWarping

import ProgressMeter
using StatsBase
using Parameters
using Distances
using FillArrays
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

include("utils.jl")

include("distance_interface.jl")

include("dtw.jl")
include("dtwnn.jl")
include("dba.jl")
include("dbaclust.jl")
include("windowed_matrix.jl")
include("fastdtw.jl")
include("datasets/datasets.jl")
include("plots.jl")

end # module
