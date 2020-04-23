module TimeWarp

import ProgressMeter
using StatsBase
using Parameters
using Distances
using Plots, Plots.PlotMeasures

export dtw,
       DTWDistance,
       DTWMethod,
       ClassicDTW,
       FastDTW,
       distpath,
       dba,
       dbaclust,
       dtw_cost_matrix,
       DBAResult,
       fastdtw

include("utils.jl")

include("distance_interface.jl")

include("dtw.jl")
include("dba.jl")
include("dbaclust.jl")
include("windowed_matrix.jl")
include("fastdtw.jl")
include("datasets/datasets.jl")
include("plots.jl")

end # module
