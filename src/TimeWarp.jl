module TimeWarp

import ProgressMeter
using StatsBase
using Parameters

# reexport to let user specify dtw(⋅,⋅,distance)
using Reexport
@reexport using Distances

export Sequence,
       SequenceArray,
       dtw,
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
include("sequence.jl")

include("distance_interface.jl")

include("dtw.jl")
include("dba.jl")
include("dbaclust.jl")
include("windowed_matrix.jl")
include("fastdtw.jl")
include("datasets.jl")
include("plots.jl")

end # module
