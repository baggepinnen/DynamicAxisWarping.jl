module TimeWarp

using Reexport
@reexport using Distances
import Compat:view

export Sequence,
       dtw,
       dba,
       dtw_cost_matrix,
       DBAResult,
       fastdtw

include("sequence.jl")
include("dtw.jl")
include("dba.jl")
include("windowed_matrix.jl")
include("fastdtw.jl")
include("datasets.jl")
include("plots.jl")

end # module
