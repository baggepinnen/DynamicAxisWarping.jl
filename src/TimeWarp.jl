module TimeWarp

using Distances
import Compat:view

export dtw,
       dba,
       DBAResult,
       fastdtw

include("dtw.jl")
include("dba.jl")
include("windowed_matrix.jl")
include("fastdtw.jl")
include("datasets.jl")

end # module
