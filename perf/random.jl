using BenchmarkTools
using TimeWarp

x = randn(5000)
y = randn(5000)

bd = @benchmark dtw($x,$y)
bf1 = @benchmark fastdtw($x,$y,1)
bf10 = @benchmark fastdtw($x,$y,10)

### HISTORY ###

## a17d716 ##

# julia> bd
# BenchmarkTools.Trial: 
#   samples:          1
#   minimum time:     12.74 s (7.60% GC)

# julia> bf1
# BenchmarkTools.Trial: 
#   samples:          123
#   minimum time:     35.41 ms (10.94% GC)

# julia> bf10
# BenchmarkTools.Trial: 
#   samples:          20
#   minimum time:     221.13 ms (6.86% GC)
