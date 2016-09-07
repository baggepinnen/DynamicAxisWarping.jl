using BenchmarkTools
using TimeWarp

srand(1234)
a = randn(400)
b = randn(400)
x = randn(5000)
y = randn(5000)

dtw_small_bench = @benchmark dtw($a,$b)
dtw_big_bench = @benchmark dtw($x,$y)
fastdtw_r1_bench = @benchmark fastdtw($x,$y,1)
fastdtw_r10_bench = @benchmark fastdtw($x,$y,10)

### HISTORY ###

##
# 1c46fca
# - added sizehint! to N⋅logN
# - very small speedup on large dtw
##

# julia> dtw_small_bench
# BenchmarkTools.Trial: 
#   samples:          61
#   minimum time:     66.47 ms (4.42% GC)

# julia> dtw_big_bench
# BenchmarkTools.Trial: 
#   samples:          1
#   minimum time:     12.85 s (8.09% GC)

# julia> fastdtw_r1_bench
# BenchmarkTools.Trial: 
#   samples:          129
#   minimum time:     35.19 ms (10.92% GC)

# julia> fastdtw_r10_bench
# BenchmarkTools.Trial: 
#   samples:          22
#   minimum time:     218.10 ms (6.72% GC)


##
# a17d716
##

# julia> dtw_small_bench
# BenchmarkTools.Trial: 
#   samples:          74
#   minimum time:     64.62 ms (4.53% GC)

# julia> dtw_big_bench
# BenchmarkTools.Trial: 
#   samples:          1
#   minimum time:     14.30 s (7.35% GC)

# julia> fastdtw_r1_bench
# BenchmarkTools.Trial: 
#   samples:          118
#   minimum time:     34.78 ms (10.94% GC)

# julia> fastdtw_r10_bench
# BenchmarkTools.Trial: 
#   samples:          22
#   minimum time:     217.50 ms (7.09% GC)
