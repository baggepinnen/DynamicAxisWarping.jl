# testing multi-threading speedup
# increase number of vectors in data
data_big = [randn(300) .+ 2(i ÷ 5) for i = 0:999]
nclust = 20
println("Single core: ")
@time dbaclust(data_big, nclust, DTW(10); n_init = 1, iterations = 10, show_progress=false, threaded=false);
println("Multi core: ")
@time dbaclust(data_big, nclust, DTW(10); n_init = 1, iterations = 10, show_progress=false, threaded=true);
