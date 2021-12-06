using BenchmarkTools
# testing multi-threading speedup
# increase number of vectors in data
data_big = [randn(300) .+ 2(i ÷ 5) for i = 0:19]
nclust = 10
result = [@benchmark dbaclust(
    $data_big,
    $nclust,
    DTW(10);
    n_init = 20,
    iterations = 10,
    n_jobs = $j,
    show_progress=false
) for j in [1, -1]]
println("Single core: ", result[1])
println("Multi core: ", result[2])
@test mean(result[1].times) > mean(result[2].times)