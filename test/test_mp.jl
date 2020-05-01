
x = randn(30)
for w = 2:10
    m,s = DynamicAxisWarping.running_mean_std(x, w)
    @test length(m) == length(s) == 30-w+1
    @test m[1] ≈ mean(x[1:w])
    @test s[1] ≈ std(x[1:w], corrected=false)

    @test m[2] ≈ mean(x[2:w+1])
    @test s[2] ≈ std(x[2:w+1], corrected=false)
end


t = range(0, stop=1, step=1/10)
y0 = sin.(2pi .* t)

T = [randn(50); y0; randn(50); y0; randn(50)]

P,I = stomp(T, length(y0))

# plot(T, layout=2)
# plot!(P, sp=2)

m = findmin(P)
@test m[1] < 1e-6
@test m[2] == 51 || m[2] == 112

# @btime stomp($randn(Float32, 10_000), 20)
# @profiler stomp(randn(10_000), 20)


# Q = randn(5)
# T = randn(10)
# d1 = DynamicAxisWarping.window_dot(Q,T)
# d2 = DynamicAxisWarping.window_dot2(Q,T)
# d3 = DynamicAxisWarping.window_dot3(Q,T)
# [d1 d2 d3]





# @time stomp(randn(Float32, 2^17), 256)
