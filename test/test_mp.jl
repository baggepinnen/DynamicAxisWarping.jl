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
