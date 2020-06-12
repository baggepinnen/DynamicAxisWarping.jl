include("gdtw_graph_implementation.jl")

@testset "GDTW: test against graph implementation" begin
    ts = range(0, stop=4π, length=50)
    rnd = () -> rand(2)
    x = LinearInterpolation(sin.(ts)' .* rnd())
    y = LinearInterpolation(sin.(1.1 .* ts)' .* rnd())

    N = 30
    M = 20

    # choices that change the result
    for metric in (Distances.SqEuclidean(), (x, y) -> norm(x - y, 1)),
        max_iters in (1, 2, 3)

        kwargs = (N=N, M=M, verbose=false, max_iters=max_iters, metric=metric)

        data = DynamicAxisWarping.prepare_gdtw(x, y; kwargs...)
        cost = DynamicAxisWarping.single_gdtw!(data)
        warp = copy(data.warp)

        # Check against graph based implementation
        for algo in (LightGraphs.dijkstra_shortest_paths,
 LightGraphs.desopo_pape_shortest_paths, LightGraphs.bellman_ford_shortest_paths)
            new_cost = DynamicAxisWarping.single_gdtw!(data, algo)
            @test new_cost ≈ cost
            @test data.warp ≈ warp
        end
    end

end


@testset "GDTW: test that refinements reduce cost" begin
    x = LinearInterpolation(rand(2, 20))
    y = LinearInterpolation(rand(2, 20))
    costs = [ gdtw(x,y; max_iters = i)[1] for i = 1:3]
    @test issorted(costs; rev=true)
end

@testset "GDTW: test that warping helps" begin
    x = LinearInterpolation(rand(20))
    y = LinearInterpolation(rand(20))
    cost, ϕ = gdtw(x, y)
    t = range(0, stop = 1, length=20)
    @test norm(x.(t) - y.(t)) >= norm( x.(ϕ.(t)) - y.(t) )
end

@testset "GDTW: Test objective function" begin
    # short time span to make it easy
    ts = range(0, stop=π/12, length=128)
    x = LinearInterpolation(sin.(ts))
    y = LinearInterpolation(sin.(1.1 .* ts))

    λinst = .1
    λcum = .1
    # We need high `N` and `M` to have low discretization error
    cost, ϕ = gdtw(x, y; N = 400, M = 300, λinst=λinst, λcum=λcum)
    ϕ′ = x -> ForwardDiff.derivative(ϕ, x)
    res = quadgk(t -> metric(x(ϕ(t)),  y(t)) + λinst*(ϕ′(t)-1)^2 + λcum*(ϕ(t) - t)^2, 0, 1 )[1]
    @test res ≈ cost rtol=1e-2
end
