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
        symmetric in (true, false)

        kwargs = (N=N, M=M, verbose=false, metric=metric, symmetric = symmetric)

        data = prepare_gdtw(x, y; kwargs...)
        cost = DynamicAxisWarping.single_gdtw!(data)
        warp = copy(data.warp)

        # Check against graph based implementation
        # Should not change the result.
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
    for symmetric in (true, false)
        costs = [ gdtw(x,y; max_iters = i, symmetric=symmetric)[1] for i = 1:3]
        @test issorted(costs; rev=true)
    end
end

@testset "GDTW: test that warping helps" begin
    x = LinearInterpolation(rand(20))
    y = LinearInterpolation(rand(20))
    for symmetric in (true, false)
        cost, ϕ, ψ = gdtw(x, y; symmetric=symmetric)
        t = range(0, stop = 1, length=20)
        @test norm(x.(t) - y.(t)) >= norm( x.(ϕ.(t)) - y.(ψ.(t)) )
    end
end

@testset "GDTW: Test objective function" begin
    # short time span to make it easy
    ts = range(0, stop=π/12, length=128)
    x = LinearInterpolation(sin.(ts))
    y = LinearInterpolation(sin.(1.1 .* ts))
    metric = (x,y) -> norm(x-y)
    λinst = .1
    λcum = .1
    for symmetric in (false, true)
        # We need high `N` and `M` to have low discretization error
        cost, ϕ, ψ = gdtw(x, y; symmetric=symmetric, N = 400, M = 300, λinst=λinst, λcum=λcum)
        ϕ′ = x -> ForwardDiff.derivative(ϕ, x)
        res = quadgk(t -> metric(x(ϕ(t)),  y(ψ(t))) + λinst*(ϕ′(t)-1)^2 + λcum*(ϕ(t) - t)^2, 0, 1 )[1]
        @test res ≈ cost rtol=1e-2
    end
end

@testset "GDTW: test symmetry" begin
    x = LinearInterpolation(rand(20))
    y = LinearInterpolation(rand(20))
    cost_xy, ϕ_xy, ψ_xy = gdtw(x, y)
    cost_yx, ϕ_yx, ψ_yx = gdtw(y, x)
    @test cost_xy ≈ cost_yx
    t = range(0, stop = 1, length=100)
    @test ϕ_xy.(t) ≈ ψ_yx.(t)
    @test ψ_xy.(t) ≈ ϕ_yx.(t)
end

@testset "GDTW: `iterative_gdtw!`" begin
    x = LinearInterpolation(rand(20))
    y = LinearInterpolation(rand(20))
    t = range(0, stop = 1, length=100)
    for symmetric in (true, false), max_iters in (1, 2, 5)
        cost, ϕ, ψ = gdtw(x, y; symmetric = symmetric, max_iters = max_iters)

        data = prepare_gdtw(x, y; symmetric = symmetric)
        local cost2
        for i = 1:max_iters
            cost2 = iterative_gdtw!(data, max_iters = i)
        end
        ϕ2, ψ2 = gdtw_warpings(data)
        @test cost ≈ cost2
        @test ϕ.(t) ≈ ϕ2.(t)
        @test ψ.(t) ≈ ψ2.(t)
    end
end
