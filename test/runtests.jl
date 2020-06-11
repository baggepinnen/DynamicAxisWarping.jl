using Test, Statistics, LinearAlgebra
using DynamicAxisWarping, SlidingDistancesBase
using Distances, Plots

@testset "DynamicAxisWarping" begin
    @info "Testing DynamicAxisWarping"

    @testset "GDTW" begin
        @info "Testing GDTW"
        include("test_gdtw.jl")
    end

    @testset "LinearInterpolation" begin
        @info "Testing LinearInterpolation"
        # Test arrays
        x = rand(20, 20, 100)
        x_interp = LinearInterpolation(x)
        @test x_interp isa Function # helps for plotting
        @test x_interp(0) ≈ x[:, :, 1]
        @test x_interp(4/99) ≈ x[:, :, 5]
        @test x_interp(98/99) ≈ x[:, :, end-1]
        @test x_interp( 4.5/99 ) ≈ (x[:, :, 5] + x[:, :, 6])/2
        @test x_interp(1) ≈ x[:, :, end]
        @test x_interp(-1) == zeros(20, 20)
        @test x_interp(2.0) == zeros(20, 20)

        # Test scalars
        x = rand(100)
        x_interp = LinearInterpolation(x)
        @test x_interp(0) == x[1]
        @test x_interp(1) == x[end]
        @test x_interp(4/99) == x[5]
    end

    @testset "Normalizers" begin
        @info "Testing Normalizers"
        a = randn(2,100)
        @test dtwnn(a,a,SqEuclidean(),3,normalizer=IsoZNormalizer).cost < 1e-20
        @test dtwnn(a,a,SqEuclidean(),3,normalizer=ZNormalizer).cost < 1e-20
    end

    @testset "Basic Dynamic Time Warping" begin
        @info "Testing Basic Dynamic Time Warping"
        a = Float64[1, 1, 1, 2, 4, 6, 5, 5, 5, 4, 4, 3, 1, 1, 1]
        b = Float64[1, 1, 2, 4, 6, 6, 6, 5, 4, 4, 4, 3, 3, 3, 1]
        cost, match1, match2 = dtw(a, b)
        @test dtw_cost(a, b, SqEuclidean(), length(a)) == cost
        @test cost == 0
        @test match1 ==
              [1, 2, 3, 4, 5, 6, 6, 6, 7, 8, 9, 10, 10, 11, 12, 12, 12, 13, 14, 15]
        @test match2 == [1, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 14, 15, 15, 15]
        @test evaluate(DTW(10), a, b) == cost

        @test @inferred(soft_dtw_cost(Float64.(a),Float64.(b), γ=0.001)) > -0.01

        a[end] += 2
        cost, match1, match2 = dtw(a, b)
        @test dtw_cost(a, b, SqEuclidean(), length(a)) == cost
        @test cost == 4
        @test evaluate(DTW(10), a, b) == cost
        cost, match1, match2 = dtw(Float64.(a), Float64.(b), transportcost=1.1)
        @test dtw_cost(a, b, SqEuclidean(), length(a)) == cost
        @test evaluate(DTW(radius=10, transportcost=1.1), a, b) == cost

        @test dtw_cost(a, b, SqEuclidean(), 0) ≈ norm(a-b)^2 # Test that radius 0 reduces to SqEuclidean distance
        @test dtw_cost(a, b, Euclidean(), 0) ≈ sum(abs, a-b) # Test that radius 0 reduces to SqEuclidean distance

        @test soft_dtw_cost(Float64.(a),Float64.(b), γ=0.001) ≈ cost rtol = 1e-2
        @test soft_dtw_cost(Float64.(a),Float64.(b), γ=0.001) == SoftDTW(0.001)(a,b)

        a = collect(1.0:10)
        b = a .+ 1
        cost, match1, match2 = dtw(a, b)
        @test dtw_cost(a, b, SqEuclidean(), length(a)) == cost
        @test cost == 2
        @test evaluate(DTW(10), a, b) == cost
        @test soft_dtw_cost(Float64.(a),Float64.(b), γ=0.01) ≈ cost rtol = 1e-2
        @test soft_dtw_cost(Float64.(a),Float64.(b), γ=0.01) == SoftDTW(0.01)(a,b) ≈ SoftDTW(0.01)(big.(a),big.(b))

        a = zeros(6)
        b = 1 .+ a
        cost, match1, match2 = dtw(a, b)
        @test dtw_cost(a, b, SqEuclidean(), length(a)) == cost
        @test cost == length(a)
        @test evaluate(DTW(10), a, b) == cost
        @test soft_dtw_cost(Float64.(a),Float64.(b), γ=0.01) ≈ cost rtol = 1e-2
        @test soft_dtw_cost(Float64.(a),Float64.(b), γ=0.01) == SoftDTW(0.01)(a,b)

        # Verify that a tie prefers diagonal moves
        a = Float64[1, 1, 1]
        b = Float64[1, 1, 1]
        cost, pa, pb = dtw(a, b)
        @test dtw_cost(a, b, SqEuclidean(), length(a)) == cost
        @test cost == 0
        @test pa == [1, 2, 3]
        @test pb == [1, 2, 3]
        @test evaluate(DTW(10), a, b) == cost
        @test soft_dtw_cost(Float64.(a),Float64.(b), γ=0.001) ≈ cost atol = 1e-2
        @test soft_dtw_cost(Float64.(a),Float64.(b), γ=0.001) == SoftDTW(0.001)(a,b)

        # Verify that trackback ends properly if it reaches an edge before reaching [1,1]
        # Also check that trackback prefers diagonal moves
        a = [0, 1, 1, 1]
        b = [0, 0, 1, 1]
        cost, pa, pb = dtw(a, b)
        @test dtw_cost(a, b, SqEuclidean(), length(a)) == cost
        @test cost == 0
        @test pa == [1, 1, 2, 3, 4]
        @test pb == [1, 2, 3, 3, 4]
        @test evaluate(DTW(10), a, b) == cost

        # test the distance api with different distances
        a, b = randn(10), randn(10)
        cost, = dtw(a, b, Euclidean())
        @test evaluate(DTW(10,Euclidean()), a, b) == cost
        @test dtw_cost(a, b, Euclidean(), length(a)) == cost
        cost, = dtw(a, b, Cityblock())
        @test evaluate(DTW(10,Cityblock()), a, b) == cost
        @test dtw_cost(a, b, Cityblock(), length(a)) == cost
        cost, = dtw(a, b, Chebyshev())
        @test evaluate(DTW(10,Chebyshev()), a, b) == cost
        @test dtw_cost(a, b, Chebyshev(), length(a)) == cost


        @test_nowarn dtwplot(a, b)
        @test_nowarn matchplot(a, b)
    end


    @testset "DTW with windows" begin
        @info "Testing DTW with windows"
        # Verify that a tie prefers diagonal moves
        dist = SqEuclidean()
        a = [1, 1, 1]
        b = [1, 1, 1]
        cost, pa, pb = dtw(a, b, dist, [1, 1, 1], [3, 3, 3])
        @test cost == 0
        @test pa == [1, 2, 3]
        @test pb == [1, 2, 3]

        # Verify that trackback ends properly if it reaches an edge before reaching [1,1]
        # Also check that trackback prefers diagonal moves
        a = [0, 1, 1, 1]
        b = [0, 0, 1, 1]
        cost, pa, pb = dtw(a, b, dist, [1, 1, 1, 1], [4, 4, 4, 4])
        @test cost == 0
        @test pa == [1, 1, 2, 3, 4]
        @test pb == [1, 2, 3, 3, 4]

        # First do the windowed test w/o windows
        a = [0, 1, 2, 3, 4, 4, 4, 4]
        b = [0, 0, 1, 2, 2, 2, 3, 4]
        best_pa = [1, 1, 2, 3, 3, 3, 4, 5, 6, 7, 8]
        best_pb = [1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8]
        cost, pa, pb = dtw(a, b)
        @test cost == 0
        @test pa == best_pa
        @test pb == best_pb

        # Wide window, not touching optimal path
        rmin = [1, 1, 1, 2, 3, 4, 5, 6]
        rmax = [4, 6, 7, 8, 8, 8, 8, 8]
        cost, pa, pb = dtw(a, b, dist, rmin, rmax)
        @test cost == 0
        @test pa == best_pa
        @test pb == best_pb

        # Bottom of window is optimal path
        rmin = [1, 3, 4, 7, 8, 8, 8, 8]
        rmax = [4, 6, 7, 8, 8, 8, 8, 8]
        cost, pa, pb = dtw(a, b, dist, rmin, rmax)
        @test cost == 0
        @test pa == best_pa
        @test pb == best_pb

        # Top of window is optimal path
        rmin = [1, 1, 1, 2, 3, 4, 5, 6]
        rmax = [2, 3, 6, 7, 8, 8, 8, 8]
        cost, pa, pb = dtw(a, b, dist, rmin, rmax)
        @test cost == 0
        @test pa == best_pa
        @test pb == best_pb

        # Top and bottom of window are optimal path
        rmin = [1, 3, 4, 7, 8, 8, 8, 8]
        rmax = [2, 3, 6, 7, 8, 8, 8, 8]
        cost, pa, pb = dtw(a, b, dist, rmin, rmax)
        @test cost == 0
        @test pa == best_pa
        @test pb == best_pb

        # Now top of window cuts into optimal path
        rmin = [1, 1, 1, 2, 3, 4, 5, 6]
        rmax = [4, 4, 5, 6, 7, 8, 8, 8]
        cost, pa, pb = dtw(a, b, dist, rmin, rmax)
        @test cost == 2
        @test pa == [1, 1, 2, 3, 3, 4, 5, 6, 7, 8]
        @test pb == [1, 2, 3, 4, 5, 6, 7, 8, 8, 8]

        # Compare windowed and regular on non-square geometry:
        seq1 = collect(1:16)
        seq2 = seq1[1:2:end]
        cost, pa, pb = dtw(seq1, seq2)
        n1, n2 = length(seq1), length(seq2)
        cost2, qa, qb = dtw(seq1, seq2, dist, fill(1, n1), fill(n2, n1))
        @test cost == cost2
        @test pa == qa
        @test pb == qb

        seq1 = [
            1.2051,
            1.7890,
            2.1314,
            2.2423,
            2.2101,
            2.0234,
            1.6898,
            1.3422,
            1.1530,
            1.2197,
            1.5435,
            2.07147,
            2.7359,
            3.469,
            4.2033,
            4.8704,
            5.4056,
            5.75243,
            5.8693,
            5.7362,
            5.35974,
            4.7761,
            4.0482,
            3.2548,
            2.4764,
        ]
        seq2 = [
            1.4927,
            2.1859,
            2.1182,
            1.5186,
            1.1858,
            1.8034,
            3.0968,
            4.5317,
            5.5763,
            5.8038,
            5.0724,
            3.6576,
            2.4823,
        ]

        cost, pa, pb = dtw(seq1, seq2)
        n1, n2 = length(seq1), length(seq2)
        cost2, qa, qb = dtw(seq1, seq2, dist, fill(1, n1), fill(n2, n1))
        @test cost == cost2
        @test pa == qa
        @test pb == qb
    end


    @testset "FastDTW compression" begin
        compress2 = DynamicAxisWarping.compress2
        s = collect(0:2:98)
        s1 = compress2(s)
        s2 = compress2(s1)
        @test s1 == float(collect(1:4:97))
        @test s2 == float(vcat(collect(3:8:91), [97]))

        s = [1]
        s1 = compress2(s)
        @test s1 == [1.0]
    end


    @testset "Window Computations" begin
        @info "Testing Window Computations"
        computewindow = DynamicAxisWarping.computewindow

        # Simplest path (along the diagonal)
        p = collect(1:8)
        rmin, rmax = computewindow(p, p, 1)
        @test rmin == [1, 1, 1, 2, 3, 4, 5, 6]
        @test rmax == [3, 4, 5, 6, 7, 8, 8, 8]

        rmin, rmax = computewindow(p, p, 2)
        @test rmin == [1, 1, 1, 1, 1, 2, 3, 4]
        @test rmax == [5, 6, 7, 8, 8, 8, 8, 8]

        # A warpy path
        pa = [1, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8]
        pb = [1, 2, 3, 3, 3, 4, 4, 5, 6, 7, 8]
        rmin, rmax = computewindow(pa, pb, 1)
        @test pa[end] == length(rmin)
        @test pa[end] == length(rmax)
        @test rmin == [1, 1, 2, 2, 2, 3, 3, 4]
        @test rmax == [4, 4, 4, 5, 5, 6, 8, 8]

        rmin, rmax = computewindow(pa, pb, 2)
        @test pa[end] == length(rmin)
        @test pa[end] == length(rmax)
        @test rmin == [1, 1, 1, 1, 1, 1, 2, 2]
        @test rmax == [5, 5, 6, 6, 7, 8, 8, 8]

        rmin, rmax = computewindow(pa, pb, 20)
        @test pa[end] == length(rmin)
        @test pa[end] == length(rmax)
        @test rmin == fill(1, 8)
        @test rmax == fill(8, 8)

        # Extreme path: follows left then upper edge
        pa = [1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8]
        pb = [1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8]
        rmin, rmax = computewindow(pa, pb, 1)
        @test pa[end] == length(rmin)
        @test pa[end] == length(rmax)
        @test rmin == [1, 1, 7, 7, 7, 7, 7, 7]
        @test rmax == [8, 8, 8, 8, 8, 8, 8, 8]

        rmin, rmax = computewindow(pa, pb, 2)
        @test pa[end] == length(rmin)
        @test pa[end] == length(rmax)
        @test rmin == [1, 1, 1, 6, 6, 6, 6, 6]
        @test rmax == [8, 8, 8, 8, 8, 8, 8, 8]


        # More columns than rows
        pa = [1, 2, 3, 4, 5, 6, 7, 8]
        pb = [1, 2, 3, 4, 4, 4, 4, 4]
        rmin, rmax = computewindow(pa, pb, 1)
        @test rmin == [1, 1, 1, 2, 3, 3, 3, 3]
        @test rmax == [3, 4, 4, 4, 4, 4, 4, 4]

        rmin, rmax = computewindow(pa, pb, 2)
        @test rmin == [1, 1, 1, 1, 1, 2, 2, 2]
        @test rmax == fill(4, 8)

        rmin, rmax = computewindow(pa, pb, 3)
        @test rmin == fill(1, 8)
        @test rmax == fill(4, 8)

        rmin, rmax = computewindow(pa, pb, 4)
        @test rmin == fill(1, 8)
        @test rmax == fill(4, 8)

        rmin, rmax = computewindow(pa, pb, 47)
        @test rmin == fill(1, 8)
        @test rmax == fill(4, 8)

        # More rows than columns
        pa, pb = pb, pa
        rmin, rmax = computewindow(pa, pb, 1)
        @test rmin == [1, 1, 1, 2]
        @test rmax == [3, 4, 8, 8]

        rmin, rmax = computewindow(pa, pb, 2)
        @test rmin == fill(1, 4)
        @test rmax == [5, 8, 8, 8]

        rmin, rmax = computewindow(pa, pb, 3)
        @test rmin == fill(1, 4)
        @test rmax == fill(8, 4)

        rmin, rmax = computewindow(pa, pb, 4)
        @test rmin == fill(1, 4)
        @test rmax == fill(8, 4)

        rmin, rmax = computewindow(pa, pb, 47)
        @test rmin == fill(1, 4)
        @test rmax == fill(8, 4)
    end


    @testset "DTW and FastDTW agreement" begin
        t = collect(1:1600)
        pktimes = [100, 300, 1000, 1300]
        x = 1 * exp.(-0.5 * ((t .- pktimes[1]) / 100) .^ 2)
        x += 2 * exp.(-0.5 * ((t .- pktimes[2]) / 150) .^ 2)
        x += 3 * exp.(-0.5 * ((t .- pktimes[3]) / 250) .^ 2)
        x += 4 * exp.(-0.5 * ((t .- pktimes[4]) / 250) .^ 2)
        y = x[1:2:end]
        cost, px, py = dtw(x, y)
        cost1, qx, qy = fastdtw(x, y, SqEuclidean(), 15)
        @test isapprox(cost, cost1)
        @test px == qx
        @test py == qy
    end

    @testset "DBA" begin
        x = [1.0, 2.0, 2.0, 3.0, 3.0, 4.0]
        y = [1.0, 3.0, 4.0]
        z = [1.0, 2.0, 2.0, 4.0]
        avg, _ = dba([x, y, z], DTW(5), init_center = z)
        @test avg == [1.0, 1.75, 2.75, 4.0]
    end


    @testset "DTW NN" begin
        @info "Testing DTW NN"
        a = randn(Float32, 100)
        b = randn(Float32, 10000)

        function naive(a, b, r=7)
            dists = map(1:length(b)-length(a)+1) do i
                dtw_cost(a, @view(b[i:i+length(a)-1]), SqEuclidean(), r)
            end
        end

        @inferred dtwnn(a, b, SqEuclidean(), 7)

        res = dtwnn(a, b, SqEuclidean(), 7)
        m = findmin(naive(a, b))
        @test m[1] ≈ res.cost
        @test m[2] == res.loc

        res = dtwnn(a, b, SqEuclidean(), 7, saveall=true)
        resn = naive(a, b)
        m = findmin(resn)
        @test m[1] ≈ res.cost
        @test m[2] == res.loc
        @test res.dists ≈ resn


        a = randn(Float64, 100)
        b = randn(Float64, 10000)

        function naive_norm(a, b)
            an = normalize(ZNormalizer, a)
            dists = map(1:length(b)-length(a)+1) do i
                bn = normalize(ZNormalizer, b[i:i+length(a)-1])
                @test mean(bn) ≈ 0 atol = 10eps(eltype(a))
                @test std(bn, corrected=false, mean=0) ≈ 1 atol = sqrt(eps(eltype(a)))
                dtw_cost(an, bn, SqEuclidean(), 7)
            end
        end

        @inferred dtwnn(a, b, SqEuclidean(), 7, normalizer=Val(ZNormalizer))
        res = dtwnn(a, b, SqEuclidean(), 7, normalizer=Val(ZNormalizer), saveall=true)
        resn = naive_norm(a, b)
        m = findmin(resn)
        @test m[1] ≈ res.cost
        @test m[2] == res.loc
        @test res.dists ≈ resn


        a = sin.(0.1 .* (1:100))    .+ 0.05 .* randn.()
        b = sin.(0.1 .* (1:10_000)) .+ 0.05 .* randn.()

        res = dtwnn(a, b, SqEuclidean(), 2, prune_endpoints = false, prune_envelope = true)
        @test res.prunestats.prune_env > 0
        resn = naive(a, b, 2)
        m = findmin(resn)
        @test m[2] == res.loc
        @test m[1] ≈ res.cost

        res = dtwnn(a, b, SqEuclidean(), 2, prune_endpoints = true, prune_envelope = true)
        @test res.prunestats.prune_env > 0
        @test res.prunestats.prune_end > 0
        resn = naive(a, b, 2)
        m = findmin(resn)
        @test m[2] == res.loc
        @test m[1] ≈ res.cost

        @test SlidingDistancesBase.value(res) == res.cost
        @test SlidingDistancesBase.location(res) == res.loc
        @test SlidingDistancesBase.payload(res) == res.dists
        @test SlidingDistancesBase.target(res) == a

        res = dtwnn(a, b, SqEuclidean(), 2, prune_endpoints = true, prune_envelope = true, saveall=true)
        rr = [res, res]
        @test_nowarn 2res
        @test minimum(rr) == SlidingDistancesBase.value(res)
        @test findmin(rr) == (res.cost, 1)
        @test maximum(rr) == SlidingDistancesBase.value(res)
        @test findmax(rr) == (res.cost, 2)


    end


    @testset "DBA clust" begin
        allsame(x) = all(==(x[1]), x)
        data = [randn(100) .+ 2(i ÷ 5) for i = 0:19]
        nclust = 4
        result = dbaclust(
            data,
            nclust,
            DTW(10);
            n_init = 20,
            iterations = 10,
        )
        inds = 1:5
        @test all([allsame(result.clustids[inds .+ 5i]) for i in 0:3])


        result = [dbaclust(
            data,
            nclust,
            FastDTW(10);
            n_init = 20,
            iterations = 10,
        ) for _ in 1:2]
        inds = 1:5
        @test any(all([allsame(result.clustids[inds .+ 5i]) for i in 0:3]) for result in result)
    end


    @testset "sparse_distmat" begin
        @info "Testing sparse_distmat"
        using Distances
        y = [sin.(0.1 .* (1:100) .+ 2pi*rand()) .+ 0.1randn(100) for _ in 1:200]
        dists,inds = DynamicAxisWarping.sparse_distmat(y, 4, SqEuclidean(), 5)
        D = [dtw_cost(y1,y2,SqEuclidean(), 5) for y1 in y, y2 in y]
        for i = 1:200
            @test partialsort(D[:,i], 2:5) == dists[i]
            @test partialsortperm(D[:,i], 2:5) == inds[i]
        end

    end

    @testset "matrix_profile" begin
        @info "Testing matrix_profile"
        include("test_matrixprofile.jl")
    end

    @testset "distance profile" begin
        @info "Testing distance profile"
        Q = randn(10)
        T = randn(100)
        d = DTW(5)
        D = distance_profile(d, Q, T)
        @test length(D) == 100-10+1
        foreach(i->@test(D[i] ≈ d(Q, getwindow(T,10,i))), 1:91)
    end


    @testset "datasets" begin
        @info "Testing datasets"

        @testset "gaussian" begin
            @info "Testing gaussian"
            data, labels = DynamicAxisWarping.Datasets.fakedata_gaussian()
        end

        @testset "ucr" begin
            @info "Testing ucr"

        end

    end

end


# using DynamicAxisWarping, Distances
# a = randn(1000)
# b = randn(1000)
# @btime dtw($a, $b)
#
# # 10.526 ms (9 allocations: 7.68 MiB)
# # add @avx
# # 2.245 ms (9 allocations: 7.68 MiB)
#
# imi,ima = radiuslimits(10,1000,1000)
# @btime dtw($a, $b, SqEuclidean(), $imi, $ima)
# # 196.412 μs (18 allocations: 193.91 KiB)
#
#
#
# s1,s2 = zeros(21),zeros(21)
# @btime dtw_cost($a, $b, SqEuclidean(), 10, cost=$s1, cost_prev=$s2)
#
# dtw(a, b, SqEuclidean(), imi, ima)
# dtw_cost(a, b, SqEuclidean(), 10)


# using DynamicAxisWarping, Distances
# a = randn(Float32, 1000)
# b = randn(Float32, 100000)
#
# function naive(a,b)
#     dists = map(1:length(b)-length(a)) do i
#         dtw_cost(a,@view(b[i:i+length(a)-1]),SqEuclidean(), 7)
#     end
# end
#
# w = DTWWorkspace(a,SqEuclidean(),7)
# @inferred dtwnn(w,b)
#
# res = dtwnn(w,b)
# m = findmin(naive(a,b))
# @test m[1] ≈ res.cost
# @test m[2] == res.loc
#
#
# @btime naive($a,$b)
# @btime dtwnn($w,$b)
#
#
# a = sin.(0.1 .* (1:100)) .+ 0.1 .* randn.()
# b = sin.(0.1 .* (1:100000)) .+ 0.1 .* randn.()
# w = DTWWorkspace(a,SqEuclidean(),7)
# res = dtwnn(w,b)
#
# plot([a b[eachindex(a) .+ (res.loc-1)]])
#
#
# res = dtwnn(w,b)
# m = findmin(naive(a,b))
# @test m[1] ≈ res.cost
# @test m[2] == res.loc


# a      = sin.((1:1000))     .+ 0.05 .* randn.()
# b      = sin.((1:1000_000)) .+ 0.05 .* randn.()
# @time dtwnn(a, b, SqEuclidean(), 7, normalizer=ZNormalizer, saveall=false)
# @btime dtwnn($a, $b, SqEuclidean(), 5, normalizer=Nothing, saveall=false)
# @btime naive_norm($a, $b)






# # Bench prune_envelope
# a = sin.(0.1 .* (1:100))    .+ 0.1 .* randn.()
# b = sin.(0.1 .* (1:1000_000)) .+ 0.1 .* randn.()
# @btime dtwnn($a, $b, SqEuclidean(), 5, prune_endpoints = false, prune_envelope = false)
# @btime dtwnn($a, $b, SqEuclidean(), 5, prune_endpoints = false, prune_envelope = true)
# @btime dtwnn($a, $b, SqEuclidean(), 5, prune_endpoints = true, prune_envelope = false)
# @btime dtwnn($a, $b, SqEuclidean(), 5, prune_endpoints = true, prune_envelope = true)
#
#
# a = sin.(0.1f0 .* (1:100))    .+ 0.1f0 .* randn.(Float32)
# b = sin.(0.1f0 .* (1:1000_000)) .+ 0.1f0 .* randn.(Float32)
# @btime dtwnn($a, $b, SqEuclidean(), 5, prune_endpoints = true, prune_envelope = true, normalizer=Val(ZNormalizer))
