using Base.Test
using TimeWarp

@testset "Sequences" begin

    # 10-dimensional time series
    X = randn(10,3)
    Sx = Sequence(X)

    @test size(Sx) == (3,)
    @test length(Sx) == 3

    for i = 1:length(Sx)
        @test Sx[i] == X[:,i]
    end

    # (10x9)-dimensional time series
    Y = randn(10,9,3)
    Sy = Sequence(Y)

    @test size(Sy) == (3,)
    @test length(Sy) == 3

    for i = 1:length(Sy)
        @test Sy[i] == Y[:,:,i]
    end

end

@testset "Basic Dynamic Time Warping" begin
    a=[1,1,1,2,4,6,5,5,5,4,4,3,1,1,1]
    b=[1,1,2,4,6,6,6,5,4,4,4,3,3,3,1]
    cost, match1, match2 = dtw(a,b)
    @test cost==0
    @test match1==[1,2,3,4,5,6,6,6,7,8,9,10,10,11,12,12,12,13,14,15]
    @test match2==[1,1,2,3,4,5,6,7,8,8,8, 9,10,11,12,13,14,15,15,15]
    @test evaluate(DTWDistance(),a,b) == cost

    a[end] += 2
    cost, match1, match2 = dtw(a,b)
    @test cost==4
    @test evaluate(DTWDistance(),a,b) == cost

    a=collect(1:10)
    b=a+1
    cost, match1, match2 = dtw(a,b)
    @test cost==2
    @test evaluate(DTWDistance(),a,b) == cost

    a=zeros(Int,6)
    b=1+a
    cost, match1, match2 = dtw(a,b)
    @test cost==length(a)
    @test evaluate(DTWDistance(),a,b) == cost

    # Verify that a tie prefers diagonal moves
    a=[1,1,1]
    b=[1,1,1]
    cost, pa, pb = dtw(a,b)
    @test cost==0
    @test pa==[1,2,3]
    @test pb==[1,2,3]
    @test evaluate(DTWDistance(),a,b) == cost

    # Verify that trackback ends properly if it reaches an edge before reaching [1,1]
    # Also check that trackback prefers diagonal moves
    a=[0,1,1,1]
    b=[0,0,1,1]
    cost, pa, pb = dtw(a,b)
    @test cost==0
    @test pa==[1,1,2,3,4]
    @test pb==[1,2,3,3,4]
    @test evaluate(DTWDistance(),a,b) == cost

    # test the distance api with different distances
    a,b = randn(10),randn(10)
    cost, = dtw(a,b,Euclidean())
    @test evaluate(DTWDistance(Euclidean()),a,b) == cost
    cost, = dtw(a,b,Cityblock())
    @test evaluate(DTWDistance(Cityblock()),a,b) == cost
    cost, = dtw(a,b,Chebyshev())
    @test evaluate(DTWDistance(Chebyshev()),a,b) == cost
end


@testset "DTW with windows" begin
    # Verify that a tie prefers diagonal moves
    a=[1,1,1]
    b=[1,1,1]
    cost, pa, pb = dtw(a,b,[1,1,1],[3,3,3])
    @test cost==0
    @test pa==[1,2,3]
    @test pb==[1,2,3]

    # Verify that trackback ends properly if it reaches an edge before reaching [1,1]
    # Also check that trackback prefers diagonal moves
    a=[0,1,1,1]
    b=[0,0,1,1]
    cost, pa, pb = dtw(a,b,[1,1,1,1],[4,4,4,4])
    @test cost==0
    @test pa==[1,1,2,3,4]
    @test pb==[1,2,3,3,4]

    # First do the windowed test w/o windows
    a=[0,1,2,3,4,4,4,4]
    b=[0,0,1,2,2,2,3,4]
    best_pa = [1,1,2,3,3,3,4,5,6,7,8]
    best_pb = [1,2,3,4,5,6,7,8,8,8,8]
    cost, pa, pb = dtw(a,b)
    @test cost == 0
    @test pa == best_pa
    @test pb == best_pb

    # Wide window, not touching optimal path
    rmin = [1,1,1,2,3,4,5,6]
    rmax = [4,6,7,8,8,8,8,8]
    cost, pa, pb = dtw(a,b,rmin,rmax)
    @test cost == 0
    @test pa == best_pa
    @test pb == best_pb

    # Bottom of window is optimal path
    rmin = [1,3,4,7,8,8,8,8]
    rmax = [4,6,7,8,8,8,8,8]
    cost, pa, pb = dtw(a,b,rmin,rmax)
    @test cost == 0
    @test pa == best_pa
    @test pb == best_pb

    # Top of window is optimal path
    rmin = [1,1,1,2,3,4,5,6]
    rmax = [2,3,6,7,8,8,8,8]
    cost, pa, pb = dtw(a,b,rmin,rmax)
    @test cost == 0
    @test pa == best_pa
    @test pb == best_pb

    # Top and bottom of window are optimal path
    rmin = [1,3,4,7,8,8,8,8]
    rmax = [2,3,6,7,8,8,8,8]
    cost, pa, pb = dtw(a,b,rmin,rmax)
    @test cost == 0
    @test pa == best_pa
    @test pb == best_pb

    # Now top of window cuts into optimal path
    rmin = [1,1,1,2,3,4,5,6]
    rmax = [4,4,5,6,7,8,8,8]
    cost, pa, pb = dtw(a,b,rmin,rmax)
    @test cost == 2
    @test pa == [1,1,2,3,3,4,5,6,7,8]
    @test pb == [1,2,3,4,5,6,7,8,8,8]

    # Compare windowed and regular on non-square geometry:
    seq1=collect(1:16)
    seq2=seq1[1:2:end]
    cost,pa,pb = dtw(seq1,seq2)
    n1,n2 = length(seq1), length(seq2)
    cost2,qa,qb = dtw(seq1, seq2, fill(1,n1), fill(n2,n1))
    @test cost==cost2
    @test pa==qa
    @test pb==qb

    seq1=[1.2051,1.7890,2.1314,2.2423,2.2101,2.0234,1.6898,1.3422,1.1530,1.2197,1.5435,2.07147,2.7359,3.469,4.2033,4.8704,5.4056,5.75243,5.8693,5.7362,5.35974,4.7761,4.0482,3.2548,2.4764]
    seq2=[1.4927,2.1859,2.1182,1.5186, 1.1858,1.8034,3.0968,4.5317,5.5763,5.8038,5.0724,3.6576,2.4823]

    cost,pa,pb = dtw(seq1,seq2)
    n1,n2 = length(seq1), length(seq2)
    cost2,qa,qb = dtw(seq1, seq2, fill(1,n1), fill(n2,n1))
    @test cost==cost2
    @test pa==qa
    @test pb==qb
end


@testset "FastDTW compression" begin
    compress2 = TimeWarp.compress2
    s=collect(0:2:98)
    s1 = compress2(s)
    s2 = compress2(s1)
    @test s1==float(collect(1:4:97))
    @test s2==float(vcat(collect(3:8:91),[97]))

    s = [1]
    s1 = compress2(s)
    @test s1==[1.0]
end


@testset "Window Computations" begin
    computewindow = TimeWarp.computewindow

    # Simplest path (along the diagonal)
    p=collect(1:8)
    rmin,rmax = computewindow(p,p,1)
    @test rmin==[1,1,1,2,3,4,5,6]
    @test rmax==[3,4,5,6,7,8,8,8]

    rmin,rmax = computewindow(p,p,2)
    @test rmin==[1,1,1,1,1,2,3,4]
    @test rmax==[5,6,7,8,8,8,8,8]

    # A warpy path
    pa=[1,1,2,3,4,5,6,7,8,8,8]
    pb=[1,2,3,3,3,4,4,5,6,7,8]
    rmin,rmax = computewindow(pa, pb, 1)
    @test pa[end]==length(rmin)
    @test pa[end]==length(rmax)
    @test rmin==[1,1,2,2,2,3,3,4]
    @test rmax==[4,4,4,5,5,6,8,8]

    rmin,rmax = computewindow(pa, pb, 2)
    @test pa[end]==length(rmin)
    @test pa[end]==length(rmax)
    @test rmin==[1,1,1,1,1,1,2,2]
    @test rmax==[5,5,6,6,7,8,8,8]

    rmin,rmax = computewindow(pa, pb, 20)
    @test pa[end]==length(rmin)
    @test pa[end]==length(rmax)
    @test rmin==fill(1,8)
    @test rmax==fill(8,8)

    # Extreme path: follows left then upper edge
    pa=[1,1,1,1,1,1,1,1,2,3,4,5,6,7,8]
    pb=[1,2,3,4,5,6,7,8,8,8,8,8,8,8,8]
    rmin,rmax = computewindow(pa, pb, 1)
    @test pa[end]==length(rmin)
    @test pa[end]==length(rmax)
    @test rmin==[1,1,7,7,7,7,7,7]
    @test rmax==[8,8,8,8,8,8,8,8]

    rmin,rmax = computewindow(pa, pb, 2)
    @test pa[end]==length(rmin)
    @test pa[end]==length(rmax)
    @test rmin==[1,1,1,6,6,6,6,6]
    @test rmax==[8,8,8,8,8,8,8,8]


    # More columns than rows
    pa=[1,2,3,4,5,6,7,8]
    pb=[1,2,3,4,4,4,4,4]
    rmin,rmax = computewindow(pa, pb, 1)
    @test rmin==[1,1,1,2,3,3,3,3]
    @test rmax==[3,4,4,4,4,4,4,4]

    rmin,rmax = computewindow(pa, pb, 2)
    @test rmin==[1,1,1,1,1,2,2,2]
    @test rmax==fill(4,8)

    rmin,rmax = computewindow(pa, pb, 3)
    @test rmin==fill(1,8)
    @test rmax==fill(4,8)

    rmin,rmax = computewindow(pa, pb, 4)
    @test rmin==fill(1,8)
    @test rmax==fill(4,8)

    rmin,rmax = computewindow(pa, pb, 47)
    @test rmin==fill(1,8)
    @test rmax==fill(4,8)

    # More rows than columns
    pa,pb = pb,pa
    rmin,rmax = computewindow(pa, pb, 1)
    @test rmin==[1,1,1,2]
    @test rmax==[3,4,8,8]

    rmin,rmax = computewindow(pa, pb, 2)
    @test rmin==fill(1,4)
    @test rmax==[5,8,8,8]

    rmin,rmax = computewindow(pa, pb, 3)
    @test rmin==fill(1,4)
    @test rmax==fill(8,4)

    rmin,rmax = computewindow(pa, pb, 4)
    @test rmin==fill(1,4)
    @test rmax==fill(8,4)

    rmin,rmax = computewindow(pa, pb, 47)
    @test rmin==fill(1,4)
    @test rmax==fill(8,4)
end


@testset "DTW and FastDTW agreement" begin
    t=collect(1:1600)
    pktimes=[100,300,1000,1300]
    x = 1*exp.(-0.5*((t-pktimes[1])/100).^2);
    x+= 2*exp.(-0.5*((t-pktimes[2])/150).^2);
    x+= 3*exp.(-0.5*((t-pktimes[3])/250).^2);
    x+= 4*exp.(-0.5*((t-pktimes[4])/250).^2);
    y = x[1:2:end];
    cost,px,py = dtw(x,y)
    cost1,qx,qy = fastdtw(x, y, 15)
    @test isapprox(cost,cost1)
    @test px == qx
    @test py == qy
end
