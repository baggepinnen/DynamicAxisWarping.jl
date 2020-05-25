
using LinearAlgebra: diagind
using SlidingDistancesBase: getwindow
using MatrixProfile
m = 10
r = 5
y = sin.(0.1 .* (1:1000) .+ 2pi*rand()) .+ 0.1randn(1000)
D = [dtw_cost(getwindow(y, m, i), getwindow(y, m, j), SqEuclidean(), r) for i in 1:length(y)-m+1, j in 1:length(y)-m+1]

for i = -r:r
    D[diagind(D,i)] .= Inf
end
profile = matrix_profile(y,m,DTW(r))
for i = 1:10:size(D,2)
    @test findmin(D[:,i]) == (profile.P[i], profile.I[i])
end

# i = 1
# findmin(D[:,i])
# (profile.P[i], profile.I[i])
